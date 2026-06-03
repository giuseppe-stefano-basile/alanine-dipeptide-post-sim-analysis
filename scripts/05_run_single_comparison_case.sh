#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/production_comparison_cases.json"
CASE_NAME=""
BOOTSTRAP_REPS="1000"
HV_FRAME_STRIDE="1"
WORKSPACE_ROOT="$REPO_DIR/workspaces"
RESULTS_ROOT="$REPO_DIR/results"
SEARCH_ROOTS=()
BOOTSTRAP_SET=0
HV_SET=0

usage() {
  cat <<EOF
Usage: $0 <case_name> --search-root DIR [options]

Required:
  <case_name> or --case NAME
  --search-root DIR             Directory to search for trajectories/logs. Repeatable.

Options:
  --cases-json FILE             Case config file. Default: $CASES_JSON
  --bootstrap-reps N            Bootstrap replicates. Default: $BOOTSTRAP_REPS
  --high-value-frame-stride N   Frame stride for high-value observables. Default: $HV_FRAME_STRIDE
  --workspace-root DIR          Where generated split inputs are written. Default: $WORKSPACE_ROOT
  --results-root DIR            Where analysis outputs are written. Default: $RESULTS_ROOT
  -h, --help                    Show this help.

Legacy positional form is still accepted:
  $0 <case_name> <bootstrap_reps> <high_value_frame_stride> --search-root DIR
EOF
}

display_path() {
  python3 - "$REPO_DIR" "$1" <<'PY'
import sys
from pathlib import Path

base = Path(sys.argv[1]).resolve()
path = Path(sys.argv[2]).expanduser().absolute()
try:
    print(path.relative_to(base))
except ValueError:
    print(path)
PY
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --case)
      CASE_NAME="${2:?missing value for --case}"
      shift 2
      ;;
    --search-root)
      SEARCH_ROOTS+=("${2:?missing value for --search-root}")
      shift 2
      ;;
    --cases-json)
      CASES_JSON="${2:?missing value for --cases-json}"
      shift 2
      ;;
    --bootstrap-reps)
      BOOTSTRAP_REPS="${2:?missing value for --bootstrap-reps}"
      BOOTSTRAP_SET=1
      shift 2
      ;;
    --high-value-frame-stride|--frame-stride)
      HV_FRAME_STRIDE="${2:?missing value for --high-value-frame-stride}"
      HV_SET=1
      shift 2
      ;;
    --workspace-root)
      WORKSPACE_ROOT="${2:?missing value for --workspace-root}"
      shift 2
      ;;
    --results-root)
      RESULTS_ROOT="${2:?missing value for --results-root}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -z "$CASE_NAME" ]]; then
        CASE_NAME="$1"
      elif [[ "$BOOTSTRAP_SET" -eq 0 ]]; then
        BOOTSTRAP_REPS="$1"
        BOOTSTRAP_SET=1
      elif [[ "$HV_SET" -eq 0 ]]; then
        HV_FRAME_STRIDE="$1"
        HV_SET=1
      else
        echo "ERROR: unknown argument: $1" >&2
        usage >&2
        exit 2
      fi
      shift
      ;;
  esac
done

if [[ -z "$CASE_NAME" ]]; then
  echo "ERROR: choose a case name" >&2
  usage >&2
  exit 2
fi
if [[ ${#SEARCH_ROOTS[@]} -eq 0 ]]; then
  echo "ERROR: choose at least one directory with --search-root" >&2
  usage >&2
  exit 2
fi
if [[ ! -f "$CASES_JSON" ]]; then
  echo "ERROR: cases file not found: $CASES_JSON" >&2
  exit 2
fi

RESOLVE_CMD=(python3 "$SCRIPT_DIR/resolve_production_inputs.py" --cases-json "$CASES_JSON" --case "$CASE_NAME" --format tsv)
for root in "${SEARCH_ROOTS[@]}"; do
  RESOLVE_CMD+=(--search-root "$root")
done

if ! RESOLVED_OUTPUT="$("${RESOLVE_CMD[@]}")"; then
  exit 2
fi
mapfile -t RESOLVED_CASES <<< "$RESOLVED_OUTPUT"
if [[ "${#RESOLVED_CASES[@]}" -ne 1 || -z "${RESOLVED_CASES[0]}" ]]; then
  echo "ERROR: expected exactly one resolved case for '$CASE_NAME'" >&2
  exit 2
fi
IFS=$'\t' read -r CASE_NAME CASE_DESC NPBC_DUMP PBC_DUMP NPBC_LOG PBC_LOG <<< "${RESOLVED_CASES[0]}"

if [[ ! -f "$NPBC_DUMP" ]]; then
  echo "ERROR: NPBC dump missing: $NPBC_DUMP"
  exit 2
fi
if [[ ! -f "$PBC_DUMP" ]]; then
  echo "ERROR: PBC dump missing: $PBC_DUMP"
  exit 2
fi

WORKSPACE_DIR="$(python3 -c 'import sys; from pathlib import Path; print(Path(sys.argv[1]).expanduser().resolve())' "$WORKSPACE_ROOT")/$CASE_NAME"
RESULT_DIR="$(python3 -c 'import sys; from pathlib import Path; print(Path(sys.argv[1]).expanduser().resolve())' "$RESULTS_ROOT")/$CASE_NAME"
MAIN_OUT="$RESULT_DIR/01_main_compare"
HV_OUT="$RESULT_DIR/02_high_value_observables"
SIG_OUT="$RESULT_DIR/03_significance"
WORKSPACE_DISPLAY="$(display_path "$WORKSPACE_DIR")"
MAIN_OUT_DISPLAY="$(display_path "$MAIN_OUT")"
HV_OUT_DISPLAY="$(display_path "$HV_OUT")"
SIG_OUT_DISPLAY="$(display_path "$SIG_OUT")"
NPBC_DUMP_DISPLAY="$(display_path "$NPBC_DUMP")"
PBC_DUMP_DISPLAY="$(display_path "$PBC_DUMP")"

mkdir -p "$WORKSPACE_DIR" "$RESULT_DIR"

echo "============================================================"
echo "Running case: $CASE_NAME"
echo "Description: $CASE_DESC"
echo "NPBC prod dump: $NPBC_DUMP"
echo "PBC prod dump:  $PBC_DUMP"
if [[ -n "${NPBC_LOG:-}" ]]; then
  echo "NPBC prod log:  $NPBC_LOG"
fi
if [[ -n "${PBC_LOG:-}" ]]; then
  echo "PBC prod log:   $PBC_LOG"
fi

echo "[1/5] Split production dumps into pseudo_eq + pseudo_prod"
python3 "$SCRIPT_DIR/01_prepare_production_chunks_for_legacy_tools.py" \
  --npbc-prod "$NPBC_DUMP" \
  --pbc-prod "$PBC_DUMP" \
  --workspace "$WORKSPACE_DIR" \
  --eq-frames 1001

echo "[2/5] Run comprehensive main comparison"
python3 "$SCRIPT_DIR/02_run_fes_basin_structural_analysis.py" \
  --workspace "$WORKSPACE_DIR" \
  --outdir "$MAIN_OUT" \
  --bootstrap-reps "$BOOTSTRAP_REPS"

echo "[3/5] Run high-value observables"
python3 "$SCRIPT_DIR/03_run_hydration_hbond_coordination_analysis.py" \
  --workspace "$WORKSPACE_DIR" \
  --outdir "$HV_OUT" \
  --frame-stride "$HV_FRAME_STRIDE"

echo "[4/5] Assess FES significance"
python3 "$SCRIPT_DIR/04_summarize_fes_statistical_significance.py" \
  --main-compare-dir "$MAIN_OUT" \
  --outdir "$SIG_OUT"

echo "[5/5] Diagnose rare-basin transitions"
python3 "$SCRIPT_DIR/06_diagnose_rare_basin_visits.py" \
  --main-compare-dir "$MAIN_OUT" \
  --outdir "$SIG_OUT"

cat > "$RESULT_DIR/COMPARISON_CASE_SUMMARY.md" <<MD
# Case Summary: $CASE_NAME

- Description: $CASE_DESC
- Workspace: $WORKSPACE_DISPLAY
- Main comparison output: $MAIN_OUT_DISPLAY
- High-value output: $HV_OUT_DISPLAY
- Significance output: $SIG_OUT_DISPLAY
- Rare-basin diagnostics: $SIG_OUT_DISPLAY/rare_basin_diagnostics.md
- NPBC production dump: $NPBC_DUMP_DISPLAY
- PBC production dump: $PBC_DUMP_DISPLAY
- Bootstrap reps: $BOOTSTRAP_REPS
- High-value frame stride: $HV_FRAME_STRIDE
MD

echo "Case complete: $CASE_NAME"
echo "Results: $RESULT_DIR"
echo "============================================================"
