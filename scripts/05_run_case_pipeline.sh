#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <case_name> [bootstrap_reps] [high_value_frame_stride]"
  exit 2
fi

CASE_NAME="$1"
BOOTSTRAP_REPS="${2:-1000}"
HV_FRAME_STRIDE="${3:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/cases.json"

read_case_field() {
  local case_name="$1"
  local field="$2"
  python3 - "$CASES_JSON" "$case_name" "$field" <<'PY'
import json
import sys
from pathlib import Path

cases = json.loads(Path(sys.argv[1]).read_text()).get("cases", [])
name = sys.argv[2]
field = sys.argv[3]
for c in cases:
    if c.get("name") == name:
        val = c.get(field)
        if val is None:
            print("")
        else:
            print(val)
        break
else:
    raise SystemExit(f"case not found: {name}")
PY
}

NPBC_DUMP_REL="$(read_case_field "$CASE_NAME" "npbc_prod_dump")"
PBC_DUMP_REL="$(read_case_field "$CASE_NAME" "pbc_prod_dump")"
CASE_DESC="$(read_case_field "$CASE_NAME" "description")"

NPBC_DUMP="$REPO_DIR/$NPBC_DUMP_REL"
PBC_DUMP="$REPO_DIR/$PBC_DUMP_REL"

if [[ ! -f "$NPBC_DUMP" ]]; then
  echo "ERROR: NPBC dump missing: $NPBC_DUMP"
  exit 2
fi
if [[ ! -f "$PBC_DUMP" ]]; then
  echo "ERROR: PBC dump missing: $PBC_DUMP"
  exit 2
fi

WORKSPACE_DIR="$REPO_DIR/workspaces/$CASE_NAME"
RESULT_DIR="$REPO_DIR/results/$CASE_NAME"
MAIN_OUT="$RESULT_DIR/01_main_compare"
HV_OUT="$RESULT_DIR/02_high_value_observables"
SIG_OUT="$RESULT_DIR/03_significance"

mkdir -p "$WORKSPACE_DIR" "$RESULT_DIR"

echo "============================================================"
echo "Running case: $CASE_NAME"
echo "Description: $CASE_DESC"
echo "NPBC prod dump: $NPBC_DUMP"
echo "PBC prod dump:  $PBC_DUMP"

echo "[1/4] Split production dumps into pseudo_eq + pseudo_prod"
python3 "$SCRIPT_DIR/01_split_production_dump.py" \
  --npbc-prod "$NPBC_DUMP" \
  --pbc-prod "$PBC_DUMP" \
  --workspace "$WORKSPACE_DIR" \
  --eq-frames 1001

echo "[2/4] Run comprehensive main comparison"
python3 "$SCRIPT_DIR/02_run_main_compare.py" \
  --workspace "$WORKSPACE_DIR" \
  --outdir "$MAIN_OUT" \
  --bootstrap-reps "$BOOTSTRAP_REPS"

echo "[3/4] Run high-value observables"
python3 "$SCRIPT_DIR/03_run_high_value_observables.py" \
  --workspace "$WORKSPACE_DIR" \
  --outdir "$HV_OUT" \
  --frame-stride "$HV_FRAME_STRIDE"

echo "[4/4] Assess FES significance"
python3 "$SCRIPT_DIR/04_assess_fes_significance.py" \
  --main-compare-dir "$MAIN_OUT" \
  --outdir "$SIG_OUT"

echo "[5/5] Diagnose rare-basin transitions"
python3 "$SCRIPT_DIR/06_rare_basin_diagnostics.py" \
  --main-compare-dir "$MAIN_OUT" \
  --outdir "$SIG_OUT"

cat > "$RESULT_DIR/CASE_SUMMARY.md" <<MD
# Case Summary: $CASE_NAME

- Description: $CASE_DESC
- Workspace: $WORKSPACE_DIR
- Main comparison output: $MAIN_OUT
- High-value output: $HV_OUT
- Significance output: $SIG_OUT
- Rare-basin diagnostics: $SIG_OUT/rare_basin_diagnostics.md
- Bootstrap reps: $BOOTSTRAP_REPS
- High-value frame stride: $HV_FRAME_STRIDE
MD

echo "Case complete: $CASE_NAME"
echo "Results: $RESULT_DIR"
echo "============================================================"
