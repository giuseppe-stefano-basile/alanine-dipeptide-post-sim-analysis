#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/production_comparison_cases.json"
BOOTSTRAP_REPS="1000"
HV_FRAME_STRIDE="1"
WORKSPACE_ROOT="$REPO_DIR/workspaces"
RESULTS_ROOT="$REPO_DIR/results"
SEARCH_ROOTS=()
CASES=()
BOOTSTRAP_SET=0
HV_SET=0

usage() {
  cat <<EOF
Usage: $0 --search-root DIR [options]

Required:
  --search-root DIR             Directory to search for trajectories/logs. Repeatable.

Options:
  --case NAME                   Run one configured case. Repeatable. Default: all cases.
  --cases-json FILE             Case config file. Default: $CASES_JSON
  --bootstrap-reps N            Bootstrap replicates. Default: $BOOTSTRAP_REPS
  --high-value-frame-stride N   Frame stride for high-value observables. Default: $HV_FRAME_STRIDE
  --workspace-root DIR          Where generated split inputs are written. Default: $WORKSPACE_ROOT
  --results-root DIR            Where analysis outputs are written. Default: $RESULTS_ROOT
  -h, --help                    Show this help.

Legacy positional form is still accepted:
  $0 <bootstrap_reps> <high_value_frame_stride> --search-root DIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --search-root)
      SEARCH_ROOTS+=("${2:?missing value for --search-root}")
      shift 2
      ;;
    --case)
      CASES+=("${2:?missing value for --case}")
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
      if [[ "$BOOTSTRAP_SET" -eq 0 ]]; then
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

if [[ ${#SEARCH_ROOTS[@]} -eq 0 ]]; then
  echo "ERROR: choose at least one directory with --search-root" >&2
  usage >&2
  exit 2
fi

RESOLVE_CMD=(python3 "$SCRIPT_DIR/resolve_production_inputs.py" --cases-json "$CASES_JSON" --format names)
for root in "${SEARCH_ROOTS[@]}"; do
  RESOLVE_CMD+=(--search-root "$root")
done
for case_name in "${CASES[@]}"; do
  RESOLVE_CMD+=(--case "$case_name")
done

if ! RESOLVED_OUTPUT="$("${RESOLVE_CMD[@]}")"; then
  exit 2
fi
mapfile -t RESOLVED_CASE_NAMES <<< "$RESOLVED_OUTPUT"
if [[ "${#RESOLVED_CASE_NAMES[@]}" -eq 0 || -z "${RESOLVED_CASE_NAMES[0]}" ]]; then
  echo "ERROR: no cases resolved from $CASES_JSON" >&2
  exit 2
fi

for c in "${RESOLVED_CASE_NAMES[@]}"; do
  CMD=(
    "$SCRIPT_DIR/05_run_single_comparison_case.sh"
    --case "$c"
    --cases-json "$CASES_JSON"
    --bootstrap-reps "$BOOTSTRAP_REPS"
    --high-value-frame-stride "$HV_FRAME_STRIDE"
    --workspace-root "$WORKSPACE_ROOT"
    --results-root "$RESULTS_ROOT"
  )
  for root in "${SEARCH_ROOTS[@]}"; do
    CMD+=(--search-root "$root")
  done
  "${CMD[@]}"
done

echo "All cases completed."
echo "See: $RESULTS_ROOT"
