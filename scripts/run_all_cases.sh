#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/cases.json"
BOOTSTRAP_REPS="${1:-1000}"
HV_FRAME_STRIDE="${2:-1}"

mapfile -t CASES < <(python3 - "$CASES_JSON" <<'PY'
import json
import sys
from pathlib import Path
cases = json.loads(Path(sys.argv[1]).read_text()).get("cases", [])
for c in cases:
    print(c["name"])
PY
)

for c in "${CASES[@]}"; do
  "$SCRIPT_DIR/05_run_case_pipeline.sh" "$c" "$BOOTSTRAP_REPS" "$HV_FRAME_STRIDE"
done

echo "All cases completed."
echo "See: $REPO_DIR/results"
