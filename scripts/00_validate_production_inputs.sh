#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/production_comparison_cases.json"

if [[ ! -f "$CASES_JSON" ]]; then
  echo "ERROR: cases file not found: $CASES_JSON"
  exit 2
fi

summarize_dump() {
  local dump_file="$1"
  awk '
    BEGIN{n=0; prev=""; uniq=0}
    /^ITEM: TIMESTEP$/ {
      getline s
      n++
      if(n==1) first=s
      last=s
      if(prev!=""){
        d=s-prev
        if(!(d in cnt)){ord[++uniq]=d}
        cnt[d]++
      }
      prev=s
    }
    END {
      printf("frames=%d first_step=%s last_step=%s", n, first, last)
      if (uniq>0) {
        printf(" strides=")
        for(i=1;i<=uniq;i++){
          d=ord[i]
          printf("%s%s:%d", (i==1?"":";"), d, cnt[d])
        }
      }
      printf("\n")
    }
  ' "$dump_file"
}

check_log_completion() {
  local log_file="$1"
  local loop_line wall_line
  loop_line="$(grep -m1 'Loop time of' "$log_file" || true)"
  wall_line="$(grep -m1 'Total wall time:' "$log_file" || true)"
  if [[ -n "$loop_line" ]]; then
    echo "    $loop_line"
  else
    echo "    WARNING: no completion line in log"
  fi
  if [[ -n "$wall_line" ]]; then
    echo "    $wall_line"
  fi
}

python - "$CASES_JSON" "$REPO_DIR" <<'PY' | while IFS=$'\t' read -r name desc npbc_dump pbc_dump npbc_log pbc_log; do
import json
import sys
from pathlib import Path

cases_file = Path(sys.argv[1])
repo = Path(sys.argv[2])
obj = json.loads(cases_file.read_text())
for c in obj.get("cases", []):
    fields = [
        c["name"],
        c.get("description", ""),
        str((repo / c["npbc_prod_dump"]).resolve()),
        str((repo / c["pbc_prod_dump"]).resolve()),
        str((repo / c.get("npbc_prod_log", "")).resolve()) if c.get("npbc_prod_log") else "",
        str((repo / c.get("pbc_prod_log", "")).resolve()) if c.get("pbc_prod_log") else "",
    ]
    print("\t".join(fields))
PY

  echo "============================================================"
  echo "CASE: $name"
  echo "Description: $desc"

  echo "  NPBC dump: $npbc_dump"
  if [[ ! -f "$npbc_dump" ]]; then
    echo "  ERROR: missing NPBC dump"
    exit 2
  fi
  echo "    $(summarize_dump "$npbc_dump")"

  echo "  PBC dump:  $pbc_dump"
  if [[ ! -f "$pbc_dump" ]]; then
    echo "  ERROR: missing PBC dump"
    exit 2
  fi
  echo "    $(summarize_dump "$pbc_dump")"

  if [[ -n "$npbc_log" ]]; then
    echo "  NPBC log:  $npbc_log"
    if [[ -f "$npbc_log" ]]; then
      check_log_completion "$npbc_log"
    else
      echo "    WARNING: NPBC log not found"
    fi
  fi
  if [[ -n "$pbc_log" ]]; then
    echo "  PBC log:   $pbc_log"
    if [[ -f "$pbc_log" ]]; then
      check_log_completion "$pbc_log"
    else
      echo "    WARNING: PBC log not found"
    fi
  fi

done

echo "============================================================"
echo "Input validation complete."
