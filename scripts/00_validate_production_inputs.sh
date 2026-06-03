#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_JSON="$REPO_DIR/configs/production_comparison_cases.json"
CASES=()
SEARCH_ROOTS=()
STRICT_LOGS=0

usage() {
  cat <<EOF
Usage: $0 --search-root DIR [options]

Required:
  --search-root DIR       Directory to search for trajectories/logs. Repeatable.

Options:
  --case NAME             Validate only one configured case. Repeatable.
  --cases-json FILE       Case config file. Default: $CASES_JSON
  --strict-logs           Treat missing optional logs as errors.
  -h, --help              Show this help.
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
    --strict-logs)
      STRICT_LOGS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ ! -f "$CASES_JSON" ]]; then
  echo "ERROR: cases file not found: $CASES_JSON"
  exit 2
fi

if [[ ${#SEARCH_ROOTS[@]} -eq 0 ]]; then
  echo "ERROR: choose at least one directory with --search-root" >&2
  usage >&2
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

RESOLVE_CMD=(python3 "$SCRIPT_DIR/resolve_production_inputs.py" --cases-json "$CASES_JSON" --format tsv)
for root in "${SEARCH_ROOTS[@]}"; do
  RESOLVE_CMD+=(--search-root "$root")
done
for case_name in "${CASES[@]}"; do
  RESOLVE_CMD+=(--case "$case_name")
done
if [[ "$STRICT_LOGS" -eq 1 ]]; then
  RESOLVE_CMD+=(--strict-logs)
fi

"${RESOLVE_CMD[@]}" | while IFS=$'\t' read -r name desc npbc_dump pbc_dump npbc_log pbc_log; do

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
