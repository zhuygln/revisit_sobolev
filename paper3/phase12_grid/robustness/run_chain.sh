#!/bin/bash
# usage: run_chain.sh <log> <model:t> [<model:t> ...]   -- serial chain runs, one cell after another
cd "$(dirname "$0")/.."
log=$1; shift
for cell in "$@"; do
  m=${cell%%:*}; t=${cell##*:}
  echo "=== $(date) start $m t=$t" >> "$log"
  ../../.venv/bin/python -u robustness.py chain --model "$m" --t "$t" --chain-max 2000 4000 8000 --wall-budget 7200 >> "$log" 2>&1
  echo "=== $(date) done $m t=$t" >> "$log"
done
