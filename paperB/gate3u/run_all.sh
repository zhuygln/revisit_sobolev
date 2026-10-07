#!/usr/bin/env bash
# The G3U production run: one process per affected frozen G3 record (the
# set is derived by analyse.py --affected, never typed), strictly
# sequential and retried (up to 4 attempts per record; a record with a
# completion marker is skipped). Starts only after the G2R chain has
# finished (one 1e6-packet job at a time on this box).
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
echo $$ > paperB/gate3u/run_all.pid
RECORDS=$($PY paperB/gate3u/analyse.py --affected | awk '/^gate3_/ {print $1}' | sort -u)
echo "G3U records: $(echo "$RECORDS" | wc -l)" > paperB/gate3u/run_all.log
try() {  # try <record>
  local rec=$1 stem=${1%.json} marker log
  marker=paperB/gate3u/${stem/gate3_/g3u_}.done; log=paperB/gate3u/run_${stem/gate3_/}.log
  for attempt in 1 2 3 4; do
    [ -f "$marker" ] && return 0
    echo "=== attempt $attempt $(date -u +%FT%TZ): $rec" >> "$log"
    $PY paperB/gate3u/run_g3u.py --record "$rec" --skip-done >> "$log" 2>&1 && [ -f "$marker" ] && return 0
    echo "=== attempt $attempt FAILED $(date -u +%FT%TZ)" >> "$log"; sleep 10
  done
  echo "GAVE UP: $rec" >> paperB/gate3u/run_all.log; return 1
}
for rec in $RECORDS; do try "$rec"; done
echo "G3U chain finished $(date -u +%FT%TZ)" > paperB/gate3u/run_all.finished
