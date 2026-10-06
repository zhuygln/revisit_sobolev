#!/usr/bin/env bash
# The G2R production run: Ce II then Nd II, N_g = 2, 4, 8, 16, 32, one
# process per (ion, N_g), strictly sequential and retried (the bit-60 host
# fault of the G3 run is still open; a crash is caught and the unit retried
# up to 4 times; a unit with a completion marker is skipped).
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
echo $$ > paperB/gate2r/run_all.pid
try() {  # try <ion> <k>
  local ion=$1 k=$2 marker=paperB/gate2r/g2r_${1}_k${2}.done log=paperB/gate2r/run_${1}_k${2}.log
  for attempt in 1 2 3 4; do
    [ -f "$marker" ] && return 0
    echo "=== attempt $attempt $(date -u +%FT%TZ): $ion k=$k" >> "$log"
    $PY paperB/gate2r/run_g2r.py --ion $ion --k $k --skip-done >> "$log" 2>&1 && [ -f "$marker" ] && return 0
    echo "=== attempt $attempt FAILED $(date -u +%FT%TZ)" >> "$log"; sleep 10
  done
  echo "GAVE UP: $ion k=$k" >> paperB/gate2r/run_all.log; return 1
}
for ion in 58CeII 60NdII; do
  for k in 2 4 8 16 32; do try $ion $k; done
done
echo "G2R chain finished $(date -u +%FT%TZ)" > paperB/gate2r/run_all.finished
