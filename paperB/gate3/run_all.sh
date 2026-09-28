#!/usr/bin/env bash
# The G3 production run, preregistered order: the decisive ions first (each:
# the anchor at theta0, then the four axes with each interior last), the
# control ion, the gray-remedy reruns at 1e6, then parts (b) and (c).
#
# STRICTLY SEQUENTIAL AND RETRIED. This 24 GB WSL box shows a repeatable
# memory fault: on 2026-09-24 four runs died with an IndexError whose index
# was a valid one with BIT 60 SET (2^60 + 1761, +596894, +472, +1103), in
# four code paths, twice with nothing else running. A crash is caught here
# and the unit retried (up to 4 attempts); a completed record is protected
# by the energy identity and kernel validation, which a flipped bit in a
# weight would break. One job at a time; a state with a completion marker
# is skipped (--skip-done), so the chain resumes where it died.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
echo $$ > paperB/gate3/run_all.pid
try() {  # try <marker> <log> <args...>: run until the marker exists, at most 4 attempts
  local marker=$1 log=$2; shift 2
  for attempt in 1 2 3 4; do
    [ -f "$marker" ] && return 0
    echo "=== attempt $attempt $(date -u +%FT%TZ): $*" >> "$log"
    $PY paperB/gate3/run_gate3.py "$@" >> "$log" 2>&1 && [ -f "$marker" ] && return 0
    echo "=== attempt $attempt FAILED $(date -u +%FT%TZ)" >> "$log"; sleep 10
  done
  echo "GAVE UP: $*" >> paperB/gate3/run_all.log; return 1
}
axis() {  # axis <ion> <axis>: resumable, at most 4 attempts
  local ion=$1 ax=$2 log=paperB/gate3/run_${1}_${2}.log
  for attempt in 1 2 3 4; do
    echo "=== attempt $attempt $(date -u +%FT%TZ)" >> "$log"
    $PY paperB/gate3/run_gate3.py --ion $ion --axis $ax --skip-done >> "$log" 2>&1 && return 0
    echo "=== attempt $attempt FAILED $(date -u +%FT%TZ); resuming" >> "$log"; sleep 10
  done
  echo "GAVE UP: $ion $ax" >> paperB/gate3/run_all.log; return 1
}
for ion in 58CeII 60NdII 57LaII; do
  try paperB/gate3/gate3_ref_${ion}.done paperB/gate3/run_${ion}_ref.log --ion $ion --axis ref
  for ax in T D J P; do axis $ion $ax; done
done
# Ce II's gray-remedy reruns (precision rule at 3e5): the whole state at 1e6
try paperB/gate3/gate3_D_D10_58CeII_n1e6.done   paperB/gate3/rerun_58CeII_D10_n1e6.log   --ion 58CeII --axis D --state 10   --n 1000000
try paperB/gate3/gate3_T_T5000_58CeII_n1e6.done paperB/gate3/rerun_58CeII_T5000_n1e6.log --ion 58CeII --axis T --state 5000 --n 1000000
try paperB/gate3/gate3_P_P1d_58CeII_n1e6.done   paperB/gate3/rerun_58CeII_P1_n1e6.log    --ion 58CeII --axis P --state 1    --n 1000000
try paperB/gate3/gate3_P_P3d_58CeII_n1e6.done   paperB/gate3/rerun_58CeII_P3_n1e6.log    --ion 58CeII --axis P --state 3    --n 1000000
try paperB/gate3/gate3_partb_blend3.done  paperB/gate3/run_partb.log --part b
try paperB/gate3/gate3_partc_p1blend.done paperB/gate3/run_partc.log --part c
echo "G3 chain finished $(date -u +%FT%TZ)" > paperB/gate3/run_all.finished
