#!/usr/bin/env bash
# The G3 production run, in the preregistered order: the decisive ions first
# (each: the anchor at theta0, then the four axes with the interior of each
# axis last), the control ion, then parts (b) and (c). One log per unit.
# STRICTLY SEQUENTIAL: two 1e6-packet jobs side by side on this 24 GB WSL box
# corrupted index arrays (bit 60 flipped on otherwise valid indices,
# 2026-09-24) and killed a third job silently. Resumable: a state with a
# completion marker is skipped, so the chain restarts where it died.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
while [ -f paperB/gate3/.wait_for ] && [ ! -f "$(cat paperB/gate3/.wait_for)" ]; do sleep 30; done
for ion in 58CeII 60NdII 57LaII; do
  [ -f paperB/gate3/gate3_ref_${ion}.done ] || $PY paperB/gate3/run_gate3.py --ion $ion --axis ref > paperB/gate3/run_${ion}_ref.log 2>&1
  for axis in T D J P; do
    $PY paperB/gate3/run_gate3.py --ion $ion --axis $axis --skip-done >> paperB/gate3/run_${ion}_${axis}.log 2>&1
  done
done
[ -f paperB/gate3/gate3_partb_blend3.done ] || $PY paperB/gate3/run_gate3.py --part b > paperB/gate3/run_partb.log 2>&1
[ -f paperB/gate3/gate3_partc_p1blend.done ] || $PY paperB/gate3/run_gate3.py --part c > paperB/gate3/run_partc.log 2>&1
echo "G3 run finished $(date -u +%FT%TZ)" > paperB/gate3/run_all.finished
