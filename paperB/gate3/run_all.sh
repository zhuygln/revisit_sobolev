#!/usr/bin/env bash
# The G3 production run, in the preregistered order: the decisive ions first
# (each: the anchor at theta0, then the four axes with the interior of each
# axis last), the control ion, then parts (b) and (c). One log per unit.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
for ion in 58CeII 60NdII 57LaII; do
  $PY paperB/gate3/run_gate3.py --ion $ion --axis ref > paperB/gate3/run_${ion}_ref.log 2>&1
  for axis in T D J P; do
    $PY paperB/gate3/run_gate3.py --ion $ion --axis $axis > paperB/gate3/run_${ion}_${axis}.log 2>&1
  done
done
$PY paperB/gate3/run_gate3.py --part b > paperB/gate3/run_partb.log 2>&1
$PY paperB/gate3/run_gate3.py --part c > paperB/gate3/run_partc.log 2>&1
echo "G3 run finished $(date -u +%FT%TZ)" > paperB/gate3/run_all.finished
