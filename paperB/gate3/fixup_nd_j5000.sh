#!/usr/bin/env bash
# Nd II's J axis gave up after 4 attempts on 2026-09-24 (the host's repeatable
# bit-60 index fault) with the interior state J5000 never completing -- the
# chain moved on to P and La II with that state missing. This fills the gap
# after the main chain finishes (never two large jobs at once), with its own
# retry loop, then leaves a marker for the analysis to notice.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
while [ ! -f paperB/gate3/run_all.finished ]; do sleep 60; done
for attempt in 1 2 3 4 5; do
  echo "=== attempt $attempt $(date -u +%FT%TZ)" >> paperB/gate3/fixup_nd_j5000.log
  $PY paperB/gate3/run_gate3.py --ion 60NdII --axis J --state 5000 >> paperB/gate3/fixup_nd_j5000.log 2>&1 \
    && [ -f paperB/gate3/gate3_J_J5000_60NdII.done ] && break
  echo "=== attempt $attempt FAILED $(date -u +%FT%TZ)" >> paperB/gate3/fixup_nd_j5000.log; sleep 10
done
echo "fixup finished $(date -u +%FT%TZ)" > paperB/gate3/fixup_nd_j5000.finished
