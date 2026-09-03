#!/bin/bash
# launched once the redo's run_grid.py (PID $1) exits: the remaining three T-scale runs in parallel
cd "$(dirname "$0")/../.."
while kill -0 "$1" 2>/dev/null; do sleep 60; done
P=../../.venv/bin/python
OMP_NUM_THREADS=1 nohup $P -u grid.py --mass 0.01 --v 0.1 --xlan 0.01 --t-scale 1.25 --out grid/tscale/model_M0.01_v0.1_X0.01_T1.25.json > grid/tscale/T1.25.log 2>&1 &
OMP_NUM_THREADS=1 nohup $P -u grid.py --mass 0.01 --v 0.1 --xlan 0.01 --t-scale 0.8 --t-scale-gas --out grid/tscale/model_M0.01_v0.1_X0.01_T0.8_gas.json > grid/tscale/T0.8_gas.log 2>&1 &
OMP_NUM_THREADS=1 nohup $P -u grid.py --mass 0.01 --v 0.1 --xlan 0.01 --t-scale 1.25 --t-scale-gas --out grid/tscale/model_M0.01_v0.1_X0.01_T1.25_gas.json > grid/tscale/T1.25_gas.log 2>&1 &
wait
