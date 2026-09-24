#!/usr/bin/env python3
"""G3: freeze one frequency support per ion over the whole state domain --
the union of every state's opacity range with G1's 0.5 % margin -- and the
three-ion union for part (b). Written to gate3_support.json BEFORE the run;
every kernel of G3 (anchor, recomputed, interpolated) is laid on it.

    .venv/bin/python paperB/gate3/support.py
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import states as S                                          # noqa: E402

OUT = HERE / "gate3_support.json"


def main():
    out = dict(rule="union over the G3 domain of [op_nu.min x 0.995, op_nu.max x 1.005]; per ion; part (b) uses the three-ion union",
               per_ion={}, per_state={})
    for ion in S.IONS:
        lo, hi = float("inf"), 0.0
        for axis, value in S.all_states():
            st = S.build_state(axis, value, ion)
            a = st["atom"]
            r = dict(nu_lo=float(a.op_nu.min()), nu_hi=float(a.op_nu.max()), n_opacity=int(a.n_opacity),
                     tau_max=float(a.op_tau.max()), n_lines=int(a.n_lines_total), T_gas=st["T_gas"], t_core=st["t_core"],
                     n_ion=st["n_ion"], coord=st["coord"], coord_value=st["coord_value"])
            out["per_state"][f"{st['label']}_{ion}"] = r
            lo = min(lo, r["nu_lo"] * 0.995); hi = max(hi, r["nu_hi"] * 1.005)
            print(f"{ion} {st['label']:7s} T {st['T_gas']:.0f} t_core {st['t_core']:.0f} n {st['n_ion']:.3e}: "
                  f"{a.n_opacity:7d} opacity lines, nu {r['nu_lo']:.3e}..{r['nu_hi']:.3e}, tau_max {r['tau_max']:.0f}", flush=True)
        out["per_ion"][ion] = dict(nu_lo=lo, nu_hi=hi)
    out["blend3"] = dict(nu_lo=min(v["nu_lo"] for v in out["per_ion"].values()), nu_hi=max(v["nu_hi"] for v in out["per_ion"].values()))
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", OUT)
    for ion, v in out["per_ion"].items():
        print(f"  {ion}: {v['nu_lo']:.4e} .. {v['nu_hi']:.4e}")
    print(f"  blend3: {out['blend3']['nu_lo']:.4e} .. {out['blend3']['nu_hi']:.4e}")


if __name__ == "__main__":
    main()
