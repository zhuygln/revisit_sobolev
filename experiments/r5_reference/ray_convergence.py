"""Why Delta_Sobolev needs both legs on the same rays -- and the cure.

Row 1 (the pathology): the analytic leg on its 200 midpoint rays against the
resolved leg on n_core in {6, 12, 24, 48, 96} -- the mismatch that made
Delta_Sob change sign in the high-beta pilot.
Row 2 (the cure): both legs on identical midpoint rays, n in {50..1600}, and
Gauss-Legendre for comparison. Delta_Sob(matched) must be flat.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference import *  # noqa

lines, n_ion = forest_lines(stim=True)
nu = nu_grid(1600)

print("La II forest, v_D = 100 km/s, band 3800-3955 A, erf resolved leg\n")
print("=== pathology: Sobolev on 200 midpoint rays vs resolved on n_core rays ===")
print(" n_core   F_sob(200)   F_res(n)    Delta_Sob")
sob200 = legs_erf(lines, n_ion, RaySet.midpoint(R_CORE, R_OUT, 200), nu=nu)["sob_first"]
for n in (6, 12, 24, 48, 96, 200):
    res = legs_erf(lines, n_ion, RaySet.midpoint(R_CORE, R_OUT, n), nu=nu)["res_first"]
    print(f"  {n:4d}     {sob200:.5f}     {res:.5f}    {100*delta(sob200, res):+7.3f}%")

print("\n=== cure: both legs on the same rays ===")
print(" rule        n     F_sob      F_res      Delta_Sob(first)   Delta_Sob(classical)   Delta_exp")
rows = []
for rule, ns in (("midpoint", (50, 100, 200, 400, 800, 1600)), ("gauss", (64, 128))):
    for n in ns:
        rs = (RaySet.midpoint if rule == "midpoint" else RaySet.gauss_legendre)(R_CORE, R_OUT, n)
        L = legs_erf(lines, n_ion, rs, nu=nu)
        d1, d0, de = delta(L["sob_first"], L["res_first"]), delta(L["sob"], L["res"]), delta(L["exp"], L["res_first"])
        rows.append((rule, n, d1))
        print(f"  {rule:9s} {n:5d}   {L['sob_first']:.5f}   {L['res_first']:.5f}      {100*d1:+8.4f}%           {100*d0:+8.4f}%         {100*de:+7.3f}%")
mid = [d for r, n, d in rows if r == "midpoint" and n >= 100]
print(f"\nmatched midpoint n>=100: Delta_Sob spread = {100*(max(mid)-min(mid)):.4f} points")
