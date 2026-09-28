# thermal_sampler unclamped searchsorted -- bug snapshot (preserved before the fix)

Pre-fix commit: see pre_fix_commit.txt (8c21e07081b702257ddd8cb77ec6ecbb5101f2f5)

## Confirmed defect

`paper2/phase1/forest_mc.py`, `ForestAtom.thermal_sampler` (line ~404-429):

    cum = np.cumsum(w / w.sum())
    def sample(u):
        return np.searchsorted(cum, u)

`cum[-1]` is not exactly 1.0 (floating-point drift over the sum of `w`,
sized to every emitting line -- 20,752,336 entries for Paper B's 13-ion
P1 blend). When a drawn `u` lands at or above `cum[-1]`, `searchsorted`
returns `len(cum)`, one past the end. Two sibling samplers in the SAME
FILE already guard exactly this case (`min(..., stop - 1)` at line 486;
`np.clip(..., 0, nb - 1)` at line 908); `thermal_sampler` alone omits it.

## Failing unit (100% reproducible, 4/4 identical attempts)

`paperB/gate3/run_gate3.py --part c` (the 13-ion P1 realistic-mixture
blend, `paperB/gate3/run_partc.log`), seeds (1, 2, 3), n = 1,000,000,
mode `sobolev_tla` (any eps leg with dead-end/thermal draws on this atom
reaches it; the crash surfaced during one of the `E{eps}` legs after
`E0.90` completed). Identical traceback and identical index on every
attempt (2026-09-24 22:28, 2026-09-25 05:57, 13:19, 20:38):

    Traceback (most recent call last):
      File ".../paperB/gate3/run_gate3.py", line 269, in <module>
        main()
      File ".../paperB/gate3/run_gate3.py", line 255, in main
        run_part_c(a.n or N_BLEND, seeds, bseeds, ng, eps_grid=eps, out=a.out); return
      File ".../paperB/gate3/run_gate3.py", line 221, in run_part_c
        row = L.run_legs(zone, atom, n, legs=specs, seeds=tuple(seeds), ng=int(ng_fine), verbose=verbose)
      File ".../paper4/phase2_energy/legs.py", line 259, in run_legs
        res = [mc(spec, s, collect_events=collect, **kw) for s in seeds_leg]
      File ".../paper4/phase2_energy/legs.py", line 168, in mc
        return run_mc(atom, zone["r_core"], zone["r_out"], zone["t_exp"], lo, hi, n, spec["mode"],
      File ".../paper2/phase1/forest_mc.py", line 1586, in run_mc
        else atom.tau_all[new_line[todo]])
    IndexError: index 20752336 is out of bounds for axis 0 with size 20752336

## Reachability audit (the mandatory check before touching the file)

`thermal_sampler` is called from four sites: the up-front `thermal =`/
`thermal_s =` construction (gated `needs_thermal = outcome in ("thermal",
"tla")` -- G3's per-state legs use `sobolev_dmacro`, so this site is NOT
reached there), and the k-packet/dead-end reprocessing inside the
`dmacro`/`macro` branch (`kpack = atom.thermal_sampler(..., weight=
"energy_beta")`, gated only on `thermal_k != "deposit"`, the project's
default is "reemit"). **This second site IS reached by every dmacro leg
with any dead-end macroatom walk** -- G1, G2, R2M and every G3 per-state
leg, all of which print nonzero `dead=` counts routinely in the hundreds
of thousands to millions.

Whether an overflow there could go UNDETECTED (a completed run silently
using a wrong line) was checked by tracing every use of `new_line`, not
assumed:
- for `outcome in ("dmacro", "macro")`, the escape-probability line that
  caught the part-c crash (1586, `atom.tau_all[new_line[todo]]`) is
  explicitly SKIPPED (`outcome not in ("dmacro", "macro")` guards it);
- but line 1624, `nu_rest[~is_bin & ~coherent] = atom.nu0_all[new_line[...]]`,
  runs UNCONDITIONALLY for every outcome and indexes `atom.nu0_all`,
  which is sized over the identical "every line" space `thermal_sampler`
  sums over.

An overflowed `new_line` value (== `len(cum)`, i.e. `n_lines_total`) is
therefore always out of range for `atom.nu0_all` too, and crashes at
1624 if not earlier -- for every outcome, not only "tla". **A run that
completed without an IndexError could not have hit this condition; no
completed unit needs to be revisited for this specific defect.**
