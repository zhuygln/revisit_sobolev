"""Freeze the Paper III analysis: regenerate every derived table, JSON and
figure from the committed transport outputs, then record what was used and
what came out in `paper3/FROZEN.json`.

    python paper3/freeze.py            regenerate everything, write FROZEN.json
    python paper3/freeze.py --check    verify the committed state against FROZEN.json
    python paper3/freeze.py --check --strict    figure-hash mismatches are errors too

The freeze imports the drivers (no subprocess) and runs with the working
directory set to paper3/phase12_grid, so the override paths recorded in
`sensitivity_chain8000.json` stay the relative strings
`robustness/chain_model_*_t[0-9].json`. Driver stdout goes to
`<scratch>/freeze_*.log`.

FROZEN.json holds: the git HEAD, whether any *input* is modified in the
working tree, the `git rev-parse HEAD:<dir>` tree hashes of the transport
output directories (identical before and after the freeze commit -- HEAD is
the parent), the UTC time, sha256 of every input (transport JSONs, passbands,
driver sources) and of every derived output, and the flat `headline` dict the
manuscript quotes through `docs/paper3/numbers.tex`.

`--check` runs four tiers: (a) input hashes; (b) hashes of the committed
outputs (hand-edited?); (c) regenerate everything into a scratch directory and
compare every JSON numerically (rel 1e-9, abs 1e-12, NaN == NaN, `cond` at
rel 1e-6); (d) recompute `headline()` from the regenerated set and compare to
`FROZEN.headline`. Figures are compared by hash and a mismatch is a warning
unless `--strict` (matplotlib output is byte-stable within one environment,
not across versions).
"""
import argparse, ast, contextlib, hashlib, importlib.util, json, math, os, subprocess, sys, tempfile
from datetime import datetime, timezone
from pathlib import Path

# byte-stable PDF/PNG output: matplotlib honours SOURCE_DATE_EPOCH for the PDF CreationDate
os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")

HERE = Path(__file__).resolve().parent          # paper3/
ROOT = HERE.parent
P12 = HERE / "phase12_grid"
P13 = HERE / "phase13_observability"
FIGDIR = HERE / "figures"
DOCS3 = ROOT / "docs" / "paper3"
FROZEN = HERE / "FROZEN.json"
TAG = "paper3-freeze"
CENTRAL = (0.01, 0.1, 0.01)
CENTRAL_KEY = str(CENTRAL)
SCRATCH = Path(os.environ.get("FREEZE_SCRATCH", tempfile.gettempdir())) / "paper3_freeze"

# derived JSONs, relative to a destination directory that mirrors phase12_grid
DERIVED = {
    "sensitivity": "sensitivity.json", "sensitivity_floored_incl": "sensitivity_floored_incl.json",
    "sensitivity_absorbing": "sensitivity_absorbing.json", "sensitivity_T1": "sensitivity_T1.json",
    "sensitivity_T2": "sensitivity_T2.json", "sensitivity_T3": "sensitivity_T3.json",
    "sensitivity_chain8000": "sensitivity_chain8000.json", "chain_table": "robustness/chain_table.json",
    "grid_table": "grid_table.json", "tscale": "tscale.json", "syserr": "syserr.json",
    "observability": "observability.json",
}
FIGS = ("fig2_bol_vs_colour", "fig4_vectors", "fig5_observability", "fig6_tangent", "fig7_tscale")
OVERRIDE_CHAIN = 8000
SOURCES = ("phase12_grid/sensitivity.py", "phase12_grid/robustness.py", "phase12_grid/grid_table.py",
           "phase12_grid/syserr.py", "phase12_grid/tscale.py", "phase12_grid/figures.py",
           "phase12_grid/grid.py", "phase12_grid/run_grid.py", "phase13_observability/observe.py")
LOOSE_KEYS = {"cond": 1e-6}          # near-singular normal matrices: BLAS-dependent


def canonical():
    """Destination paths of the committed derived files."""
    d = {k: P12 / v for k, v in DERIVED.items()}
    d["observability"] = P13 / "observability.json"
    d["figdir"] = FIGDIR
    return d


def scratch_dest(root):
    root = Path(root)
    d = {k: root / v for k, v in DERIVED.items()}
    d["figdir"] = root / "figures"
    for p in d.values():
        p.parent.mkdir(parents=True, exist_ok=True)
    d["figdir"].mkdir(exist_ok=True)
    return d


# ---------------------------------------------------------------- generation
def _drivers():
    sys.path.insert(0, str(P12)); sys.path.insert(0, str(P13)); sys.path.insert(0, str(ROOT))
    import sensitivity, robustness, grid_table, syserr, tscale, figures, observe   # noqa: E402
    return sensitivity, robustness, grid_table, syserr, tscale, figures, observe


def generate(dest, log_dir=None):
    """Run every driver into `dest` (a dict from `canonical()` or `scratch_dest()`)."""
    log_dir = Path(log_dir or SCRATCH); log_dir.mkdir(parents=True, exist_ok=True)
    S, R, G, Y, T, F, O = _drivers()
    grid_dir = P12 / "grid"
    chain_files = sorted(Path("robustness").glob("chain_model_*_t[0-9].json"))   # relative: cwd is P12
    assert len(chain_files) == 4, chain_files
    with open(log_dir / "freeze_drivers.log", "w") as log, contextlib.redirect_stdout(log):
        S.main(grid_dir, out=dest["sensitivity"])
        S.main(grid_dir, out=dest["sensitivity_floored_incl"], floored="include")
        S.main(grid_dir, out=dest["sensitivity_absorbing"], core="absorbing")
        for tg in ("T1", "T2", "T3"):
            S.main(grid_dir, out=dest[f"sensitivity_{tg}"], tangent=tg)
        S.main(grid_dir, out=dest["sensitivity_chain8000"], override=[str(f) for f in chain_files],
               override_chain=OVERRIDE_CHAIN)
        Path(dest["chain_table"]).write_text(json.dumps(R.chain_summary(), indent=1))
        Path(dest["grid_table"]).write_text(json.dumps(G.summary(G.load(grid_dir)), indent=1))
        obs = O.main(grid_dir, dest["observability"])
        O.fig5(obs, dest["figdir"], leg="C_both", central=CENTRAL_KEY)
        ts = T.main(grid_dir, grid_dir / "tscale", dest["tscale"])
        T.fig7(ts, dest["figdir"], noise_floor=ts.get("noise_floor"))
        F.fig2(grid_dir / f"{F.model_name(*CENTRAL)}.json", dest["figdir"])
        F.fig4(Path(dest["sensitivity"]), CENTRAL_KEY, dest["figdir"])
        F.fig6(Path(dest["sensitivity"]).parent, CENTRAL_KEY, dest["figdir"])
        Y.main(dest["sensitivity"], dest["syserr"], residual=dest["sensitivity_T1"])
    h = headline(dest)
    paper = paper_items(dest, h, log_dir)
    return h, paper


def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def paper_items(dest, h, log_dir):
    """Commit 10+: LaTeX number macros / tables and display items, if docs/paper3 exists.
    Returns the list of files written (for the output manifest)."""
    written = []
    lt, di = DOCS3 / "latex_tables.py", DOCS3 / "display_items.py"
    out_dir = DOCS3 if dest is canonical_dest else Path(dest["figdir"]).parent / "paper3"
    if lt.exists():
        written += _load_module(lt, "paper3_latex_tables").main(h, dest, out_dir)
    if di.exists():
        with open(Path(log_dir) / "freeze_display_items.log", "w") as log, contextlib.redirect_stdout(log):
            written += _load_module(di, "paper3_display_items").main(dest, out_dir / "figures")
    return [Path(w) for w in written]


canonical_dest = None   # set in main() so paper_items can tell the real run from a check


# ------------------------------------------------------------------ headline
def _j(p):
    return json.loads(Path(p).read_text())


def _by_X(points, leg="C_both", cls="C-B"):
    """`cls` count per X, as {X: [k, n]} from the per-point classes (summary has no by-X)."""
    out = {}
    for key, rec in points.items():
        x = ast.literal_eval(key)[2]
        r = rec["legs"][leg]
        k, n = out.setdefault(f"{x:g}", [0, 0])
        out[f"{x:g}"] = [k + (r.get("cls") == cls), n + (r.get("status") == "ok" and not r.get("underdetermined"))]
    return out


def headline(dest):
    """The flat dict of every number the manuscript quotes (docs/paper3/numbers.tex)."""
    T0, T1, T2, T3 = (_j(dest[f"sensitivity{s}"]) for s in ("", "_T1", "_T2", "_T3"))
    g, ct, ob, ts, se = (_j(dest[k]) for k in ("grid_table", "chain_table", "observability", "tscale", "syserr"))
    ov = _j(dest["sensitivity_chain8000"])
    h = {}
    s0 = T0["summary"]["C_both"]
    h.update({"gate2.T0.cb": s0["C-B"], "gate2.T0.n": s0["n"], "gate2.T0.median_R": s0["median_R"],
              "gate2.T0.median_chi2_res_dof": s0["median_chi2_res_dof"], "gate2.T0.cb_by_X": _by_X(T0["points"]),
              "gate2.T0.C_binned.cb": T0["summary"]["C_binned"]["C-B"]})
    s1 = T1["summary"]["C_both"]
    h.update({"gate2.T1.cb": s1["C-B"], "gate2.T1.cb_by_X": _by_X(T1["points"]),
              "gate2.T1.median_chi2_res_dof": s1["median_chi2_res_dof"], "gate2.T1.leftover_gt4": s1["leftover_chi2_res_dof_gt_4"],
              "gate2.T1.n_determined": s1["n_determined"],
              "gate2.T1.abs_aL_median": s1["abs_a_L"]["median"], "gate2.T1.abs_aL_max": s1["abs_a_L"]["max"]})
    s2 = T2["summary"]["C_both"]
    h.update({"gate2.T2.cb": s2["C-B"], "gate2.T2.cb_by_X": _by_X(T2["points"]),
              "gate2.T2.median_chi2_res_dof": s2["median_chi2_res_dof"], "gate2.T2.leftover_gt4": s2["leftover_chi2_res_dof_gt_4"],
              "gate2.T2.n_determined": s2["n_determined"], "gate2.T2.aT_median": s2["a_T_bb"]["median"]})
    s3 = T3["summary"]["C_both"]
    lin = [r["legs"]["C_both"].get("lin_2c") for r in T3["points"].values()]
    lin = [v for v in lin if v is not None]
    h.update({"gate2.T3.cb": s3["C-B"], "gate2.T3.cb_by_X": _by_X(T3["points"]),
              "gate2.T3.median_chi2_res_dof": s3["median_chi2_res_dof"], "gate2.T3.leftover_gt4": s3["leftover_chi2_res_dof_gt_4"],
              "gate2.T3.n_determined": s3["n_determined"], "gate2.T3.underdetermined": s3["underdetermined"],
              "gate2.T3.lin2c_gt1": [sum(v > 1 for v in lin), len(lin)]})
    h.update({"control.A_redist.cc": T0["summary"]["A_redist"]["C-C"], "control.A_redist.n": T0["summary"]["A_redist"]["n"],
              "control.B_opacity.cb": T0["summary"]["B_opacity"]["C-B"],
              "thresholds.chi2_small": T0["thresholds"]["chi2_small"], "thresholds.R_max": T0["thresholds"]["R_max"],
              "thresholds.signif_min": T0["thresholds"]["signif_min"]})
    be = g["band_extremes"]["C_both"]
    h.update({"grid.n_models": g["n_models"], "grid.n_cells": len(g["cells"]),
              "grid.n_cells_ran": sum(c["ran"] for c in g["cells"]), "grid.n_redone_cells": g["n_redone_cells"],
              "colour.worst_per_model_range": g["worst_dcolor_range"], "colour.g_dm_range": [be["g"]["min"], be["g"]["max"]],
              "colour.g_n": be["g"]["n"], "colour.K_dm_max": be["K"]["max"], "colour.K_n": be["K"]["n"],
              "colour.nir_negative": g["nir_negative"],
              "floor.n_well": g["floor"]["n_well"], "floor.n_min": g["floor"]["n_min"], "floor.median": g["floor"]["median"],
              "floor.p90": g["floor"]["p90"], "floor.max": g["floor"]["max"], "floor.redone_range": g["floor"]["redone_range"]})
    sc = ob["summary"]["C_both"]
    for scn in ("dense", "sparse", "optical"):
        for tg in ("T0", "T1", "T2", "T3"):
            r = sc[scn][tg]
            h[f"gate3.{scn}.{tg}.survives"] = r["survives"]; h[f"gate3.{scn}.{tg}.eligible"] = r["eligible"]
            h[f"gate3.{scn}.{tg}.by_X"] = {x: [q["survives"], q["eligible"]] for x, q in r["by_X"].items()}
    h["gate3.dense.T0.median_nir_share"] = sc["dense"]["T0"]["median_nir_share"]
    from sobolev import photometry as _phot
    h["gate3.distance_Mpc"] = round(_phot.D_40MPC / 3.085677581e24)   # sobolev.photometry.D_40MPC: the AT2017gfo distance
    for tag, short in (("illumination_only", "illum"), ("with_T_gas", "gas")):
        c = ts["variants"][tag]["compare"]["C_both"]
        h[f"tscale.{short}.cos"] = c["cos"]; h[f"tscale.{short}.norm_ratio"] = c["norm_ratio"]
    cs = ct["summary"]
    h.update({"chain.n_cells": cs["n_cells"], "chain.base_cap": int(cs["base_cap"]), "chain.top_cap": int(cs["top_cap"]),
              "chain.trapped_range": cs["trapped_range_base"], "chain.trapped_range_top": cs["trapped_range_top"],
              "chain.dm_change_8000_range": cs["C_both"]["dm_change_range_top"],
              "chain.dm_ref_change_8000_range": cs["dm_ref_change_range_top"],
              "chain.signs_kept": cs["C_both"]["signs_kept_top"], "chain.criterion_met": cs["C_both"]["criterion_met_top"],
              "chain.rel_max": cs["rel_max"],
              "chain.C_binned.signs_kept": cs["C_binned"]["signs_kept_top"], "chain.C_binned.criterion_met": cs["C_binned"]["criterion_met_top"],
              "chain.cb_with_override": ov["summary"]["C_both"]["C-B"], "chain.n_with_override": ov["summary"]["C_both"]["n"],
              "chain.median_chi2_res_dof_override": ov["summary"]["C_both"]["median_chi2_res_dof"],
              "chain.median_R_override": ov["summary"]["C_both"]["median_R"]})
    for leg in ("C_both", "C_binned", "A_redist"):
        L = se["legs"][leg]
        h[f"syserr.one_mode.{leg}"] = L["one_mode"]["f1"]
        h[f"syserr.svd_filled.{leg}"] = L["one_mode"]["svd_filled"]
        h[f"syserr.null_median.{leg}"] = L["null_sign_scramble"]["median"]
        h[f"syserr.null_p95.{leg}"] = L["null_sign_scramble"]["p95"]
        h[f"syserr.frac_gt_0p5.{leg}"] = L["frac_gt"]["0.5"]
        h[f"syserr.frac_gt_1.{leg}"] = L["frac_gt"]["1"]
        h[f"syserr.sign_pattern.{leg}"] = L["sign_pattern_gK"]
        h[f"syserr.chi2_equiv_1mag_median.{leg}"] = L["chi2_equiv"]["1"]["median"]
        h[f"syserr.chi2_equiv_1mag_range.{leg}"] = [L["chi2_equiv"]["1"]["min"], L["chi2_equiv"]["1"]["max"]]
        h[f"syserr.chi2_equiv_0p5mag_median.{leg}"] = L["chi2_equiv"]["0.5"]["median"]
        h[f"syserr.n_live.{leg}"] = L["n_live"]
    h["syserr.n_points"] = se["legs"]["C_both"]["n_points"]; h["syserr.n_keys"] = se["legs"]["C_both"]["n_keys"]
    h["syserr.mp_scale"] = se["legs"]["C_both"]["one_mode"]["mp_scale"]
    h["syserr.svd_n_keys"] = se["legs"]["C_both"]["one_mode"]["svd_n_keys"]
    if "residual" in se:
        Rr = se["residual"]["legs"]["C_both"]
        h.update({"syserr.residual_T1.one_mode": Rr["one_mode"]["f1"], "syserr.residual_T1.frac_gt_0p5": Rr["frac_gt"]["0.5"],
                  "syserr.residual_T1.chi2_equiv_1mag_median": Rr["chi2_equiv"]["1"]["median"]})
    return h


# ------------------------------------------------------------------ manifests
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rel(p):
    return str(Path(p).resolve().relative_to(ROOT))


def input_files():
    fs = sorted((P12 / "grid").glob("model_*.json")) + sorted((P12 / "grid" / "tscale").glob("model_*.json"))
    fs += sorted((P12 / "robustness").glob("chain_model_*.json"))
    fs += sorted((ROOT / "data" / "filters").glob("*.dat"))
    fs += [HERE / s for s in SOURCES] + [ROOT / "sobolev" / "photometry.py"]
    return fs


def output_files(dest, extra=()):
    fs = [Path(dest[k]) for k in DERIVED]
    fs += [Path(dest["figdir"]) / f"{n}.{ext}" for n in FIGS for ext in ("png", "pdf")]
    return fs + list(extra)


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def git_state():
    inputs = [rel(f) for f in input_files()]
    dirty = git("status", "--porcelain", "--", *inputs)
    trees = {d: git("rev-parse", f"HEAD:paper3/phase12_grid/{d}") for d in ("grid", "grid/tscale", "robustness")}
    trees["data/filters"] = git("rev-parse", "HEAD:data/filters")
    return {"head": git("rev-parse", "HEAD"), "inputs_dirty": bool(dirty), "dirty_inputs": dirty.splitlines(),
            "tree_dirty": bool(git("status", "--porcelain")), "trees": trees, "tag": TAG}


# ------------------------------------------------------------------ comparison
def close(a, b, rel_tol=1e-9, abs_tol=1e-12):
    if isinstance(a, bool) != isinstance(b, bool):
        return False
    if isinstance(a, float) or isinstance(b, float):
        if not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
            return False
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        return math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol)
    return a == b


def compare(a, b, path="", loose=LOOSE_KEYS, diffs=None, limit=20):
    """Recursive numeric comparison of two JSON trees; returns the list of differences."""
    diffs = [] if diffs is None else diffs
    if len(diffs) >= limit:
        return diffs
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            diffs.append(f"{path}: keys differ {sorted(set(a) ^ set(b))[:6]}"); return diffs
        for k in a:
            compare(a[k], b[k], f"{path}/{k}", loose, diffs, limit)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            diffs.append(f"{path}: length {len(a)} vs {len(b)}"); return diffs
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, f"{path}[{i}]", loose, diffs, limit)
    else:
        key = path.rsplit("/", 1)[-1].split("[")[0]
        tol = loose.get(key, 1e-9)
        if not close(a, b, rel_tol=tol):
            diffs.append(f"{path}: {a!r} vs {b!r}")
    return diffs


def check_inputs(fz, root=ROOT):
    """Tier (a): every input in the manifest exists and hashes as frozen; no new input."""
    fails = []
    for p, h in fz["inputs"].items():
        f = root / p
        if not f.exists():
            fails.append(f"input missing: {p}")
        elif sha(f) != h:
            fails.append(f"input changed: {p}")
    for p in sorted({rel(f) for f in input_files()} - set(fz["inputs"])):
        fails.append(f"input not in manifest: {p}")
    return fails


def check_outputs(fz, root=ROOT, strict=False):
    """Tier (b): the committed outputs hash as frozen (a hand-edited table would not).
    Figures are warnings unless strict. Returns (fails, warns)."""
    fails, warns = [], []
    for p, h in fz["outputs"].items():
        f = root / p
        if not f.exists():
            fails.append(f"output missing: {p}")
        elif sha(f) != h:
            (warns if f.suffix in (".png", ".pdf") and not strict else fails).append(f"output differs from FROZEN: {p}")
    return fails, warns


def check_regenerated(fz, dest, sd, h_new, strict=False):
    """Tiers (c) and (d): the regenerated JSONs equal the committed ones numerically
    and the headline recomputed from them equals FROZEN.headline."""
    fails, warns = [], []
    for k in DERIVED:
        fails += [f"{DERIVED[k]}{x}" for x in compare(_j(dest[k]), _j(sd[k]))]
    for n in FIGS:
        for ext in ("png", "pdf"):
            a, b = Path(dest["figdir"]) / f"{n}.{ext}", Path(sd["figdir"]) / f"{n}.{ext}"
            if sha(a) != sha(b):
                (fails if strict else warns).append(f"figure regenerated differently: {n}.{ext}")
    fails += [f"headline{x}" for x in compare(fz["headline"], h_new)]
    return fails, warns


def check(strict=False, scratch=None):
    fz = _j(FROZEN)
    dest = canonical()
    fails = check_inputs(fz)
    f, warns = check_outputs(fz, strict=strict)
    fails += f
    scratch = Path(scratch or SCRATCH / "check")
    sd = scratch_dest(scratch)
    h_new, _ = generate(sd, log_dir=scratch)
    f, w = check_regenerated(fz, dest, sd, h_new, strict)
    fails += f; warns += w
    for w in warns:
        print("WARN", w)
    for f in fails:
        print("FAIL", f)
    print(f"freeze check: {len(fails)} failures, {len(warns)} warnings "
          f"(HEAD {fz['git']['head'][:10]}, frozen {fz['utc']})")
    return not fails


def main():
    global canonical_dest
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--strict", action="store_true", help="figure-hash mismatches fail the check")
    ap.add_argument("--scratch", default=None, help="where --check regenerates (default: FREEZE_SCRATCH or tmp)")
    a = ap.parse_args()
    os.chdir(P12)
    canonical_dest = canonical()
    if a.check:
        sys.exit(0 if check(a.strict, a.scratch) else 1)
    dest = canonical_dest
    h, paper = generate(dest, log_dir=a.scratch)
    import numpy, matplotlib
    fz = {"tag": TAG, "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "git": git_state(),
          "versions": {"python": sys.version.split()[0], "numpy": numpy.__version__, "matplotlib": matplotlib.__version__},
          "source_date_epoch": os.environ["SOURCE_DATE_EPOCH"],
          "inputs": {rel(f): sha(f) for f in input_files()},
          "outputs": {rel(f): sha(f) for f in output_files(dest, paper)}, "headline": h}
    FROZEN.write_text(json.dumps(fz, indent=1))
    print(f"wrote {FROZEN}: {len(fz['inputs'])} inputs, {len(fz['outputs'])} outputs, {len(h)} headline numbers; "
          f"HEAD {fz['git']['head'][:10]}, inputs dirty: {fz['git']['inputs_dirty']}")
    for k in ("gate2.T0.cb", "gate2.T1.cb_by_X", "colour.worst_per_model_range", "floor.median", "gate3.dense.T0.survives",
              "chain.dm_change_8000_range", "syserr.one_mode.C_both", "syserr.null_median.C_both"):
        print(f"  {k}: {h[k]}")


if __name__ == "__main__":
    main()
