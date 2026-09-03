"""Run the 27-point (M_ej, v_ej, X_lan) grid as a pool of `grid.py` subprocesses.

Heaviest models first (sorted by M X / v^3, the shell's n_ion at fixed epoch),
one process per model with OMP_NUM_THREADS=1, restartable: a model whose JSON
carries `complete: true` is skipped. `--dry-run` builds every atom at every
epoch without transport and prints n_ion, n_opacity, tau_max, S and the RSS --
run it first; the Nd II list alone is 3.3 M lines per process.

Redo (§4.41): `--redo over_budget --budget 5400` re-runs every cell whose row
is `over_budget`, one process per (model, epoch) into `grid/redo/<model>_t<t>.json`
(logs appended there); the model JSONs are untouched until `--merge-redo`
replaces the rows and tags them `redo = {budget_s, git, file}`.

Usage: python run_grid.py [--workers 8] [--n 300000] [--budget 1500] [--dry-run]
                          [--timeout 28800] [--only M,v,X]
       python run_grid.py --redo over_budget --budget 5400 --workers 8
       python run_grid.py --merge-redo
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
GRID_DIR = HERE / "grid"
M_GRID = (0.003, 0.01, 0.03)
V_GRID = (0.05, 0.1, 0.2)
X_GRID = (1e-3, 1e-2, 1e-1)
PY = sys.executable


def points():
    pts = [(m, v, x) for m in M_GRID for v in V_GRID for x in X_GRID]
    return sorted(pts, key=lambda p: -p[0] * p[2] / p[1] ** 3)


def is_complete(path):
    try:
        return json.loads(path.read_text()).get("complete", False)
    except Exception:
        return False


def run_one(m, v, x, n, budget, dry, timeout, out_dir):
    from grid import model_name
    name = model_name(m, v, x)
    out = out_dir / f"{name}.json"
    log = out_dir / f"{name}.log"
    if not dry and is_complete(out):
        return name, "skipped (complete)", 0.0
    cmd = [PY, str(HERE / "grid.py"), "--mass", str(m), "--v", str(v), "--xlan", str(x),
           "--n", str(n), "--budget", str(budget), "--out", str(out)] + (["--dry-run"] if dry else [])
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    t0 = time.time()
    with open(log, "w") as f:
        try:
            rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env,
                                timeout=timeout).returncode
            status = "ok" if rc == 0 else f"rc={rc}"
        except subprocess.TimeoutExpired:
            status = "timeout"
    return name, status, time.time() - t0


def redo_cells(out_dir, status="over_budget"):
    """(m, v, x, t_d, projected_s) for every row of the given status, heaviest first."""
    cells = []
    for p in sorted(out_dir.glob("model_*.json")):
        d = json.loads(p.read_text())
        if d.get("t_scale", 1.0) != 1.0:
            continue
        for r in d["rows"]:
            if r.get("status") == status:
                cells.append((d["m_ej_msun"], d["v_ej_c"], d["x_lan"], r["t_d"],
                              r.get("projected_s", 0.0)))
    return sorted(cells, key=lambda c: -c[4])


def run_cell(m, v, x, t, n, budget, timeout, out_dir):
    """One (model, epoch) into out_dir/redo, log appended."""
    from grid import model_name
    name = f"{model_name(m, v, x)}_t{t:g}"
    redo_dir = out_dir / "redo"; redo_dir.mkdir(exist_ok=True)
    out, log = redo_dir / f"{name}.json", redo_dir / f"{name}.log"
    cmd = [PY, str(HERE / "grid.py"), "--mass", str(m), "--v", str(v), "--xlan", str(x),
           "--n", str(n), "--budget", str(budget), "--epochs", str(t), "--out", str(out)]
    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    t0 = time.time()
    with open(log, "a") as f:
        f.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(cmd)}\n"); f.flush()
        try:
            rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env,
                                timeout=timeout).returncode
            status = "ok" if rc == 0 else f"rc={rc}"
        except subprocess.TimeoutExpired:
            status = "timeout"
    return name, status, time.time() - t0


def merge_redo(out_dir, verbose=True):
    """Replace each model's row by the transported row of its redo file."""
    merged = []
    for p in sorted((out_dir / "redo").glob("model_*_t*.json")):
        d = json.loads(p.read_text())
        rows = [r for r in d["rows"] if r.get("status") in ("ok", "reduced_n")]
        if not rows:
            continue
        model = out_dir / f"{p.stem.rsplit('_t', 1)[0]}.json"
        md = json.loads(model.read_text())
        for r in rows:
            r = dict(r, redo={"budget_s": d["budget_s"], "git": d.get("git"), "file": p.name,
                              "previous_status": next((q.get("status") for q in md["rows"]
                                                       if q["t_d"] == r["t_d"]), None)})
            md["rows"] = [r if q["t_d"] == r["t_d"] else q for q in md["rows"]]
            merged.append((model.name, r["t_d"], r["status"], r["n_used"]))
        model.write_text(json.dumps(md, indent=1))
    if verbose:
        for m in merged:
            print("  merged %-32s t=%-4g %-10s n_used=%d" % m)
    return merged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--n", type=float, default=3e5)
    ap.add_argument("--budget", type=float, default=1500.0)
    ap.add_argument("--timeout", type=float, default=8 * 3600)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", default=None, help="comma list M,v,X to run a single point")
    ap.add_argument("--out-dir", default=str(GRID_DIR))
    ap.add_argument("--redo", default=None, help="re-run every cell with this row status")
    ap.add_argument("--merge-redo", action="store_true")
    a = ap.parse_args()
    out_dir = Path(a.out_dir); out_dir.mkdir(exist_ok=True)
    if a.merge_redo:
        merge_redo(out_dir)
        return
    if a.redo:
        cells = redo_cells(out_dir, a.redo)
        print(f"{len(cells)} {a.redo} cells, {a.workers} workers, n={int(a.n)}, "
              f"budget={a.budget}s", flush=True)
        t0 = time.time()
        with ThreadPoolExecutor(a.workers) as ex:
            futs = [ex.submit(run_cell, m, v, x, t, int(a.n), a.budget, a.timeout, out_dir)
                    for m, v, x, t, _ in cells]
            for f in as_completed(futs):
                name, status, dt = f.result()
                print(f"  {name:40s} {status:12s} {dt:7.0f}s   [{time.time()-t0:.0f}s elapsed]",
                      flush=True)
        print(f"finished {len(cells)} cells in {time.time()-t0:.0f}s", flush=True)
        return
    pts = points()
    if a.only:
        m, v, x = (float(s) for s in a.only.split(","))
        pts = [(m, v, x)]
    print(f"{len(pts)} models, {a.workers} workers, n={int(a.n)}, budget={a.budget}s"
          f"{', DRY RUN' if a.dry_run else ''}", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(run_one, m, v, x, int(a.n), a.budget, a.dry_run, a.timeout, out_dir)
                for m, v, x in pts]
        for f in as_completed(futs):
            name, status, dt = f.result()
            print(f"  {name:32s} {status:20s} {dt:7.0f}s   [{time.time()-t0:.0f}s elapsed]", flush=True)
    print(f"finished {len(pts)} models in {time.time()-t0:.0f}s", flush=True)
    if a.dry_run:
        summarize_dry(out_dir)


def summarize_dry(out_dir):
    print(f"\n{'model':30s}{'t':>5}{'T_eff':>7}{'v_ph':>6}{'n_ion(Ce)':>11}{'n_op':>9}"
          f"{'tau_max':>9}{'S':>10}{'RSS MB':>8}")
    for p in sorted(out_dir.glob("model_*.json")):
        d = json.loads(p.read_text())
        for r in d["rows"]:
            print(f"{p.stem:30s}{r['t_d']:5.1f}{r['T_gas']:7.0f}{r['v_core']/d['v_ej_c']:6.3f}"
                  f"{r['n_ion'].get('58CeII', 0):11.2e}{r['n_opacity']:9d}{r['tau_max']:9.1e}"
                  f"{r.get('band_S_band', float('nan')):10.1f}{r['rss_mb_atom']:8.0f}")


if __name__ == "__main__":
    main()
