"""Golden packet histories of `run_mc`, every mode, taken before Paper IV
touched the transport (WP0 of the Paper IV plan).

`paper3/freeze.py --check` regenerates derived tables from *committed*
transport outputs and never re-runs `run_mc`, so until now nothing guarded
the transport numerically except the physics tests, which are statistical.
Paper IV adds energy packets, a downward macroatom, a re-emitting core and a
CSR branching sampler to `forest_mc.py`; each must leave every existing mode
bit-identical when its keyword sits at the default. This file pins that:
SHA-256 of the packet histories (escape frequencies, fates, weights, event
and first-line counters) for every mode on two toy atoms, classical and
worldline, plus the line-memory variants of the grouped legs.

Regenerate ONLY when a change is meant to move the histories (state why in
the commit):
  2026-09-10  the bin legs (expansion_*, binned_*, dual_*) re-pinned after
              the off-by-one in `nu_of_G` was fixed (Paper IV Phase 8,
              results_report 4.55); the Sobolev legs did not move. The
              pre-fix digests are kept in golden_run_mc_prefix_2026-09-10.json
              for provenance and are tested by nothing.

    .venv/bin/python tests/test_golden_run_mc.py --regen
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "tests", ROOT / "paper2/phase1", ROOT / "paper2/phase0/three_level_atom",
          ROOT / "paper3", ROOT / "paper3/synthetic"):
    sys.path.insert(0, str(p))

from forest import synthetic_forest            # noqa: E402
from forest_mc import MODES, run_mc            # noqa: E402

# The Paper II/III modes whose histories are pinned. Paper IV modes (dmacro,
# binned/dual thermal) are tested in their own files; they never enter here.
PAPER3_MODES = ("sobolev_group", "sobolev_absorb", "expansion_absorb", "sobolev_thermal",
                "expansion_thermal", "sobolev_branch", "sobolev_tla", "expansion_tla",
                "expansion_branch", "binned_group", "binned_absorb", "expansion_group",
                "dual_branch", "dual_group", "dual_absorb")
assert set(PAPER3_MODES) <= set(MODES)
from redistribution import RedistributionKernel  # noqa: E402
import test_forest_mc as tfm                   # noqa: E402

GOLDEN = ROOT / "tests/data/golden_run_mc.json"
N_PACKETS = 10000
HASHED = ("nu_out_all", "fate", "w", "n_events", "first_line")


def _three_level_case():
    fa, _ = tfm.three_level(1.5)
    fa.temperature = 3000.0
    fa.emis_w = np.array([1.0, 1.0])       # the toy's upper level is unpopulated
    lo, hi = tfm.pump_band()
    return dict(atom=fa, r_core=tfm.R_CORE, r_out=tfm.R_OUT, t_exp=tfm.T_EXP,
                lo=lo, hi=hi, t_core=None, kernel=tfm._toy_kernel())


def _synthetic_case():
    r_core, r_out, t_exp, t_core = 8.64e13, 2.592e14, 86400.0, 6000.0
    atom, _ = synthetic_forest(n_lines=60, tau=6.0, span=0.25, n_exit=2, dlnlam=0.05, seed=3)
    lo, hi = atom.op_nu.min() * 0.99, atom.op_nu.max() * 1.01
    ref = run_mc(atom, r_core, r_out, t_exp, lo, hi, 20000, "sobolev_branch",
                 seed=1, t_core=t_core, collect_events=True)
    e = ref["events"]
    kern = RedistributionKernel.from_branching_mc(e[0], e[1], np.ones(e[0].size), 16,
                                                  nu_lo=lo, nu_hi=hi)
    return dict(atom=atom, r_core=r_core, r_out=r_out, t_exp=t_exp,
                lo=lo, hi=hi, t_core=t_core, kernel=kern)


CASES = {"three_level": _three_level_case, "synthetic": _synthetic_case}


def variants():
    """(case, mode, relativity, line_memory) for every pinned run."""
    out = []
    for case in CASES:
        for mode in PAPER3_MODES:
            for rel in (None, "worldline"):
                out.append((case, mode, rel, False))
                if mode.endswith("_group") and rel is None:
                    out.append((case, mode, rel, 1))
    return out


def _digest(res):
    h = hashlib.sha256()
    for key in HASHED:
        h.update(np.ascontiguousarray(res[key]).tobytes())
    return h.hexdigest()


def run_variant(setup, mode, rel, mem):
    kw = {"kernel": setup["kernel"]} if mode.endswith("_group") else {}
    res = run_mc(setup["atom"], setup["r_core"], setup["r_out"], setup["t_exp"],
                 setup["lo"], setup["hi"], N_PACKETS, mode, seed=1, eps=0.5,
                 t_core=setup["t_core"], relativity=rel, line_memory=mem, **kw)
    a = res["accounting"]
    return dict(sha256=_digest(res), steps=int(res["steps"]),
                n_escaped=res["n_escaped"], n_core=res["n_core"],
                n_absorbed=res["n_absorbed"], n_interactions=res["n_interactions"],
                E_esc_over_inj=a["E_esc"] / a["E_inj"],
                E_dep_lab_over_inj=a["E_dep_lab"] / a["E_inj"])


def key_of(case, mode, rel, mem):
    return f"{case}|{mode}|{rel or 'classical'}|mem={int(mem)}"


def compute():
    setups = {c: f() for c, f in CASES.items()}
    return {key_of(*v): run_variant(setups[v[0]], *v[1:]) for v in variants()}


@pytest.fixture(scope="module")
def golden():
    if not GOLDEN.exists():
        pytest.skip("golden_run_mc.json not present")
    return json.loads(GOLDEN.read_text())


@pytest.fixture(scope="module")
def setups():
    return {c: f() for c, f in CASES.items()}


@pytest.mark.parametrize("case,mode,rel,mem", variants())
def test_history_matches_golden(golden, setups, case, mode, rel, mem):
    got = run_variant(setups[case], mode, rel, mem)
    exp = golden["runs"][key_of(case, mode, rel, mem)]
    assert got["sha256"] == exp["sha256"], (case, mode, rel, mem, got, exp)
    assert got["steps"] == exp["steps"]


def test_golden_covers_every_mode(golden):
    assert set(golden["runs"]) == {key_of(*v) for v in variants()}
    assert golden["n_packets"] == N_PACKETS and tuple(golden["hashed"]) == HASHED


if __name__ == "__main__":
    if "--regen" not in sys.argv:
        sys.exit("pass --regen to rewrite tests/data/golden_run_mc.json")
    import subprocess
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                          text=True).stdout.strip()
    runs = compute()
    GOLDEN.write_text(json.dumps(dict(git_head=head, n_packets=N_PACKETS, hashed=HASHED,
                                      seed=1, eps=0.5, runs=runs), indent=1) + "\n")
    print(f"wrote {GOLDEN} ({len(runs)} runs at {head[:10]})")
