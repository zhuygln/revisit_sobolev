"""The SEDONA launch helpers: defaults must reproduce the historical behaviour.

These pin the contract the experiment drivers rely on. The defaults matter as
much as the overrides -- every committed SEDONA result was produced with one
rank and unscaled timeouts, so a change in either default would silently
invalidate comparisons against them.
"""
import os

import pytest

from sobolev.sedona import (n_ranks, sedona_cmd, sedona_exe, sedona_home,
                            sedona_timeout)

ENV = ("SEDONA_HOME", "SEDONA_EXE", "SEDONA_NP", "SEDONA_MPIRUN",
       "SEDONA_TIMEOUT_SCALE")


@pytest.fixture
def clean_env(monkeypatch):
    for k in ENV:
        monkeypatch.delenv(k, raising=False)


def test_defaults_are_the_historical_behaviour(clean_env):
    assert sedona_home() == os.path.expanduser("~/personal/pubsed")
    assert sedona_exe().endswith("/personal/pubsed/src/sedona6.ex")
    assert n_ranks() == 1
    assert sedona_cmd() == [sedona_exe(), "param.lua"]
    assert sedona_timeout(5000) == 5000.0


def test_home_and_exe_overrides(clean_env, monkeypatch):
    monkeypatch.setenv("SEDONA_HOME", "/opt/pubsed")
    assert sedona_exe() == "/opt/pubsed/src/sedona6.ex"
    # SEDONA_EXE wins over the home-derived path
    monkeypatch.setenv("SEDONA_EXE", "/elsewhere/sedona6.ex")
    assert sedona_exe() == "/elsewhere/sedona6.ex"
    assert sedona_cmd() == ["/elsewhere/sedona6.ex", "param.lua"]


def test_timeout_scales(clean_env, monkeypatch):
    monkeypatch.setenv("SEDONA_TIMEOUT_SCALE", "2.5")
    assert sedona_timeout(2000) == 5000.0
    assert sedona_timeout(560) == 1400.0


def test_multi_rank_is_opt_in(clean_env, monkeypatch):
    monkeypatch.setenv("SEDONA_NP", "8")
    assert n_ranks() == 8
    assert sedona_cmd() == ["mpirun", "-n", "8", sedona_exe(), "param.lua"]
    monkeypatch.setenv("SEDONA_MPIRUN", "/opt/ompi/bin/mpirun")
    assert sedona_cmd()[0] == "/opt/ompi/bin/mpirun"


def test_rank_count_is_floored_at_one(clean_env, monkeypatch):
    """A stray SEDONA_NP=0 must not produce `mpirun -n 0`."""
    monkeypatch.setenv("SEDONA_NP", "0")
    assert n_ranks() == 1
    assert sedona_cmd() == [sedona_exe(), "param.lua"]


def test_param_file_is_passed_through(clean_env):
    assert sedona_cmd("other.lua")[-1] == "other.lua"


def test_every_driver_uses_the_helper():
    """No driver may reconstruct the path or a raw timeout on its own."""
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    offenders = []
    for p in (root / "experiments").rglob("*.py"):
        s = p.read_text()
        if "subprocess.run" not in s or "sedona6.ex" not in s and "sedona_cmd" not in s:
            continue
        if "sedona_cmd" in s:
            if "yozhuz" in s or 'expanduser("~/personal/pubsed")' in s:
                offenders.append(f"{p.name}: hardcoded path")
            for line in s.splitlines():
                if "timeout=" in line and "sedona_timeout" not in line:
                    offenders.append(f"{p.name}: raw timeout -- {line.strip()[:60]}")
    assert not offenders, offenders
