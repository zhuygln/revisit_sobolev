"""Locating and launching SEDONA from the experiment drivers.

SEDONA lives outside this repo, and both *where* it is and *how long* it
takes are machine facts, not properties of an experiment. Before this module
the drivers hardcoded one machine's home directory (fixed 2026-08-29) and one
machine's wall-clock budgets, and the second of those is the subtler failure:
a timeout tuned on a fast box does not announce itself on a slow one, it just
kills the run partway with a TimeoutExpired and no result.

Environment:

    SEDONA_HOME           default ~/personal/pubsed
    SEDONA_EXE            default $SEDONA_HOME/src/sedona6.ex
    SEDONA_TIMEOUT_SCALE  multiplies every driver's timeout; default 1.0
    SEDONA_NP             MPI ranks; default 1
    SEDONA_MPIRUN         launcher used when SEDONA_NP > 1; default "mpirun"

The defaults reproduce the historical behaviour exactly: one rank, unscaled
timeouts, the binary under ~/personal/pubsed.

On ranks: the committed results are single-rank, and SEDONA seeds its RNG per
rank, so raising SEDONA_NP changes the random stream. Multi-rank runs are
comparable to the committed numbers only within Monte Carlo noise, never
bit-for-bit. That is why the default is 1 and why this is opt-in.
"""
import os

DEFAULT_HOME = "~/personal/pubsed"


def sedona_home():
    return os.environ.get("SEDONA_HOME", os.path.expanduser(DEFAULT_HOME))


def sedona_exe():
    return os.environ.get("SEDONA_EXE", f"{sedona_home()}/src/sedona6.ex")


def n_ranks():
    return max(1, int(os.environ.get("SEDONA_NP", "1")))


def sedona_cmd(param="param.lua"):
    """Argv for one SEDONA run, single-rank unless SEDONA_NP says otherwise."""
    n = n_ranks()
    if n == 1:
        return [sedona_exe(), param]
    return [os.environ.get("SEDONA_MPIRUN", "mpirun"), "-n", str(n),
            sedona_exe(), param]


def sedona_timeout(default):
    """A driver's own timeout, scaled by SEDONA_TIMEOUT_SCALE.

    The in-script values stay as written -- they carry each experiment's
    relative cost, which is real information -- and the scale factor adapts
    the whole set to a machine. Measured throughput on the 2026-08-29 box was
    ~695 s per 1e6 packets single-rank, which put vc_control's 1e7-packet run
    at ~6950 s against its written 5000 s budget: the first case that needed
    this.
    """
    return float(default) * float(os.environ.get("SEDONA_TIMEOUT_SCALE", "1"))
