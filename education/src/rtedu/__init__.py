"""rtedu: the toy radiative-transfer code of the education book.

Level-2 code: small, reusable, transparent, and deliberately independent of
the production package `sobolev/`. Nothing here is publication evidence; the
demos teach and validate concepts that the production code implements with
real atomic data. Only the bridge chapter (E4) imports production code.

Every chapter's notebook seeds its random generator from SEEDS so that the
figures, the GIFs and data/results.json are byte-stable across builds.
"""
from pathlib import Path

EDU_DIR = Path(__file__).resolve().parents[2]      # education/
DATA_DIR = EDU_DIR / "data"
FIG_DIR = EDU_DIR / "figures"
VIDEO_DIR = EDU_DIR / "videos"

SEEDS = {"ch00": 100, "ch01": 101, "ch02": 102, "ch03": 103, "ch04": 104}

# toy constants (CGS), written out so the toy formulas read without imports
C = 2.99792458e10          # cm s^-1
E_ESU = 4.80320425e-10     # esu
M_E = 9.1093837e-28        # g
SIGMA_CLASSICAL = 3.141592653589793 * E_ESU ** 2 / (M_E * C)   # pi e^2 / (m_e c), cm^2 s^-1
DAY = 86400.0
