"""The education tests run under `make education-test` (pytest education/tests)
and stay out of the research suite (root testpaths = tests/). rtedu is not
installed: it runs from education/src (PI amendment 5)."""
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
