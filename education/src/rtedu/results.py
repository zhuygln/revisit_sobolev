"""One source of numbers for the book (PI amendment 4: transactional).

Each notebook ends with `record("chNN", {...})`, which writes
data/generated/chNN.json; `merge()` builds data/results.json from those
files deterministically; the chapters read only results.json through
`load("chNN")`. `make education-notebooks` deletes both before executing
chapter 1 onward, so results.json is always the current build.
"""
import json
from pathlib import Path

from . import DATA_DIR

GENERATED = DATA_DIR / "generated"
RESULTS = DATA_DIR / "results.json"
SIG = 6


def _round(obj):
    """Round every float to SIG significant digits so the JSON is stable
    against last-bit noise in BLAS; lists and dicts are walked."""
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, (int, str)) or obj is None:
        return obj
    if isinstance(obj, float):
        return float(f"{obj:.{SIG}g}") if obj == obj else None
    if isinstance(obj, dict):
        return {str(k): _round(v) for k, v in obj.items()}
    if hasattr(obj, "tolist"):
        return _round(obj.tolist())
    if isinstance(obj, (list, tuple)):
        return [_round(v) for v in obj]
    return obj


def record(chapter, values, generated=GENERATED):
    """Write the chapter's generated values (rounded, sorted keys)."""
    generated = Path(generated); generated.mkdir(parents=True, exist_ok=True)
    p = generated / f"{chapter}.json"
    p.write_text(json.dumps(_round(values), indent=1, sort_keys=True) + "\n")
    return p


def merge(generated=GENERATED, results=RESULTS):
    """Build results.json from every data/generated/*.json, in name order."""
    generated, results = Path(generated), Path(results)
    out = {}
    for p in sorted(generated.glob("*.json")):
        out[p.stem] = json.loads(p.read_text())
    results.parent.mkdir(parents=True, exist_ok=True)
    results.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    return out


def load(chapter, results=RESULTS, generated=GENERATED):
    """The chapter's values from results.json (the merged, current build).
    During a build the merge has not happened yet, so a later notebook that
    needs an earlier one's values reads that chapter's generated file."""
    results = Path(results)
    if results.exists():
        d = json.loads(results.read_text())
        if chapter in d:
            return d[chapter]
    return json.loads((Path(generated) / f"{chapter}.json").read_text())
