#!/usr/bin/env python3
"""Two builds must agree: results.json, figures and GIFs byte for byte; the
HTML animations by their embedded frames (PI amendment 6). Run after a
second `make education` with the first build's hashes saved by
`--save`: `check_determinism.py --save` then rebuild then `check_determinism.py`."""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "src"))
from rtedu.visualization import animation_frame_hashes   # noqa: E402

STATE = HERE / "data" / ".determinism.json"


def snapshot():
    out = {}
    for p in sorted(list((HERE / "figures").glob("*")) + list((HERE / "videos").glob("*.gif")) + [HERE / "data" / "results.json"]):
        out[str(p.relative_to(HERE))] = hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted((HERE / "videos").glob("*.html")):
        out[str(p.relative_to(HERE))] = "frames:" + hashlib.sha256("".join(animation_frame_hashes(p)).encode()).hexdigest()
    return out


def main():
    now = snapshot()
    if "--save" in sys.argv:
        STATE.write_text(json.dumps(now, indent=1, sort_keys=True)); print(f"saved {len(now)} hashes"); return 0
    if not STATE.exists():
        print("no saved hashes: run with --save first, rebuild, then run again"); return 0
    old = json.loads(STATE.read_text())
    diff = [k for k in sorted(set(old) | set(now)) if old.get(k) != now.get(k)]
    if diff:
        print("DETERMINISM CHECK FAILED:", file=sys.stderr)
        for k in diff:
            print("  - " + k, file=sys.stderr)
        return 1
    print(f"determinism check OK: {len(now)} artefacts identical (HTML animations by frames)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
