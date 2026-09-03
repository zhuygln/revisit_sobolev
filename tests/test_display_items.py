"""docs/paper3/display_items.py: every display item renders from the frozen
JSONs into a temporary directory at a journal column width, without a
timestamp in the PDF."""
import importlib.util, os, sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS3 = ROOT / "docs" / "paper3"
FROZEN = ROOT / "paper3" / "FROZEN.json"
needs_frozen = pytest.mark.skipif(not FROZEN.exists(), reason="FROZEN.json not present")


def _load():
    spec = importlib.util.spec_from_file_location("paper3_display_items", DOCS3 / "display_items.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pdf_size_mm(path):
    import re
    m = re.search(rb"/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]", path.read_bytes())
    assert m, "no MediaBox"
    w = (float(m.group(3)) - float(m.group(1))) * 25.4 / 72
    return w


@needs_frozen
@pytest.mark.parametrize("which", ["fig1", "fig2", "fig3", "fig4", "edfig1", "edfig2", "edfig3", "edfig4"])
def test_item_renders_at_column_width(tmp_path, which):
    os.environ.setdefault("SOURCE_DATE_EPOCH", "1700000000")
    sys.path.insert(0, str(ROOT / "paper3"))
    import freeze
    di = _load()
    written = di.main(freeze.canonical(), tmp_path, which=[which])
    pdfs = [Path(p) for p in written if str(p).endswith(".pdf")]
    pngs = [Path(p) for p in written if str(p).endswith(".png")]
    assert len(pdfs) == 1 and len(pngs) == 1
    data = pdfs[0].read_bytes()
    assert b"/CreationDate" not in data and b"/ModDate" not in data
    # saved with a tight bounding box, so the page is the content width: never
    # wider than the 180 mm double column, and drawn at that scale
    w = _pdf_size_mm(pdfs[0])
    assert 130 < w <= 180.5, w


def test_every_manuscript_figure_has_a_generator():
    import re
    di = _load()
    tex = (DOCS3 / "manuscript.tex").read_text()
    used = set(re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", tex))
    generated = {di.NAMES[k] for k in di.ITEMS} if hasattr(di, "NAMES") else None
    if generated is None:
        pytest.skip("display_items has no NAMES map")
    assert used == generated
