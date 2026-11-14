"""Verify that the published and relabeled figures have identical vector plots.

Run with a Python interpreter that provides PyMuPDF.
"""

from pathlib import Path
import hashlib
import json
import sys
import fitz

qa = Path(sys.argv[1])
records = {}
for revised in sorted((qa / "figures").glob("*.pdf")):
    name = revised.name
    with fitz.open(qa / "before" / name) as old, fitz.open(revised) as new:
        assert old[0].rect == new[0].rect

        # The four legend sample lines occupy the final ten percent of the page.
        def paths(page):
            return [
                str(p)
                for p in page.get_drawings()
                if p["rect"].y0 < 0.90 * page.rect.height
            ]

        before = paths(old[0])
        after = paths(new[0])
        assert before == after, name
        label_found = False
        for block in new[0].get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    if "Linear ANCOVA (RCT)" in span["text"]:
                        assert new[0].rect.contains(fitz.Rect(span["bbox"]))
                        label_found = True
        assert label_found, name
        records[name] = {
            "nonlegend_vector_paths_identical": True,
            "path_count": len(before),
            "paths_sha256": hashlib.sha256("\n".join(before).encode()).hexdigest(),
            "new_label_fully_inside_page": True,
        }
(qa / "vector-validation.json").write_text(json.dumps(records, indent=2) + "\n")
print(
    "Verified unchanged plot vectors and complete legend labels for",
    len(records),
    "figures.",
)
