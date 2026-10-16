#!/usr/bin/env python3
"""Check figure dimensions, fonts, and text placement, and save PNG previews.

Requires pypdf, pdfplumber, and Poppler."""
from pathlib import Path
import hashlib
import json
import subprocess
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import pdfplumber
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis/outputs/figure-audit"
receipt = json.loads((OUT / "generation-receipt.json").read_text())
names = list(receipt["figures"])


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check(name):
    spec = receipt["figures"][name]
    f = ROOT / spec["pdf"]
    proof = ROOT / "figures/png" / name
    subprocess.run(
        ["pdftoppm", "-singlefile", "-r", "144", "-png", str(f), str(proof)], check=True
    )
    font_output = subprocess.check_output(["pdffonts", str(f)], text=True)
    rows = font_output.splitlines()[2:]
    assert rows and all(
        row.split()[-5] == "yes" for row in rows
    ), f"Unembedded font: {name}"
    font_names = [row.split()[0].split("+")[-1] for row in rows]
    assert any(n.startswith("CMMI") for n in font_names), f"Missing CM math: {name}"
    assert "Symbol" not in font_names, f"Legacy Symbol math fallback: {name}"
    assert any(
        n.startswith("SourceSans3") for n in font_names
    ), f"Missing Source Sans 3: {name}"
    allowed = ("CMMI", "CMR", "CMSY", "CMEX", "MSBM", "SourceSans3")
    assert all(n.startswith(allowed) for n in font_names), (name, font_names)
    pdf = PdfReader(f)
    assert len(pdf.pages) == 1
    p = pdf.pages[0]
    width, height = float(p.mediabox.width), float(p.mediabox.height)
    assert abs(width - 468) < 1e-6 and abs(height / 72 - spec["height_inches"]) < 1e-6
    assert height <= 6.3 * 72 + 1e-6
    sizes = Counter()

    def visitor(text, cm, tm, font, size):
        if text.strip():
            # Account for the text and graphics matrices, including rotated labels.
            scale = math.hypot(
                tm[0] * cm[0] + tm[1] * cm[2], tm[0] * cm[1] + tm[1] * cm[3]
            )
            sizes[round(size * scale, 3)] += len(text.strip())

    p.extract_text(visitor_text=visitor)
    # Native TeX uses 72.27 pt/in; PDF uses 72 bp/in. Standard 12 pt
    # text/math and TeX's native 8 pt scripts are 11.955 and 7.970 bp.
    allowed_sizes = {round(pt * 72 / 72.27, 3) for pt in (8, 12)}
    assert set(sizes).issubset(allowed_sizes), (name, sizes)
    assert round(12 * 72 / 72.27, 3) in sizes, (name, sizes)
    with pdfplumber.open(f) as d:
        page = d.pages[0]
        outside = [
            c
            for c in page.chars
            if c["x0"] < -0.05
            or c["x1"] > width + 0.05
            or c["top"] < -0.05
            or c["bottom"] > height + 0.05
        ]
        assert not outside, f"Text outside media box: {name}"
    (OUT / f"{name}-fonts.txt").write_text(font_output)
    return dict(
        name=name,
        width_inches=width / 72,
        height_inches=height / 72,
        pdf_sha256=sha(f),
        proof_png=str(proof.relative_to(ROOT)) + ".png",
        proof_sha256=sha(proof.with_suffix(".png")),
        embedded_fonts=True,
        text_outside_media_box=0,
        font_names=font_names,
        math_family="Computer Modern",
        font_sizes_pt=dict(sorted(sizes.items())),
        page_count=1,
    )


with ThreadPoolExecutor(max_workers=4) as ex:
    results = list(ex.map(check, names))
assert all(
    sha(ROOT / p) == v for p, v in receipt["sources"].items()
), "Source input changed"
result = dict(
    status="PASS",
    figures=results,
    proof_dpi=144,
    base_text_tex_pt=12,
    mathematical_scripts_tex_pt=8,
    note="Font-size checks use PDF text operators and convert TeX points to PDF points; rotated-label bounding-box height is not font size. Source Sans 3 ordinary text and native CM mathematical fonts are required.",
)
(OUT / "pdf-qa.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        dict(
            status="PASS",
            figure_count=len(results),
            width_inches=6.5,
            all_fonts_embedded=True,
            text_outside_media_box=0,
        )
    )
)
