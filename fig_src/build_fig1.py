"""Build Figure 1 with Source Sans 3 text and Computer Modern mathematics."""

from pathlib import Path
import os
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "fig_src/fig1_geometry_cm.tex"
BUILD = ROOT / "fig_src/build-fig1-cm"
BUILD.mkdir(exist_ok=True)
for executable in ("latexmk", "xelatex", "pdffonts", "pdftoppm"):
    if shutil.which(executable) is None:
        raise SystemExit(
            f"{executable} not found; install TeX Live/MacTeX or Poppler and add its binaries to PATH"
        )
(ROOT / "figures/png").mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
command = [
    "latexmk",
    "-norc",
    "-xelatex",
    "-interaction=nonstopmode",
    "-halt-on-error",
    "-outdir=" + str(BUILD),
    str(SOURCE),
]
result = subprocess.run(
    command,
    cwd=ROOT,
    env=env,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
(BUILD / "compile-console.txt").write_text(result.stdout)
if result.returncode:
    raise SystemExit(result.stdout[-6000:])
pdf = BUILD / "fig1_geometry_cm.pdf"
fonts = subprocess.check_output(["pdffonts", str(pdf)], text=True)
assert "CMMI" in fonts and "CMR" in fonts and "CMSY" in fonts, fonts
assert "SourceSans3-Regular" in fonts and "SourceSans3-Bold" in fonts, fonts
assert not any(name in fonts for name in ["Arial", "NimbusSans", "+Symbol "]), fonts
(BUILD / "fonts.txt").write_text(fonts)
shutil.copy2(pdf, ROOT / "figures/fig1_geometry.pdf")
subprocess.run(
    [
        "pdftoppm",
        "-singlefile",
        "-r",
        "144",
        "-png",
        str(pdf),
        str(ROOT / "figures/png/fig1_geometry"),
    ],
    check=True,
)
print(
    "Built Figure 1: Source Sans 3 text and Computer Modern mathematics."
)
