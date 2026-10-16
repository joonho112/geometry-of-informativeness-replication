"""Compile one R-generated TikZ figure with Source Sans 3 text and CM mathematics."""

from pathlib import Path
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
source = Path(sys.argv[1]).resolve()
build = source.parent / "build"
build.mkdir(exist_ok=True)
for executable in ("latexmk", "xelatex", "pdffonts", "pdftoppm"):
    if shutil.which(executable) is None:
        raise SystemExit(
            f"{executable} not found; install TeX Live/MacTeX or Poppler and add its binaries to PATH"
        )
(ROOT / "figures/png").mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
result = subprocess.run(
    [
        "latexmk",
        "-norc",
        "-xelatex",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-outdir=" + str(build),
        source.name,
    ],
    cwd=source.parent,
    env=env,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
(build / (source.stem + "-compile.txt")).write_text(result.stdout)
if result.returncode:
    raise SystemExit(result.stdout[-6000:])
pdf = build / (source.stem + ".pdf")
fonts = subprocess.check_output(["pdffonts", str(pdf)], text=True)
assert "CMMI" in fonts, fonts
assert (
    "SourceSans3" in fonts and "Arial" not in fonts and "NimbusSan" not in fonts
), fonts
assert "+Symbol " not in fonts, fonts
(build / (source.stem + "-fonts.txt")).write_text(fonts)
shutil.copy2(pdf, ROOT / "figures" / pdf.name)
print("Built " + pdf.name + ": Source Sans 3 text and Computer Modern mathematics.")
