# Figures 1–8

Figure 1 is drawn in TikZ. Figures 2–8 use R and `tikzDevice`. All figures use
Source Sans 3 for text and Computer Modern for mathematics. Font files and
their license are included in `fonts/`; `fontspec` uses `no-math` to preserve
the mathematical fonts.

From the repository directory, run:

```sh
./run_all.sh --figures
```

The PDFs are written to `figures/fig1_geometry.pdf` through
`figures/fig8_ocrs.pdf`. Requirements and optional PDF checks are described in
the [running guide](../docs/running-the-analyses.md#figures-and-output-comparisons).

## Inputs

The figure scripts use saved simulation and aggregate empirical results,
including the files in `analysis/inputs/figure-reference/`. They do not refit
models or read participant records. Running the analyses writes new results to
`analysis/outputs/` and does not replace these plotting inputs.

The generator checks input files against `manifest/figure-source-lock.json`
and compares the plotted values with their expected summaries. It stops if
an input has changed. If replacing an input deliberately, first check it
against the paper's results and then update its checksum. The comparison
files are written to `analysis/outputs/figure-audit/`.

## Editing a figure

Edit `analysis/scripts/07_figures.R` for Figures 2–8 and
`fig_src/fig1_geometry_cm.tex` for Figure 1. Files in `fig_src/generated-cm/`
are generated and will be overwritten. To build selected figures:

```sh
Rscript analysis/scripts/07_figures.R --figures2-8
python3 fig_src/build_fig1.py
```

The R generator can also be called by its absolute path from another working
directory. It uses `python3` unless `PYTHON` specifies another executable.
TeX and Poppler commands are found through `PATH`.

`./run_all.sh --figure-qa` checks embedded fonts, page dimensions, text size,
and text outside the page, and renders previews in `figures/png/`. All figures
are 6.5 inches wide, with heights set for the manuscript. Inspect the previews
after editing; the numerical and PDF checks cannot detect every layout problem.
