# The Geometry of Informativeness: replication code

Replication code for *The Geometry of Informativeness: A Bregman-Centroid
Diagnostic for Partial Pooling in Multisite Trials* (Lee, 2026).

The scripts reproduce the simulations, Online Credit Recovery Study (OCRS)
analyses, and figures.
The companion R package [poolinfo](https://github.com/joonho112/poolinfo) provides
the general methods; it is not needed to run these scripts.

## Requirements

The analyses were tested with R 4.6.0 on macOS. You will need:

- R packages `jsonlite`, `digest`, `ggplot2`, `patchwork`, and `tikzDevice`;
- Python 3;
- TeX Live or MacTeX with `latexmk`, `xelatex`, TikZ, `standalone`, `preview`,
  `fontspec`, `xunicode`, `amsmath`, and `amssymb`;
- Poppler's `pdftoppm` and `pdffonts` for the figures.

The TeX and Poppler commands must be on `PATH`. Source Sans 3 fonts are included.
To check R dependencies, run `Rscript setup.R --check`. Use
`Rscript setup.R --install` to install missing R packages. The analysis scripts
do not install packages or download data. For the tested package versions, see
[`renv.lock`](renv.lock) and the [setup details](docs/running-the-analyses.md#setup).

## Run the analyses

Open a terminal in the repository directory. To check that the supplied inputs
match those used for the paper, run:

```sh
./run_all.sh --verify
```

This check needs only Python 3. To reproduce the numerical checks, simulation
summaries, empirical posterior and sensitivity summaries, and all eight figures,
run:

```sh
./run_all.sh --reproduce
```

The command uses saved simulation losses and aggregate OCRS inputs. It does
not need participant data or refit all 810 simulation cells. It compares the
recomputed results with the paper's saved results and stops if a check fails.
Figures use saved plotting inputs, so running the figure step alone does not
repeat the analyses.

Results are written to `analysis/outputs/`, figure PDFs to `figures/`, and
comparison reports to `verification/generated/`. The manuscript is not included
or compiled. The [running guide](docs/running-the-analyses.md) explains how to run
individual steps or generate fresh simulations. The
[table and figure index](verification/reproduction-map.csv) links each paper
result to its script and output file.

## Participant data

The participant-level OCRS file is not redistributed. The supplied aggregate
inputs are sufficient for the default run, but cannot reconstruct participant
eligibility or the first-stage regressions. For those steps, obtain the original
file and run `./run_all.sh --empirical` as described in
[Data access](data/DATA_ACCESS.md). This also repeats the participant bootstrap
and subset analyses.

Empirical posterior summaries condition on the supplied sampling variances.
The [running guide](docs/running-the-analyses.md#interpretation) describes the
scope of the simulation comparisons and empirical uncertainty summaries.

## Files

| Directory | Contents |
|---|---|
| `analysis/inputs/` | Saved simulation losses, aggregate data, and plotting inputs |
| `analysis/scripts/` | Analysis scripts and the original simulation routines |
| `analysis/outputs/` | Recomputed tables and numerical checks |
| `reference/` | Paper results used for comparison |
| `fig_src/`, `figures/` | Figure sources and generated PDFs |
| `manifest/`, `verification/` | Input checksums and result comparison code |

Keep the supplied inputs when rerunning the analyses. Fresh simulations write
to new directories under `analysis/generated/`.

## Citation and license

Lee, J. (2026). *The Geometry of Informativeness: A Bregman-Centroid Diagnostic
for Partial Pooling in Multisite Trials*. Citation metadata are in
[`CITATION.cff`](CITATION.cff).

Code is available under the [MIT license](LICENSE). Bundled fonts and simulation
routines retain their own license notices; OCRS data retain their original
release terms.

JoonHo Lee, The University of Alabama — [jlee296@ua.edu](mailto:jlee296@ua.edu).
