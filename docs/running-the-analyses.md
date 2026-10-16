# Running the analyses

Run commands from the repository directory. `./run_all.sh --reproduce` runs the
numerical checks, saved-loss simulation analysis, empirical posterior and
sensitivity calculations, and figures. Each step can also be run separately.

## Setup

The [README](../README.md#requirements) lists dependencies. To use the versions
recorded in `renv.lock` in a local R library:

```sh
Rscript -e 'install.packages("renv", repos="https://cloud.r-project.org")'
Rscript -e 'renv::restore(lockfile="renv.lock", library=".R-library", prompt=FALSE)'
export R_LIBS_USER="$PWD/.R-library"
```

Set `RSCRIPT` or `PYTHON` to select a different executable. The scripts were
tested on macOS. Windows users need a Bash environment for `run_all.sh` or can
run the individual R and Python scripts. Fresh simulations support one worker
on Windows.

## Individual steps

| Command | Calculation |
|---|---|
| `./run_all.sh --verify` | Check supplied input files against their SHA-256 checksums |
| `./run_all.sh --theory` | Numerical checks of the formulas, including the effective site count |
| `./run_all.sh --simulation` | Summarize saved losses, fit metamodels, and calculate Monte Carlo standard errors (MCSEs) |
| `./run_all.sh --posterior` | Calculate six OCRS posterior specifications by adaptive quadrature |
| `./run_all.sh --sensitivity` | Refit saved aggregate bootstrap and subset inputs by quadrature |
| `./run_all.sh --figures` | Draw Figures 1–8 from saved plotting inputs |
| `./run_all.sh --empirical` | Rebuild first-stage estimates and resamples from external participant data |
| `./run_all.sh --smoke` | Run small fresh simulations to check the code |
| `./run_all.sh --full-simulation` | Rerun both complete simulation experiments |

`--quick` is an alias for `--theory`. The numerical formula checks supplement
the proofs; they do not verify them. For the participant data, see
[Data access](../data/DATA_ACCESS.md).

The saved-loss analysis runs four scripts from `analysis/`, in this order:

```sh
cd analysis
Rscript scripts/05_simulation_postprocess.R
Rscript scripts/05_evidence_audit.R
Rscript scripts/05_fixedN_design_audit.R
Rscript scripts/05_paired_component_mcse.R
cd ..
```

The final script replaces preliminary standard errors with MCSEs that account
for shared random seeds. Use the driver to avoid accidentally stopping before
this step. It also compares the results with those in `reference/`.

## Fresh simulations

The stress experiment has 702 cells, 500 replications per cell, and six
procedures. The fixed-enrollment experiment has 108 cells, 500 replications per
cell, and five procedures. The default reproduction uses the included losses
from these experiments.

To rerun one complete experiment:

```sh
Rscript analysis/scripts/run_historical.R --arm=stress --cores=1
Rscript analysis/scripts/run_historical.R --arm=fixedN --cores=1
```

Each run creates a new directory under `analysis/generated/`. An optional
`--workdir=/new/path` must name a directory that does not yet exist. Full runs
use 500 replications and 2,000 posterior draws; they can take considerable time.
They use the original 250-node grid/jitter posterior algorithm and the included
simulation routines. A current `multisiteDGP` installation is not required.
New results do not replace the supplied inputs automatically.

Add `--smoke` to either command for three families, two profiles per family,
two replications, and 200 posterior draws. These small runs check execution,
not the paper's estimates. `--spot-check` instead retains the original fitting
settings for selected cases and compares their losses with the saved results.

## Interpretation

Simulation metamodel predictors `I`, `H`, `V`, and `S_bar` come from the true
generating reliability profile. Five-fold prediction is evaluated within the
chosen cell grid, with related cells in the training folds; it does not assess
prediction in new trials.

Stress summaries are geometric means of paired loss ratios. Fixed-enrollment
summaries average ratios of arithmetic cell means. The reported K=9
component-MCSE bounds do not apply to the K=3 summaries conditional on the
effect distribution or to the fixed-enrollment aggregate ratios.

The simulation evaluates the implemented PM, CB, and GR procedures. Its grid
fits are not replaced by the quadrature calculations used for the OCRS
analysis. The six synthetic numerical examples check PM and implemented CB;
they do not establish accuracy for every simulation cell or for GR.

OCRS inference uses a normal working model with a reference measure flat in
the mean and standard deviation, conditional on the supplied sampling
variances. Posterior means and intervals come from adaptive quadrature;
cached draws are used for displays and numerical comparisons. Variation across
bootstrap refits is distinct from conditional posterior uncertainty. Neither
summary establishes latent-effect coverage. The shift in the credits bootstrap
and dependence due to shared students remain limitations of the analysis.

## Figures and output comparisons

The default run compares regenerated results with saved paper results before
drawing the figures. The figures read saved inputs in
`analysis/inputs/figure-reference/`, not the latest analysis outputs. Thus
`--figures` reproduces the paper's displays without refitting models.

To check the generated PDFs, install the optional Python packages and run:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-qa.txt
PYTHON="$PWD/.venv/bin/python" ./run_all.sh --figure-qa
```

This checks fonts, dimensions, and text boundaries and creates previews in
`figures/png/`. See the [figure notes](../fig_src/README-current-figures.md) for
source files. Numerical comparisons use absolute tolerance 1e-10 and relative
tolerance 1e-8 unless a specific check states otherwise. PDF checksums may change
between builds because of metadata.

Saved tables in `analysis/inputs/tables/` are comparison inputs; reported
tables are regenerated in `analysis/outputs/`.

After running `--reproduce` and `--empirical`, check the empirical fits,
first-stage estimates, resamples, and saved results with:

```sh
Rscript verification/empirical-verify.R reference
```
