# Analysis scripts

Run the analyses with `./run_all.sh` from the repository directory; see the
[README](../README.md) and [running guide](../docs/running-the-analyses.md).

- `inputs/` contains the supplied simulation losses and aggregate empirical data.
- `scripts/` contains the analysis code; `scripts/legacy/` contains the original
  simulation and estimation routines.
- `outputs/` contains regenerated tables and numerical checks.
- `generated/` holds fresh simulation runs, each in a separate directory.

The four saved-loss simulation scripts run from this directory in the order
shown in the running guide. Run `05_paired_component_mcse.R` last because it
replaces the preliminary standard errors. The driver handles this order.
Empirical and figure scripts locate their inputs relative to the script file.
Main figures are written to `../figures/`.
