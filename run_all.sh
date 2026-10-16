#!/usr/bin/env bash
# Reproduce the paper's results using the supplied inputs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"
mode="${1:---reproduce}"
if [ "$#" -gt 1 ]; then echo 'Choose one mode; use --help.' >&2; exit 2; fi
help_text() {
  cat <<'HELP'
Usage: ./run_all.sh [MODE]
  --verify            Check the supplied input files (Python standard library).
  --theory, --quick   Run the numerical checks in a separate working directory.
  --simulation       Recompute saved-loss summaries and final component MCSEs.
  --posterior        Recompute six empirical posterior specifications and numerical checks.
  --sensitivity      Recompute quadrature summaries of saved aggregate resamples.
  --figures          Draw Figures 1–8 from the supplied results.
  --figure-qa        Render and check the generated PDFs (Poppler + Python packages).
  --empirical        Rebuild first stage and all resampling from external OCRS_DATA.
  --smoke            Small fresh simulation runs; not the manuscript simulation results.
  --full-simulation  Rerun both full historical simulation arms in separate directories.
  --reproduce        Default: verify, theory, saved simulation, posterior, saved sensitivity,
                     numerical comparisons, and figures. No participant data needed.
  --help             Show this help.
Set PYTHON and RSCRIPT to choose runtimes. See README.md for dependencies and scope.
Results are written to analysis/outputs/ and figures/. Input files are kept unchanged.
HELP
}
case "$mode" in
  --help|-h) help_text; exit 0 ;;
  --verify|--theory|--quick|--simulation|--posterior|--sensitivity|--figures|--figure-qa|--empirical|--smoke|--full-simulation|--reproduce) ;;
  *) echo "Unknown mode: $mode" >&2; help_text >&2; exit 2 ;;
esac
verify() { "$PYTHON" manifest/verify.py; }
compare() { "$PYTHON" verification/compare_outputs.py --group "$1"; }
theory() { "$RSCRIPT" analysis/scripts/run_historical.R --arm=theory; }
simulation() {
  (cd analysis
   "$RSCRIPT" scripts/05_simulation_postprocess.R
   "$RSCRIPT" scripts/05_evidence_audit.R
   "$RSCRIPT" scripts/05_fixedN_design_audit.R
   "$RSCRIPT" scripts/05_paired_component_mcse.R)
  compare simulation
}
posterior() { "$RSCRIPT" analysis/scripts/05_posterior_stability.R; compare posterior; }
sensitivity() { "$RSCRIPT" analysis/scripts/05_ocrs_audit.R --aggregate-reference; compare sensitivity; }
figures() { "$RSCRIPT" analysis/scripts/07_figures.R; }
verify
case "$mode" in
  --verify) ;;
  --theory|--quick) theory ;;
  --simulation) simulation ;;
  --posterior) posterior ;;
  --sensitivity) sensitivity ;;
  --figures) figures ;;
  --figure-qa) "$PYTHON" analysis/scripts/07_figure_qa.py ;;
  --empirical)
    "$RSCRIPT" analysis/scripts/05_ocrs_audit.R --rerun-sensitivity --reference-comparisons --final-checks
    compare firststage; compare sensitivity ;;
  --smoke)
    "$RSCRIPT" analysis/scripts/run_historical.R --arm=stress --smoke
    "$RSCRIPT" analysis/scripts/run_historical.R --arm=fixedN --smoke ;;
  --full-simulation)
    "$RSCRIPT" analysis/scripts/run_historical.R --arm=stress
    "$RSCRIPT" analysis/scripts/run_historical.R --arm=fixedN ;;
  --reproduce) theory; simulation; posterior; sensitivity; figures ;;
esac
printf '\nCompleted %s. See analysis/outputs/, figures/, and verification/generated/.\n' "$mode"
