#!/usr/bin/env Rscript
# Check empirical outputs after running --reproduce and --empirical.
# Usage: Rscript verification/empirical-verify.R [reference_directory]
# The optional directory contains numerical-stability/ and ocrs-audit/ CSV files.
# If numerical-stability/ contains ocrs_fits_corrected.rds, compare its retained fits.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop('Supply at most one reference directory')
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
root <- dirname(dirname(script))
setwd(root)
stopifnot(requireNamespace('jsonlite', quietly = TRUE), requireNamespace('digest', quietly = TRUE))
numerics <- 'analysis/outputs/numerical-stability'
empirical <- 'analysis/outputs/ocrs-audit'
read <- function(folder, name) read.csv(file.path(folder, name), check.names = FALSE)
fit <- readRDS(file.path(numerics, 'ocrs_fits_corrected.rds'))
summ <- read(numerics, 'ocrs-corrected-summaries.csv')
allfits <- c(fit$fits, list(logodds = fit$logodds_fit))
stopifnot(length(allfits) == 6L, nrow(summ) == 168L)
for (nm in names(allfits)) {
  z <- allfits[[nm]]
  rr <- summ[summ$fit == nm, ]
  stopifnot(
    all(z$v > 0), length(z$v) == length(z$y),
    max(abs(unname(z$reference_summary[rr$metric]) - rr$value)) < 1e-12,
    !'theta' %in% names(z$fit), is.null(z$pdg$draws),
    length(z$pdg$plot_draws$I) == 20000L
  )
}

stopifnot(
  abs(fit$fits$pooled$pdg$I['mean'] - .670895621552224) < 1e-12,
  abs(fit$fits$credits$pdg$I['mean'] - .318137143980031) < 1e-12
)
input <- read(empirical, 'saved-input-comparison.csv')
mat <- read(empirical, 'design-matrix-audit.csv')
stopifnot(
  nrow(input) == 6L, nrow(mat) == 96L,
  max(input$max_y_diff) < 1e-10, max(input$max_v_diff) < 1e-10,
  max(mat$independent_coef_absdiff) < 1e-9,
  max(mat$independent_v_absdiff) < 1e-9
)
boot <- readRDS(file.path(empirical, 'bootstrap-replication.rds'))
sub <- readRDS(file.path(empirical, 'subset-replication.rds'))
stopifnot(
  identical(boot$participant_data_included, FALSE),
  nrow(boot$passing$y) == 500L, nrow(boot$credits$y) == 500L,
  nrow(sub$keep) == 256L,
  nrow(read(empirical, 'bootstrap-reference-passing.csv')) == 500L,
  nrow(read(empirical, 'bootstrap-reference-credits.csv')) == 500L,
  nrow(read(empirical, 'subset-reference-replicates.csv')) == 256L,
  nrow(read(empirical, 'bootstrap-floor-reference-replicates.csv')) == 3000L
)

# Check saved objects for functions, personal paths, and participant fields.
scan <- function(x) {
  if (is.environment(x) || is.function(x)) stop('Environment or function in exported object')
  if (is.character(x) && any(grepl('/Users/|/home/', x))) stop('Private absolute path in exported object')
  if (any(names(x) %in% c('stdids', 'ragrp', 'scm_passcrs', 'rablock')))
    stop('Participant-level schema field in exported object')
  if (is.list(x)) invisible(lapply(x, scan))
  invisible(TRUE)
}

rds <- c(
  'analysis/inputs/ocrs_fits.rds',
  list.files(numerics, '\\.rds$', full.names = TRUE),
  list.files(empirical, '\\.rds$', full.names = TRUE)
)
for (p in rds) scan(readRDS(p))
stopifnot(length(list.files('analysis',
  pattern = 'OCRS_publicuse_impactanalysisdata\\.csv$',
  recursive = TRUE
)) == 0L)

# Compare all CSV files in the supplied reference directory.
comparisons <- list()
compact_checks <- list()
if (length(args)) {
  source <- normalizePath(args[1], mustWork = TRUE)
  for (folder in c('numerical-stability', 'ocrs-audit')) {
    files <- list.files(file.path(source, folder), '\\.csv$')
    if (!length(files)) stop('No reference CSV files in ', file.path(source, folder))
    for (f in files) {
      a <- read(file.path(source, folder), f)
      path <- file.path('analysis/outputs', folder, f)
      stopifnot(file.exists(path))
      z <- read.csv(path, check.names = FALSE)
      stopifnot(identical(names(a), names(z)), nrow(a) == nrow(z))
      numeric <- vapply(a, is.numeric, logical(1))
      mx <- 0
      if (any(numeric)) {
        stopifnot(identical(is.na(a[, numeric, drop = FALSE]), is.na(z[, numeric, drop = FALSE])))
        dif <- abs(as.matrix(a[, numeric, drop = FALSE]) - as.matrix(z[, numeric, drop = FALSE]))
        mx <- max(c(0, dif[is.finite(dif)]))
        stopifnot(mx < 1e-9)
      }
      if (any(!numeric)) stopifnot(identical(a[, !numeric, drop = FALSE], z[, !numeric, drop = FALSE]))
      comparisons[[length(comparisons) + 1L]] <- list(
        path = path,
        rows = nrow(z),
        max_numeric_abs_difference = mx,
        pass = TRUE
      )
    }
  }
  fit_reference <- file.path(source, 'numerical-stability/ocrs_fits_corrected.rds')
  if (file.exists(fit_reference)) {
    original <- readRDS(fit_reference)
    compact <- function(z) {
      z$fit <- z$fit[c('pm', 'pv')]
      z$pdg$draws <- NULL
      z$pdg$plot_draws <- z$pdg$plot_draws['I']
      z
    }
    oldfits <- c(original$fits, list(logodds = original$logodds_fit))
    for (nm in names(allfits)) {
      stopifnot(isTRUE(all.equal(compact(oldfits[[nm]]), allfits[[nm]], tolerance = 1e-12)))
      compact_checks[[nm]] <- TRUE
    }
  }
}

report <- list(
  status = 'pass', checked_utc = format(Sys.time(), tz = 'UTC', usetz = TRUE),
  six_fits = TRUE, first_stage_saved_inputs = TRUE, independent_HC2 = TRUE,
  bootstrap_refits_per_outcome = 500L, subsets = 256L, floor_cases = 3000L,
  rds_structural_privacy_checks = length(rds), participant_csv_distributed = FALSE,
  reference_numeric_comparisons = comparisons, compact_fit_equivalence = compact_checks
)
dir.create('verification/generated', showWarnings = FALSE, recursive = TRUE)
jsonlite::write_json(report,
  'verification/generated/empirical-verification.json',
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
cat(
  'Empirical verification passed:',
  length(allfits),
  'fits;',
  length(comparisons),
  'CSV comparisons;',
  length(rds),
  'aggregate-only RDS files.\n'
)
