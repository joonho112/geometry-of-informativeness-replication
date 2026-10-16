# Recreate the fixed-N site sizes and compare their diagnostics with saved values.
# Run from analysis/. This does not refit the hierarchical models.

source('scripts/legacy/00_utils.R')
source('scripts/dependencies/historical-dgp.R')
dir.create('outputs/evidence-audit', recursive = TRUE, showWarnings = FALSE)
SIGMA_T <- .2
REPS <- 500L
# Read the original effect and site-size generators without running the fits.
ex <- parse('scripts/legacy/08_sim_fixedN.R')
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) &&
      is.symbol(e[[2]]) &&
      as.character(e[[2]]) %in% c('gen_tau_matrix', 'draw_sizes')) eval(e)
}
fam <- expand.grid(
  J = c(25L, 50L, 100L),
  nbar = c(20, 50, 100),
  G = c('Gaussian', 'Mixture', 'ALD'),
  stringsAsFactors = FALSE
)
rows <- list()
k <- 0
for (i in seq_len(nrow(fam))) {
  g <- fam[i, ]
  set.seed(20260731L + i)
  # Consume the same latent draws as the original run before drawing site sizes.
  z <- gen_tau_matrix(g$G, g$J, REPS)
  eps <- matrix(rnorm(REPS * g$J), REPS, g$J)
  # cv controls the parent Gamma law; clipping and rescaling change the realized CV.
  for (cv in c(0, .25, .5, .75)) {
    n <- draw_sizes(g$J, g$nbar, cv)
    dg <- diagnostics_profile(.04, 4 / n)
    k <- k + 1
    rows[[k]] <- data.frame(
      family = i,
      J = g$J,
      nbar = g$nbar,
      G = g$G,
      CV = cv,
      total = sum(n),
      min_size = min(n),
      max_size = max(n),
      realized_CV = sqrt(mean((n - mean(n))^2)) / mean(n),
      I = dg$I,
      V_logit = dg$V
    )
  }
}
z <- do.call(rbind, rows)
old <- unique(read.csv('inputs/sim_fixedN_cells.csv')[c('J', 'nbar', 'G', 'CV', 'I', 'V_logit')])
m <- merge(z, old, by = c('J', 'nbar', 'G', 'CV'), suffixes = c('', '.saved'))
stopifnot(
  nrow(m) == 108,
  max(abs(m$I - m$I.saved)) < 1e-12,
  max(abs(m$V_logit - m$V_logit.saved)) < 1e-12,
  max(abs(m$total - m$J * m$nbar)) < 1e-10
)
write.csv(z, 'outputs/evidence-audit/fixedN-realized-designs.csv', row.names = FALSE)
jsonlite::write_json(
  list(
    design_diagnostic_reproduction = 'PASS',
    min_actual_size = min(z$min_size),
    designs_below_5 = sum(z$min_size < 5 - 1e-12),
    target_total_maxabs = max(abs(z$total - z$J * z$nbar)),
    note = 'Parent Gamma CV is a tuning parameter; clipped and twice rescaled sizes are continuous, not integer enrollments; last rescaling may fall below 5.'
  ),
  'outputs/evidence-audit/fixedN-design-summary.json',
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
print(range(z$min_size))
print(table(z$CV, z$min_size < 5 - 1e-12))
cat('Effect generators: multisiteDGP 0.1.0 numeric kernels\n')
