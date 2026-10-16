# Check simulation summaries against the saved losses and tables.
# Run from analysis/ after 05_simulation_postprocess.R; no models are refitted.

source('scripts/legacy/00_utils.R')
suppressPackageStartupMessages(library(jsonlite))
dir.create('outputs/evidence-audit', recursive = TRUE, showWarnings = FALSE)
writeout <- function(x, n) write.csv(
  x,
  file.path('outputs/evidence-audit', paste0(n, '.csv')),
  row.names = FALSE
)
r <- readRDS('inputs/sim_main_reps.rds')
f <- readRDS('inputs/sim_fixedN_reps.rds')
agg <- function(x) aggregate(x[c('rmse', 'mselp', 'isel')], x[c('cell', 'est')], mean)
cells <- agg(r)
fcells <- agg(f)
compare_cells <- function(a, path) {
  b <- read.csv(path)
  m <- merge(a, b, by = c('cell', 'est'), suffixes = c('.recalc', '.saved'))
  stopifnot(nrow(m) == nrow(a), !anyDuplicated(m[c('cell', 'est')]))
  dif <- sapply(
    c('rmse', 'mselp', 'isel'),
    function(g) max(abs(m[[paste0(g, '.recalc')]] - m[[paste0(g, '.saved')]]))
  )
  stopifnot(max(dif) < 1e-12)
  dif
}
cellcheck <- list(
  stress = compare_cells(cells, 'inputs/sim_main_cells.csv'),
  fixedN = compare_cells(fcells, 'inputs/sim_fixedN_cells.csv')
)
# Compare table values. The paired MCSE script updates the SE columns separately.
tabcheck <- list()
for (n in c(
  'tbl_fair',
  'tbl_metamodel',
  'tbl_mcse',
  'tbl_ratios',
  'tbl_ratios_shape',
  'tbl_oracle'
)) {
  a <- read.csv(paste0('inputs/tables/', n, '.csv'))
  b <- read.csv(paste0('outputs/tables/', n, '.csv'))
  cols <- intersect(names(a), names(b))
  if (n %in% c('tbl_ratios', 'tbl_ratios_shape', 'tbl_oracle')) cols <- cols[!grepl('^se_', cols)]
  stopifnot(nrow(a) == nrow(b))
  errors <- sapply(
    cols,
    function(z) if (is.numeric(a[[z]])) max(abs(a[[z]] - b[[z]])) else as.numeric(!identical(a[[z]], b[[z]]))
  )
  stopifnot(max(errors) < 1e-10)
  tabcheck[[n]] <- errors
}
# Geometric loss ratios by effect distribution, averaging over J.
b <- subset(r, V == 0, c(family, rep, est, rmse, mselp, isel))
names(b)[4:6] <- paste0(names(b)[4:6], '0')
m <- merge(subset(r, V == 2 & shape == 'even'), b, by = c('family', 'rep', 'est'))
for (g in c('rmse', 'mselp', 'isel')) m[[g]] <- log(m[[g]] / m[[paste0(g, '0')]])
rob <- aggregate(m[c('rmse', 'mselp', 'isel')], m[c('est', 'I', 'G')], function(x) exp(mean(x)))
writeout(rob, 'stress-ratios-by-G')
# Check positive losses before calculating log ratios and rank improvements.
stopifnot(all(r$rmse > 0), all(r$mselp > 0), all(r$isel > 0))
g <- merge(
  subset(r, shape == 'even' & V > 0 & est == 'gr'),
  subset(r, shape == 'even' & V > 0 & est == 'ml', c(family, combo, rep, mselp)),
  by = c('family', 'combo', 'rep'),
  suffixes = c('', '.ml')
)
g$gainlog <- log(g$mselp / g$mselp.ml)
gains <- aggregate(g['gainlog'], g[c('I', 'V')], function(x) 1 - exp(mean(x)))
writeout(gains, 'gr-rank-gains-even')
# Best observed mean loss in each cell; this is a Monte Carlo comparison.
wins <- do.call(rbind, lapply(c('rmse', 'mselp', 'isel'), function(goal) {
  z <- subset(cells, est %in% c('ml', 'pm', 'cb', 'gr'))
  d <- do.call(
    rbind,
    lapply(
      split(z, z$cell),
      function(q) data.frame(
        cell = q$cell[1],
        goal,
        winner = q$est[which.min(q[[goal]])],
        pm_is_min = q[[goal]][q$est == 'pm'] <= min(q[[goal]]) + 1e-12
      )
    )
  )
  d
}))
writeout(wins, 'stress-cell-winners')
# Fixed-N ratios: divide cell mean losses, then average across families.
fmeta <- unique(f[c(
  'cell',
  'family',
  'est',
  'J',
  'nbar',
  'G',
  'CV',
  'I',
  'V_logit',
  'H',
  'S_bar'
)])
fz <- merge(fcells, fmeta, by = c('cell', 'est'))
fz0 <- subset(fz, CV == 0, c(family, est, rmse, mselp, isel))
names(fz0)[3:5] <- paste0(names(fz0)[3:5], '0')
fm <- merge(subset(fz, CV > 0), fz0, by = c('family', 'est'))
for (g in c('rmse', 'mselp', 'isel')) fm[[paste0(g, '_ratio')]] <- fm[[g]] / fm[[paste0(g, '0')]]
writeout(fm, 'fixedN-family-cellmean-ratios')
fs <- aggregate(fm[c('rmse_ratio', 'mselp_ratio', 'isel_ratio')], fm[c('est', 'CV')], mean)
writeout(fs, 'fixedN-mean-cell-ratios')
fd <- aggregate(fz[c('I', 'V_logit')], fz['CV'], mean)
writeout(fd, 'fixedN-diagnostics')
# Check whether observed mean losses rise with CV in each family.
fixclaims <- do.call(rbind, lapply(split(fz, list(fz$family, fz$est), drop = TRUE), function(z) {
  z <- z[order(z$CV), ]
  data.frame(
    family = z$family[1],
    est = z$est[1],
    goal = c('rmse', 'mselp', 'isel'),
    monotone = sapply(z[c('rmse', 'mselp', 'isel')], function(x) all(diff(x) >= 0)),
    any_below_balanced = sapply(z[c('rmse', 'mselp', 'isel')], function(x) any(x[-1] < x[1]))
  )
}))
writeout(fixclaims, 'fixedN-monotonicity')
# Compare fitted-PM loss changes with the affine-risk direction.
# The exact affine-risk identity does not apply to fitted PM.
mm <- merge(subset(r, V > 0 & est == 'pm'), subset(b, est == 'pm'), by = c('family', 'rep', 'est'))
mm$lr <- log(mm$rmse / mm$rmse0)
dirs <- aggregate(mm['lr'], mm[c('family', 'I', 'V', 'shape', 'S_bar')], mean)
dirs$agrees_affine_risk <- sign(dirs$lr) == sign(dirs$I - dirs$S_bar)
writeout(dirs, 'pm-family-direction')
# Total-precision multiplier for each sampling-variance profile.
prof <- unique(r[c('J', 'I', 'V', 'shape')])
prof$precision_multiplier <- mapply(function(J, I, V, shape) {
  gm <- .04 * (1 - I) / I
  v <- make_profile(shape, J, gm, V)
  mean(gm / v)
}, prof$J, prof$I, prof$V, prof$shape)
writeout(prof, 'profile-precision-multipliers')
summary <- list(
  cell_reproduction_maxabs = cellcheck,
  table_reproduction_maxabs = tabcheck,
  stress_rows = nrow(r),
  stress_cells = length(unique(r$cell)),
  fixed_rows = nrow(f),
  fixed_cells = length(unique(f$cell)),
  pm_RMSE_best_cells = sum(subset(wins, goal == 'rmse')$pm_is_min),
  total_stress_cells = 702L,
  fixedN_nonmonotone_rows = sum(!fixclaims$monotone),
  fixedN_any_below_balanced = sum(fixclaims$any_below_balanced),
  fixedN_checks = nrow(fixclaims),
  pm_profile_direction_disagreements = sum(!dirs$agrees_affine_risk),
  pm_profile_direction_total = nrow(dirs),
  precision_multiplier_V2 = range(subset(prof, V == 2)$precision_multiplier),
  GR_rank_gain_V2 = range(subset(gains, V == 2)$gainlog),
  scope = 'Recomputed existing saved losses and output tables; no full simulation rerun, no population robustness or coverage claim.'
)
write_json(
  summary,
  'outputs/evidence-audit/summary.json',
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
writeLines(capture.output(sessionInfo()), 'outputs/evidence-audit/sessionInfo.txt')
print(summary)
