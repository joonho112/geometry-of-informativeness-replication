#!/usr/bin/env Rscript
# Draw Figures 1–8 from the results supplied with the paper.
# Run: Rscript analysis/scripts/07_figures.R
# Figure 1 is drawn in TikZ; the others use R with TeX labels.
# All figures use Source Sans 3 text and Computer Modern mathematics.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})
arg <- commandArgs(FALSE)
this <- sub('^--file=', '', arg[grepl('^--file=', arg)])
ROOT <- normalizePath(file.path(dirname(this), '../..'))
setwd(ROOT)
source('fig_src/cm_math_device.R')
OUT <- 'analysis/outputs/figure-audit'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create('figures/png', recursive = TRUE, showWarnings = FALSE)
source('analysis/scripts/legacy/00_utils.R')
assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)
close_to <- function(x, y, tol = 1e-10) all(is.finite(x)) && all(is.finite(y)) && max(abs(x - y)) < tol
sha <- function(f) digest::digest(file = f, algo = 'sha256')
inputs <- c(
  'manifest/figure-contract.json',
  'analysis/scripts/legacy/00_utils.R',
  'fig_src/fonts/SourceSans3-Regular.otf',
  'fig_src/fonts/SourceSans3-Bold.otf',
  'analysis/inputs/sim_main_reps.rds',
  'analysis/inputs/precomputed/sim_landscape.rds',
  'analysis/inputs/precomputed/sim_stress_cells.rds',
  'analysis/inputs/tables/tbl_ocrs_sites.csv',
  'analysis/inputs/figure-reference/tbl_ratios_shape.csv',
  'analysis/inputs/figure-reference/tbl_fair.csv',
  'analysis/inputs/figure-reference/design-illustration-data.csv',
  'analysis/inputs/figure-reference/ocrs_fits_corrected.rds',
  'analysis/inputs/figure-reference/ocrs-corrected-summaries.csv'
)
source_hashes <- setNames(vapply(inputs, sha, ''), inputs)
# Check that the input files match those used for the paper.
lock <- 'manifest/figure-source-lock.json'
assert(file.exists(lock), 'Missing figure input checksums')
old <- unlist(jsonlite::read_json(lock))
assert(
  identical(source_hashes, old),
  'An input checksum changed. Check the file before updating the recorded checksum.'
)
COL <- c(
  blue = '#0072B2', sky = '#56B4E9', orange = '#E69F00', vermillion = '#D55E00',
  green = '#009E73', purple = '#CC79A7', ink = '#1A1A1A', bar = '#7FA8C9', red = '#B2182B'
)
FAMILY <- 'Source Sans 3'
PT <- 12
SZ <- PT / ggplot2::.pt
theme_paper <- function() theme_minimal(base_size = PT, base_family = FAMILY) +
  theme(
    text = element_text(size = PT, color = unname(COL['ink'])),
    axis.text = element_text(size = PT, color = unname(COL['ink'])), axis.title = element_text(size = PT),
    strip.text = element_text(size = PT, face = 'bold', hjust = 0),
    plot.title = element_text(size = PT, face = 'bold', hjust = 0),
    plot.subtitle = element_text(size = PT), plot.caption = element_text(size = PT),
    legend.text = element_text(size = PT), legend.title = element_text(size = PT),
    panel.grid.minor = element_blank(), panel.grid.major = element_line(color = 'grey90', linewidth = .25),
    plot.margin = margin(6, 8, 6, 6), legend.position = 'bottom', legend.box = 'vertical',
    panel.spacing = unit(9, 'pt')
  )
# Collect PDF dimensions and fonts for the figure checks.
specs <- list()
savefig <- function(p, name, h) {
  assert(h <= 6.3, 'Figure exceeds the maximum height')
  f <- file.path('figures', paste0(name, '.pdf'))
  save_cm_figure(p, name, h)
  specs[[name]] <<- list(
    pdf = f, width_inches = 6.5, height_inches = h, base_font_pt = PT,
    text_family = 'Source Sans 3', math_family = 'Computer Modern',
    authoring_source = paste0('fig_src/generated-cm/', name, '.tex')
  )
  message('Wrote ', f, ' (6.5 x ', h, ' inches)')
}
# Drawing helpers for the workflow diagram; coordinates are in inches.
tt <- function(x, y, label, bold = FALSE, col = unname(COL['ink']), just = 'centre', ...) grid.text(label,
  x = unit(x, 'inches'), y = unit(y, 'inches'), just = just, gp = gpar(
    fontfamily = FAMILY, fontsize = PT,
    fontface = if (bold) 'bold' else 'plain', col = col
  ), ...
)
ln <- function(x0, y0, x1, y1, col = unname(COL['ink']), arrow = FALSE, lty = 1, lwd = 1) grid.lines(
  x = unit(c(x0, x1), 'inches'), y = unit(c(y0, y1), 'inches'),
  gp = gpar(col = col, lty = lty, lwd = lwd), arrow = if (arrow) grid::arrow(length = unit(.08, 'inches'), type = 'closed') else NULL
)
bx <- function(x, y, w, h, fill = 'white', col = 'grey55') grid.roundrect(
  x = unit(x, 'inches'), y = unit(y, 'inches'),
  width = unit(w, 'inches'), height = unit(h, 'inches'), r = unit(.04, 'inches'), gp = gpar(fill = fill, col = col, lwd = .8)
)
# Figure 1: centered L2 projection and the Bernoulli–KL centroid.
if (!('--figures2-8' %in% commandArgs(TRUE))) {
  fig1_status <- system2(FIGURE_PYTHON, shQuote(file.path(ROOT, 'fig_src/build_fig1.py')))
  assert(fig1_status == 0L, 'Native LaTeX Figure 1 build failed')
}
specs[['fig1_geometry']] <- list(
  pdf = 'figures/fig1_geometry.pdf',
  width_inches = 6.5, height_inches = 4.8, base_font_pt = 12,
  authoring_source = 'fig_src/fig1_geometry_cm.tex', text_family = 'Source Sans 3',
  math_family = 'Computer Modern'
)
# Figure 2: three reliability profiles with the same informativeness.
gm <- .2^2 * (1 - .3) / .3
vs <- list(
  'Homogeneous (V = 0)' = profile_even(40, gm, 0),
  'Even spread (V = 1.5)' = profile_even(40, gm, 1.5),
  'Concentrated precise (V = 1.5)' = profile_concentrated(40, gm, 1.5)
)
d2 <- do.call(rbind, lapply(names(vs), function(nm) {
  d <- diagnostics_profile(.2^2, vs[[nm]])
  assert(close_to(d$I, .3), 'Figure 2: I target failure')
  data.frame(profile = nm, site = 1:40, S = sort(d$S), I = d$I, Sbar = d$S_bar, H = d$H, V = d$V, Jeff = d$J_eff)
}))
d2$profile <- factor(d2$profile, levels = names(vs))
l2 <- unique(d2[c('profile', 'I', 'Sbar', 'H', 'V', 'Jeff')])
l2$label <- sprintf('$H=%.3f\\quad\\bar S-I=%+.3f\\quad J_{\\mathrm{eff}}=%.1f$', pmax(l2$H, 0), l2$Sbar - l2$I, l2$Jeff)
write.csv(d2, file.path(OUT, 'fig2-profile-values.csv'), row.names = FALSE)
p2 <- ggplot(d2, aes(site, S)) +
  geom_col(fill = unname(COL['bar']), width = .85) +
  geom_hline(aes(yintercept = I, linetype = 'I'), color = unname(COL['ink']), linewidth = .6) +
  geom_hline(aes(yintercept = Sbar, linetype = 'Mean reliability'), color = unname(COL['vermillion']), linewidth = .6) +
  geom_text(data = l2, aes(x = 1, y = 1.20, label = label), inherit.aes = FALSE, hjust = 0, size = SZ, family = FAMILY, parse = FALSE) +
  facet_wrap(~profile, ncol = 1, labeller = as_labeller(c(
    'Homogeneous (V = 0)' = 'Homogeneous ($V=0$)',
    'Even spread (V = 1.5)' = 'Even spread ($V=1.5$)',
    'Concentrated precise (V = 1.5)' = 'Concentrated precise ($V=1.5$)'
  ))) +
  scale_y_continuous(limits = c(0, 1.28), breaks = c(0, .5, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40), expand = expansion(mult = c(.01, .01))) +
  scale_linetype_manual(
    values = c(I = 'solid', 'Mean reliability' = '22'),
    labels = c(I = '$I$', 'Mean reliability' = 'Mean reliability'), name = NULL
  ) +
  labs(x = 'Site, sorted by reliability', y = '$S_j$') +
  theme_paper() +
  theme(legend.margin = margin(0, 0, 0, 0), legend.key.width = unit(24, 'pt'))
savefig(p2, 'fig2_profiles', 6.2)
# Check the saved posterior summaries against the table values before plotting.
fits <- readRDS('analysis/inputs/figure-reference/ocrs_fits_corrected.rds')$fits
ref <- read.csv('analysis/inputs/figure-reference/ocrs-corrected-summaries.csv')
for (nm in names(fits)) for (metric in names(fits[[nm]]$reference_summary)) {
  z <- ref$value[ref$fit == nm & ref$metric == metric]
  assert(length(z) == 1 && close_to(z, unname(fits[[nm]]$reference_summary[metric])), paste('Empirical source mismatch', nm, metric))
}
po <- fits$pooled
cr <- fits$credits
assert(close_to(po$pdg$I['mean'], .670895621552224), 'Pooled I anchor changed')
assert(close_to(cr$pdg$I['mean'], .318137143980031, 1e-9), 'Credits I anchor changed')
fmt <- function(x, k = 3) sprintf(paste0('%.', k, 'f'), x)
ival <- function(x, k = 3) paste0(fmt(x['mean'], k), '\\;[', fmt(x['2.5%'], k), ', ', fmt(x['97.5%'], k), ']')
# Figure 3: posterior calculations for the OCRS example.
f3 <- function() {
  grid.newpage()
  tt(.16, 6.05, 'Empirical computation', TRUE, just = 'left')
  left <- 1.91
  w <- 3.45
  yy <- c(5.30, 4.15, 3.00, 1.85, .70)
  heads <- c('1. Inputs', '2. Working normal model', '3. Posterior integration', '4. Reliability profile', '5. Report')
  body <- c(
    'Site contrasts and estimated variances',
    'Flat reference prior on mean and SD\nCondition on positive site variances',
    'Integrate the common mean analytically\nAdaptive quadrature over $\\mathrm{SD}\\ge0$',
    '\n$V$ is fixed by the variance profile',
    'Posterior means and 95\\% credible intervals\nCompare summaries for the intended goal'
  )
  for (i in seq_along(yy)) {
    bx(left, yy[i], w, .93, if (i == 4) '#EFF5FA' else 'white')
    tt(left, yy[i] + .25, heads[i], TRUE)
    tt(left, yy[i] - .13, body[i])
    if (i == 4) tt(left, yy[i] - .015, 'Evaluate $I$, $\\bar S$, $H$ and $J_{\\mathrm{eff}}$')
    if (i < 5) ln(left, yy[i] - .48, left, yy[i + 1] + .49, arrow = TRUE, col = 'grey50')
  }
  bx(5.06, 3.00, 2.48, 5.78, fill = '#F7F7F5')
  tt(5.06, 5.62, 'OCRS: course passing', TRUE)
  tt(5.06, 5.27, '32 school-by-subject sites')
  tt(3.98, 4.79, 'Mean [95\\% credible interval]', just = 'left')
  tt(3.98, 4.38, paste0('$I = ', ival(po$pdg$I), '$'), just = 'left')
  tt(3.98, 3.97, paste0('$H = ', ival(po$pdg$H), '$'), just = 'left')
  tt(3.98, 3.56, paste0('$\\bar S = ', fmt(po$pdg$S_bar['mean']), '$'), just = 'left')
  tt(3.98, 3.15, paste0('$J_{\\mathrm{eff}} = ', fmt(po$pdg$J_eff['mean'], 1), ';\\ V = ', fmt(po$pdg$V['mean']), '$'), just = 'left')
  tt(3.98, 2.63, sprintf('$\\Pr(I > 0.5) = %.3f$', po$p_I_gt_half), just = 'left')
  tt(3.98, 2.22, sprintf('Plug-in $I = %.3f$', po$I_plug), just = 'left')
  tt(5.06, 1.51, 'Credits in year 4', TRUE)
  tt(3.98, 1.10, paste0('$I = ', ival(cr$pdg$I), '$'), just = 'left')
  tt(5.06, .48, 'Reported values use\nadaptive quadrature.')
}
savefig(f3, 'fig3_workflow', 6.3)
# Figure 4: the lognormal design relation and its local approximation.
design <- read.csv('analysis/inputs/figure-reference/design-illustration-data.csv')
assert(close_to(design$lognormal_population_ratio, (1 + design$CV_population_or_local_scenario^2)^(-.5)), 'Exact lognormal ratio mismatch')
assert(close_to(design$local_ratio, exp(-design$CV_population_or_local_scenario^2 / 2)), 'Local ratio mismatch')
d4 <- unique(design[c('CV_population_or_local_scenario', 'lognormal_population_ratio', 'local_ratio')])
d4 <- rbind(data.frame(CV = d4[[1]], ratio = d4[[2]], method = 'Exact lognormal'), data.frame(CV = d4[[1]], ratio = d4[[3]], method = 'Local approximation'))
p4a <- ggplot(d4, aes(CV, ratio, color = method, linetype = method)) +
  geom_line(linewidth = .85) +
  scale_color_manual(values = c('Exact lognormal' = unname(COL['blue']), 'Local approximation' = unname(COL['vermillion'])), name = NULL) +
  scale_linetype_manual(values = c('Exact lognormal' = 'solid', 'Local approximation' = '22'), name = NULL) +
  scale_x_continuous(breaks = seq(0, 1.2, .3), expand = expansion(mult = c(.01, .01))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) +
  annotate('text', x = .60, y = .27, label = '$(1+\\mathrm{CV}^2)^{-1/2}$', parse = FALSE, color = unname(COL['blue']), size = SZ, family = FAMILY) +
  annotate('text', x = .60, y = .12, label = '$\\exp(-\\mathrm{CV}^2/2)$', parse = FALSE, color = unname(COL['vermillion']), size = SZ, family = FAMILY) +
  labs(x = NULL, y = 'Geometric / arithmetic mean', title = '(a) Hypothetical site-size imbalance') +
  theme_paper() +
  theme(legend.margin = margin(0, 0, 0, 0))
assert(identical(sort(unique(design$sigma)), c(.1, .2, .3, .4)) && all(design$kappa == 4) && all(design$mean_size == 50), 'Figure 4 design constants changed')
assert(close_to(design$I_lognormal_population, design$sigma^2 / (design$sigma^2 + 4 / (50 * design$lognormal_population_ratio))), 'Figure 4 exact I mismatch')
assert(close_to(design$I_local_approximation, design$sigma^2 / (design$sigma^2 + 4 / (50 * design$local_ratio))), 'Figure 4 local I mismatch')
id4 <- rbind(
  data.frame(CV = design[[1]], I = design$I_lognormal_population, sigma = design$sigma, method = 'Exact lognormal'),
  data.frame(CV = design[[1]], I = design$I_local_approximation, sigma = design$sigma, method = 'Local approximation')
)
id4$sigma <- factor(sprintf('%.1f', id4$sigma), levels = c('0.1', '0.2', '0.3', '0.4'))
p4b <- ggplot(id4, aes(CV, I, color = sigma, linetype = method)) +
  geom_line(linewidth = .75) +
  scale_color_manual(values = c('0.1' = '#7298B5', '0.2' = '#3C7FA8', '0.3' = '#0E587F', '0.4' = '#003047'), name = '$\\sigma$') +
  scale_linetype_manual(values = c('Exact lognormal' = 'solid', 'Local approximation' = '22'), guide = 'none') +
  scale_x_continuous(breaks = seq(0, 1.2, .3), expand = expansion(mult = c(.01, .01))) +
  scale_y_continuous(limits = c(0, .7), breaks = c(0, .2, .4, .6)) +
  labs(
    x = 'Hypothetical coefficient of variation of site sizes', y = 'Informativeness $I$',
    title = '(b) Consequence for informativeness'
  ) +
  theme_paper() +
  theme(legend.margin = margin(0, 0, 0, 0))
p4 <- p4a / p4b + plot_layout(heights = c(1, 1.05))
write.csv(d4, file.path(OUT, 'fig4-design-values.csv'), row.names = FALSE)
write.csv(id4, file.path(OUT, 'fig4-informativeness-values.csv'), row.names = FALSE)
savefig(p4, 'fig4_design', 6.3)
# Figure 5: geometric mean loss ratios across reliability profiles.
r <- readRDS('analysis/inputs/sim_main_reps.rds')
base <- subset(r, V == 0 & est == 'pm', c(family, rep, rmse))
names(base)[3] <- 'rmse0'
assert(!anyDuplicated(base[c('family', 'rep')]), 'Duplicate matched V=0 rows')
shape_names <- c(even = 'Even', bimodal = 'Bimodal', concentrated = 'Concentrated precise', concentrated_noisy = 'Concentrated noisy')
land <- do.call(rbind, lapply(names(shape_names), function(sh) {
  z <- merge(subset(r, V > 0 & shape == sh & est == 'pm'), base, by = c('family', 'rep'))
  do.call(rbind, lapply(split(z, list(z$I, z$V), drop = TRUE), function(x)
    data.frame(I = x$I[1], V = x$V[1], shape = sh, pct = 100 * (exp(mean(log(x$rmse / x$rmse0))) - 1))))
}))
pre <- readRDS('analysis/inputs/precomputed/sim_landscape.rds')
check <- merge(land, pre, by = c('I', 'V', 'shape'))
assert(nrow(check) == 72 && close_to(check$pct.x, check$pct.y), 'Figure 5 precomputed mismatch')
tab <- subset(read.csv('analysis/inputs/figure-reference/tbl_ratios_shape.csv'), est == 'pm')
check <- merge(subset(land, V == 2), tab, by = c('I', 'shape'))
assert(close_to(check$pct, 100 * (check$r_rmse - 1)), 'Figure 5 ratio table mismatch')
write.csv(land, file.path(OUT, 'fig5-landscape-values.csv'), row.names = FALSE)
land$shape <- factor(shape_names[land$shape], levels = shape_names)
land$Ilab <- sprintf('%.2f', land$I)
land$tile <- ifelse(round(land$pct) == 0, '0', sprintf('%+d', round(land$pct)))
assert(max(abs(land$pct)) < 17, 'Figure 5 color scale would clip')
p5 <- ggplot(land, aes(factor(Ilab), factor(V), fill = pct)) +
  geom_tile(color = 'white', linewidth = .6) +
  geom_text(aes(label = tile), size = SZ, family = FAMILY) +
  scale_fill_gradient2(
    low = unname(COL['blue']), mid = 'white', high = unname(COL['red']), midpoint = 0, limits = c(-17, 17),
    name = 'PM RMSE change versus $V = 0$ (\\%)', breaks = c(-15, -10, -5, 0, 5, 10, 15)
  ) +
  facet_wrap(~shape, ncol = 2) +
  coord_cartesian(expand = FALSE) +
  labs(x = 'Informativeness $I$', y = 'Dispersion $V$') +
  theme_paper() +
  theme(panel.grid = element_blank(), legend.title.position = 'top', legend.key.width = unit(23, 'pt'), legend.key.height = unit(9, 'pt')) +
  guides(fill = guide_colorbar(barwidth = unit(4.7, 'inches'), barheight = unit(.12, 'inches'), title.position = 'top'))
savefig(p5, 'fig5_landscape', 5.5)
# Figure 6: held-out prediction error for each generating-profile summary.
fair <- read.csv('analysis/inputs/figure-reference/tbl_fair.csv')
assert(nrow(fair) == 27 && !anyDuplicated(fair[c('est', 'goal', 'summary')]), 'Figure 6 needs 9 x 3 rows')
win <- do.call(rbind, lapply(split(fair, list(fair$goal, fair$est), drop = TRUE), function(x) x[which.max(x$cvR2), ]))
assert(sum(win$summary == 'S_bar') == 8 && any(win$goal == 'isel' & win$est == 'gr' & win$summary == 'H'), 'Figure 6 winner mismatch')
fair$residual_pct <- 100 * (1 - fair$cvR2)
fair$row <- paste(toupper(fair$goal), toupper(fair$est), sep = ' / ')
fair$row <- factor(fair$row, levels = rev(as.vector(t(outer(c('RMSE', 'MSELP', 'ISEL'), c('PM', 'CB', 'GR'), paste, sep = ' / ')))))
fair$summary <- factor(fair$summary, levels = c('H', 'V_logit', 'S_bar'))
write.csv(fair, file.path(OUT, 'fig6-crossvalidation-values.csv'), row.names = FALSE)
p6 <- ggplot(fair, aes(residual_pct, row, color = summary, shape = summary)) +
  geom_line(aes(group = row), color = 'grey75', linewidth = .5) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(H = unname(COL['green']), V_logit = unname(COL['orange']), S_bar = unname(COL['blue'])), name = 'Added summary', labels = c('$H$', '$V$', '$\\bar S$')) +
  scale_shape_manual(values = c(H = 16, V_logit = 17, S_bar = 15), name = 'Added summary', labels = c('$H$', '$V$', '$\\bar S$')) +
  scale_x_continuous(limits = c(0, 1.10), breaks = seq(0, 1, .2), expand = c(0, 0)) +
  labs(x = '$100\\times(1-\\mathrm{CV}\\,R^2)$ (lower is better)', y = 'Loss / estimator') +
  theme_paper() +
  theme(panel.grid.major.y = element_blank(), legend.position = 'top')
savefig(p6, 'fig6_metamodel', 5.2)
# Figure 7: mean losses and Monte Carlo error for selected Gaussian cells.
z <- subset(r, J == 50 & G == 'Gaussian' & shape == 'even' & I %in% c(.15, .45, .75))
agg <- do.call(rbind, lapply(split(z, list(z$est, z$I, z$V), drop = TRUE), function(x) {
  assert(nrow(x) == 500 && length(unique(x$rep)) == 500, 'Figure 7 cell count mismatch')
  data.frame(est = x$est[1], I = x$I[1], V = x$V[1], mselp = mean(x$mselp), mselp_se = sd(x$mselp) / sqrt(500), isel = mean(x$isel), isel_se = sd(x$isel) / sqrt(500))
}))
pre <- readRDS('analysis/inputs/precomputed/sim_stress_cells.rds')
chk <- merge(agg, pre, by = c('est', 'I', 'V'))
for (n in c('mselp', 'mselp_se', 'isel', 'isel_se')) assert(close_to(chk[[paste0(n, '.x')]], chk[[paste0(n, '.y')]]), 'Figure 7 precomputed mismatch')
write.csv(agg, file.path(OUT, 'fig7-cell-values.csv'), row.names = FALSE)
estcols <- c(ML = 'grey45', PM = unname(COL['blue']), CB = unname(COL['sky']), GR = unname(COL['vermillion']))
estshapes <- c(ML = 17, PM = 16, CB = 15, GR = 18)
make7 <- function(goal, ests, title) {
  d <- subset(agg, est %in% ests)
  d$est <- factor(toupper(d$est), levels = names(estcols))
  d$mean <- d[[goal]]
  d$lo <- d$mean - 1.96 * d[[paste0(goal, '_se')]]
  d$hi <- d$mean + 1.96 * d[[paste0(goal, '_se')]]
  ggplot(d, aes(V, mean, color = est, shape = est)) +
    geom_line(linewidth = .65) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .035, linewidth = .45) +
    geom_point(size = 2.0, show.legend = TRUE) +
    facet_wrap(~I, ncol = 3, labeller = labeller(I = function(x) paste0('$I=', x, '$')), scales = 'free_y') +
    scale_color_manual(values = estcols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = estshapes, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = c(0, .5, 1, 2), labels = c('0', '.5', '1', '2')) +
    labs(x = 'Dispersion $V$', y = toupper(goal), title = title) +
    theme_paper() +
    theme(
      legend.position = 'bottom', panel.spacing.x = unit(10, 'pt'),
      axis.title.y = element_text(margin = margin(r = 6))
    )
}
p7a <- make7('mselp', c('ml', 'pm', 'gr'), '(a) Ranking loss') + theme(legend.position = 'none', axis.title.x = element_blank())
p7b <- make7('isel', c('pm', 'cb', 'gr'), '(b) Empirical-distribution loss')
p7 <- p7a / p7b + plot_layout(heights = c(1, 1.18))
savefig(p7, 'fig7_stress', 6.3)
rm(r, z)
gc(verbose = FALSE)
# Figure 8: OCRS site estimates and the posterior distribution of informativeness.
sites <- read.csv('analysis/inputs/tables/tbl_ocrs_sites.csv')
assert(nrow(sites) == length(po$y) && close_to(sites$tau_hat, po$y) && close_to(sites$se^2, po$v), 'OCRS site ordering mismatch')
ca <- data.frame(y = po$y, se = sqrt(po$v), pm = po$fit$pm, subject = sites$subject)
ca <- ca[order(ca$y), ]
ca$ord <- seq_len(nrow(ca))
write.csv(ca, file.path(OUT, 'fig8-site-values.csv'), row.names = FALSE)
subcols <- c('Algebra 1' = unname(COL['orange']), 'English 9' = unname(COL['purple']))
subshapes <- c('Algebra 1' = 16, 'English 9' = 17)
p8a <- ggplot(ca, aes(ord, y, color = subject, shape = subject)) +
  geom_hline(yintercept = 0, color = 'grey60', linewidth = .4) +
  geom_errorbar(aes(ymin = y - 1.96 * se, ymax = y + 1.96 * se), width = 0, linewidth = .45) +
  geom_point(size = 1.7) +
  geom_point(aes(y = pm), color = 'black', shape = 18, size = 1.7) +
  scale_color_manual(values = subcols, name = NULL) +
  scale_shape_manual(values = subshapes, name = NULL) +
  scale_x_continuous(breaks = c(1, 10, 20, 32)) +
  labs(x = 'Site, sorted by estimate', y = 'Course pass-rate difference', title = '(a) Site contrasts and PM') +
  theme_paper() +
  theme(legend.position = 'top', legend.margin = margin(0, 0, 0, 0), legend.key.width = unit(10, 'pt')) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
sr <- data.frame(site = 1:length(po$y), S = sort(po$pdg$S_mean))
p8b <- ggplot(sr, aes(site, S)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = po$pdg$I['2.5%'], ymax = po$pdg$I['97.5%'], fill = unname(COL['blue']), alpha = .13) +
  geom_col(fill = unname(COL['bar']), width = .85) +
  geom_hline(yintercept = po$pdg$I['mean'], color = unname(COL['blue']), linewidth = .7) +
  scale_x_continuous(breaks = c(1, 10, 20, 32)) +
  scale_y_continuous(limits = c(0, 1.08), breaks = c(0, .5, 1), expand = c(0, 0)) +
  annotate('text', x = 1, y = 1.01, label = '$I = .671\\;[.497, .812]$', hjust = 0, size = SZ, family = FAMILY) +
  labs(x = 'Site, sorted by reliability', y = 'Posterior mean reliability', title = '(b) Reliability profile') +
  theme_paper()
outcols <- c(Passing = unname(COL['blue']), Credits = unname(COL['vermillion']))
# Use readable outcome names in the plot, keeping the original data labels.
outlabels <- c(Passing = 'Course passing', Credits = 'Credits')
di <- rbind(data.frame(outcome = 'Passing', I = po$pdg$plot_draws$I), data.frame(outcome = 'Credits', I = cr$pdg$plot_draws$I))
assert(nrow(di) == 40000, 'Unexpected deterministic display cache size')
plug <- data.frame(outcome = names(outcols), post = c(po$pdg$I['mean'], cr$pdg$I['mean']), plugin = c(po$I_plug, cr$I_plug))
# Reflected Gaussian KDE avoids leaking probability outside the [0,1] support.
# The density is smoothed for display; reported summaries use quadrature.
den <- do.call(rbind, lapply(split(di, di$outcome), function(d) {
  bw <- bw.nrd0(d$I)
  a <- density(c(d$I, -d$I, 2 - d$I), bw = bw, n = 1024, from = 0, to = 1)
  data.frame(outcome = d$outcome[1], I = a$x, density = 3 * a$y, bw = bw)
}))
write.csv(den, file.path(OUT, 'fig8-density-display.csv'), row.names = FALSE)
p8c <- ggplot(den, aes(I, density, color = outcome, fill = outcome)) +
  geom_area(alpha = .16, position = 'identity') +
  geom_line(linewidth = .65) +
  geom_vline(data = plug, aes(xintercept = post, color = outcome), linewidth = .6) +
  geom_vline(data = plug, aes(xintercept = plugin, color = outcome), linetype = '22', linewidth = .6) +
  scale_color_manual(values = outcols, labels = outlabels, name = NULL) +
  scale_fill_manual(values = outcols, labels = outlabels, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1), expand = c(0, 0)) +
  labs(x = 'Informativeness $I$', y = 'Posterior density', title = '(c) Outcome comparison') +
  theme_paper() +
  theme(legend.position = 'top', legend.margin = margin(0, 0, 0, 0), legend.key.width = unit(10, 'pt'))
pos <- data.frame(
  outcome = names(outcols), I = plug$post,
  lo = c(po$pdg$I['2.5%'], cr$pdg$I['2.5%']), hi = c(po$pdg$I['97.5%'], cr$pdg$I['97.5%']),
  V = c(po$pdg$V['mean'], cr$pdg$V['mean'])
)
write.csv(pos, file.path(OUT, 'fig8-profile-position-values.csv'), row.names = FALSE)
p8d <- ggplot(pos, aes(I, V, color = outcome)) +
  geom_vline(xintercept = .5, linetype = '22', color = 'grey55', linewidth = .5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = 'y', width = .018, linewidth = .65) +
  geom_point(size = 2.7) +
  geom_text(aes(label = outlabels[outcome]), nudge_y = .048, size = SZ, family = FAMILY, show.legend = FALSE) +
  scale_color_manual(values = outcols, guide = 'none') +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(limits = c(0, .5), breaks = c(0, .25, .5)) +
  labs(x = 'Informativeness $I$', y = 'Fixed dispersion $V$', title = '(d) Profile position') +
  theme_paper()
p8 <- (p8a | p8b) / (p8c | p8d) + plot_layout(heights = c(1, 1))
savefig(p8, 'fig8_ocrs', 6.3)
assert(identical(vapply(inputs, sha, ''), unname(source_hashes)) ||
  identical(vapply(inputs, sha, ''), source_hashes), 'Inputs changed during generation')
jsonlite::write_json(list(
  script = list(path = 'analysis/scripts/07_figures.R', sha256 = sha('analysis/scripts/07_figures.R')),
  sources = as.list(source_hashes), figures = specs,
  assertions = 'Passed: input checksums; Figure 2 targets; Figure 4 identities; Figure 5 all 72 cells and 24 V=2 comparisons; Figure 6 27 rows, nine comparisons, and selected predictors; Figure 7 means and MCSEs from 500 replications; Figure 8 summary values and site ordering.',
  generated = format(Sys.time(), tz = 'UTC', usetz = TRUE)
), file.path(OUT, 'generation-receipt.json'), pretty = TRUE, auto_unbox = TRUE)
message('All eight PDFs available; run analysis/scripts/07_figure_qa.py for PDF QA.')
