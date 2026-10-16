# Recompute simulation tables from the saved replication losses.
# Run from analysis/: Rscript scripts/05_simulation_postprocess.R
# Then run 05_evidence_audit.R, 05_fixedN_design_audit.R, and
# 05_paired_component_mcse.R, which replaces the preliminary paired SEs.
#
# Metamodel predictors describe the true generating profiles, not estimates.
# Panel A measures residual variation explained on the simulation grid.
# Panel B compares summaries using the same model form and five-fold CV.
# Plots here are working plots; 07_figures.R produces manuscript Figures 1–8.

source("scripts/legacy/00_utils.R")
d <- read.csv("inputs/sim_main_cells.csv")
r <- readRDS("inputs/sim_main_reps.rds")
dir.create("outputs/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/figures", showWarnings = FALSE, recursive = TRUE)

GOALS <- c("rmse", "mselp", "isel")
FEAS <- c("ml", "pm", "cb", "gr")

# Equal sampling variances give the same MSELP for ML, PM, and CB.

h <- subset(r, V == 0 & est %in% c("ml", "pm", "cb"))
hw <- reshape(h[, c("cell", "rep", "est", "mselp")],
  idvar = c("cell", "rep"),
  timevar = "est", direction = "wide"
)
tie_pm <- max(abs(hw$mselp.pm - hw$mselp.ml))
tie_cb <- max(abs(hw$mselp.cb - hw$mselp.ml))
cat(sprintf(
  "V=0 rank-tie invariant: max|pm-ml| = %.2e, max|cb-ml| = %.2e\n",
  tie_pm, tie_cb
))
if (tie_pm > 1e-12 || tie_cb > 1e-12) stop("homogeneous rank-tie invariant violated")

# Metamodels

# I, H, V_logit and S_bar describe the true generating profile (sigma2=0.04
# and known sampling variances), not estimates from the simulated observations.
# The separately saved I_hat column is not a predictor in these models.
# Cross-validation splits cells within the stated grid; it is not a test of
# extrapolation or of prediction using estimated profile diagnostics.

# Panel A: residual variation explained by dispersion and profile shape.
panelA <- do.call(rbind, lapply(GOALS, function(goal) {
  do.call(rbind, lapply(FEAS, function(es) {
    dd <- subset(d, est == es)
    y <- log(dd[[goal]])
    m0 <- lm(y ~ factor(I) * factor(J) * factor(G), data = dd)
    m1 <- lm(y ~ factor(I) * factor(J) * factor(G) +
      factor(V) * factor(shape), data = dd)
    r0 <- summary(m0)$r.squared
    r1 <- summary(m1)$r.squared
    data.frame(goal,
      est = es, R2_base = r0, R2_profile = r1,
      resid_share = (r1 - r0) / (1 - r0)
    )
  }))
}))
write.csv(panelA, "outputs/tables/tbl_metamodel.csv", row.names = FALSE)
cat("\n== Panel A: profile-beyond-I residual share ==\n")
print(transform(panelA,
  R2_base = round(R2_base, 4), R2_profile = round(R2_profile, 4),
  resid_share = round(resid_share, 3)
), row.names = FALSE)

# Panel B: compare single summaries using the same model form.
cv_r2 <- function(dd, y, form, K = 5, seed = 1) {
  set.seed(seed)
  fold <- sample(rep_len(seq_len(K), nrow(dd)))
  press <- 0
  for (k in seq_len(K)) {
    m <- lm(form, data = transform(dd[fold != k, ], .y = y[fold != k]))
    press <- press + sum((y[fold == k] -
      predict(m, newdata = dd[fold == k, ]))^2)
  }
  1 - press / sum((y - mean(y))^2)
}
panelB <- do.call(rbind, lapply(GOALS, function(goal) {
  do.call(rbind, lapply(c("pm", "cb", "gr"), function(es) {
    dd <- subset(d, est == es)
    y <- log(dd[[goal]])
    dd$.y <- y
    do.call(rbind, lapply(c("H", "V_logit", "S_bar"), function(x) {
      form <- as.formula(paste(
        ".y ~ factor(I)*factor(J)*factor(G) +", x,
        "+", x, ":factor(I)"
      ))
      m <- lm(form, data = dd)
      data.frame(goal,
        est = es, summary = x,
        adjR2 = summary(m)$adj.r.squared, AIC = AIC(m),
        cvR2 = cv_r2(dd, y, form)
      )
    }))
  }))
}))
write.csv(panelB, "outputs/tables/tbl_fair.csv", row.names = FALSE)
cat("\n== Panel B: fair contest (same structure, same df) ==\n")
print(transform(panelB,
  adjR2 = round(adjR2, 5), AIC = round(AIC, 1),
  cvR2 = round(cvR2, 5)
), row.names = FALSE)

# Paired loss ratios: V = 2 relative to V = 0.

base0 <- subset(r, V == 0)[, c("family", "rep", "est", "rmse", "mselp", "isel")]
names(base0)[4:6] <- paste0(c("rmse", "mselp", "isel"), "0")

paired_ratios <- function(shp) {
  m2 <- merge(subset(r, V == 2 & shape == shp),
    base0,
    by = c("family", "rep", "est")
  )
  out <- do.call(rbind, lapply(
    split(m2, list(m2$est, m2$I), drop = TRUE),
    function(x) {
      lr <- cbind(
        log(x$rmse / x$rmse0), log(x$mselp / x$mselp0),
        log(x$isel / x$isel0)
      )
      # average paired log-ratios within family first, then across families
      fam <- rowsum(lr, x$family) / as.vector(table(x$family))
      mns <- colMeans(fam)
      ses <- sqrt(Reduce(`+`, lapply(split(as.data.frame(lr), x$family), function(z) apply(z, 2, var) / nrow(z)))) / nrow(fam)
      # These preliminary SEs omit shared fit seeds; the paired MCSE script replaces them.
      heterogeneity <- apply(fam, 2, sd)
      data.frame(
        est = x$est[1], I = x$I[1], shape = shp,
        r_rmse = exp(mns[1]), se_rmse = exp(mns[1]) * ses[1],
        r_mselp = exp(mns[2]), se_mselp = exp(mns[2]) * ses[2],
        r_isel = exp(mns[3]), se_isel = exp(mns[3]) * ses[3],
        logratio_family_sd_rmse = heterogeneity[1],
        logratio_family_sd_mselp = heterogeneity[2],
        logratio_family_sd_isel = heterogeneity[3],
        uncertainty_status = "independent-replication approximation; shared-seed adjustment follows"
      )
    }
  ))
  out[order(out$est, out$I), ]
}
ratios_even <- paired_ratios("even")
ratios_shape <- do.call(rbind, lapply(
  c("even", "bimodal", "concentrated", "concentrated_noisy"), paired_ratios
))
write.csv(ratios_even, "outputs/tables/tbl_ratios.csv", row.names = FALSE)
write.csv(ratios_shape, "outputs/tables/tbl_ratios_shape.csv", row.names = FALSE)
cat("\n== Paired loss ratios V=2/V=0 (even; provisional independent-rep SEs in parens) ==\n")
pr <- subset(ratios_even, est %in% c("pm", "gr", "opm"))
print(
  data.frame(
    est = pr$est, I = pr$I,
    rmse = sprintf("%.3f (%.3f)", pr$r_rmse, pr$se_rmse),
    mselp = sprintf("%.3f (%.3f)", pr$r_mselp, pr$se_mselp),
    isel = sprintf("%.3f (%.3f)", pr$r_isel, pr$se_isel)
  ),
  row.names = FALSE
)

# Compare known, partly estimated, and fully estimated hyperparameters.

oracle <- subset(ratios_even, est %in% c("opm", "spm", "pm"),
  select = c(est, I, r_rmse, se_rmse)
)
write.csv(oracle, "outputs/tables/tbl_oracle.csv", row.names = FALSE)
cat("\n== RMSE ratio decomposition (even, V=2 vs V=0) ==\n")
print(reshape(oracle[, c("est", "I", "r_rmse")],
  idvar = "I", timevar = "est",
  direction = "wide"
), row.names = FALSE, digits = 3)

# Relative Monte Carlo SEs of cell mean losses.

mcse <- do.call(rbind, lapply(
  split(r, list(r$cell, r$est), drop = TRUE),
  function(x) data.frame(
    cell = x$cell[1], est = x$est[1],
    rel_mcse_rmse = sd(x$rmse) / sqrt(nrow(x)) / mean(x$rmse),
    rel_mcse_mselp = sd(x$mselp) / sqrt(nrow(x)) / mean(x$mselp),
    rel_mcse_isel = sd(x$isel) / sqrt(nrow(x)) / mean(x$isel)
  )
))
mc_sum <- data.frame(
  metric = c("rmse", "mselp", "isel"),
  median_rel_mcse = c(
    median(mcse$rel_mcse_rmse), median(mcse$rel_mcse_mselp),
    median(mcse$rel_mcse_isel)
  ),
  max_rel_mcse = c(
    max(mcse$rel_mcse_rmse), max(mcse$rel_mcse_mselp),
    max(mcse$rel_mcse_isel)
  )
)
write.csv(mc_sum, "outputs/tables/tbl_mcse.csv", row.names = FALSE)
cat("\n== relative MCSE of cell means ==\n")
print(mc_sum, row.names = FALSE)
