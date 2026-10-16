# Grid-based analysis of the Online Credit Recovery Study.
# Supplies HC2 helper functions and fits for numerical comparisons.
# The paper uses adaptive quadrature: see 05_posterior_stability.R
# and 05_ocrs_audit.R. Use --historical-replay to run the grid comparison.
#
# Site effects are treatment coefficients adjusted for randomization blocks,
# with HC2 standard errors. A site is a school-by-subject combination.
# Single-arm blocks are absorbed by the block fixed effects.
# Sensitivities examine one site per school, log odds, and first-stage uncertainty
# from 500 student bootstraps within site, treatment arm, and randomization block.
# Results are written under analysis/outputs/historical-ocrs/.


if (sys.nframe() == 0L && !"--historical-replay" %in% commandArgs(trailingOnly = TRUE))
  stop("Use 05_ocrs_audit.R for the reported analysis. To run the grid comparison, add --historical-replay.")
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
analysis_dir <- dirname(dirname(dirname(script_path)))
setwd(analysis_dir)
source(file.path("scripts", "legacy", "00_utils.R"))
set.seed(20260731)
DATA <- Sys.getenv(
  "OCRS_DATA",
  file.path("..", "data", "OCRS_publicuse_impactanalysisdata.csv")
)
stopifnot(file.exists(DATA))
if (digest::digest(DATA, file = TRUE, algo = "sha256") !=
  "d20cbbaaf07ee8e7aad7e5d0e516b74f0a1114617fb3c7db9a4c7c2afcb681c7")
  stop("OCRS_DATA SHA256 differs from the manuscript input")
d <- read.csv(DATA, stringsAsFactors = FALSE)
d$subject <- ifelse(d$crssub == "A1", "Algebra 1", "English 9")
d$site <- paste(d$studyschids, d$subject, sep = " x ")
d$credits_y4 <- ifelse(d$subject == "Algebra 1", d$cby_y4_math, d$cby_y4_eng)
stopifnot(!anyNA(d$scm_passcrs), !anyNA(d$ragrp), !anyNA(d$credits_y4))

# Block-adjusted site effects
hc2_ols <- function(X, y, coef_name = "ragrp") {
  qrX <- qr(X)
  keep <- qrX$pivot[seq_len(qrX$rank)]
  Xk <- X[, keep, drop = FALSE]
  XtXi <- chol2inv(chol(crossprod(Xk)))
  bhat <- XtXi %*% crossprod(Xk, y)
  e <- as.numeric(y - Xk %*% bhat)
  h <- rowSums((Xk %*% XtXi) * Xk)
  meat <- crossprod(Xk * (e^2 / pmax(1 - h, 1e-8)), Xk)
  Vb <- XtXi %*% meat %*% XtXi
  i <- match(coef_name, colnames(Xk))
  c(est = bhat[i], se = sqrt(Vb[i, i]))
}

site_effect_block <- function(dd, outcome, covars = NULL) {
  # Include an intercept, treatment, block indicators, and optional covariates.
  X <- cbind(`(Intercept)` = 1, ragrp = dd$ragrp)
  if (length(unique(dd$rablock)) > 1) {
    B <- model.matrix(~ factor(rablock), data = dd)[, -1, drop = FALSE]
    X <- cbind(X, B)
  }
  if (!is.null(covars)) X <- cbind(X, as.matrix(dd[, covars, drop = FALSE]))
  hc2_ols(X, dd[[outcome]])
}

sites <- do.call(rbind, lapply(split(d, d$site), function(dd) {
  un <- site_effect_block(dd, "scm_passcrs")
  ad <- site_effect_block(
    transform(dd, pre = if (dd$subject[1] == "Algebra 1") state_math_imp else state_ela_imp),
    "scm_passcrs",
    covars = c("pre", "state_miss", "gpa9th_imp", "gpa9th_miss")
  )
  cr <- site_effect_block(dd, "credits_y4")

  # Compare with the unadjusted course-passing difference.
  y1 <- dd$scm_passcrs[dd$ragrp == 1]
  y0 <- dd$scm_passcrs[dd$ragrp == 0]
  data.frame(
    site = dd$site[1], school = dd$studyschids[1], subject = dd$subject[1],
    n1 = length(y1), n0 = length(y0),
    p1 = mean(y1), p0 = mean(y0),
    tau_hat = un["est"], se = un["se"],
    tau_hat_adj = ad["est"], se_adj = ad["se"],
    tau_hat_cr = cr["est"], se_cr = cr["se"],
    tau_hat_pool = mean(y1) - mean(y0),
    se_pool = sqrt(var(y1) / length(y1) + var(y0) / length(y0)),
    n_blocks = length(unique(dd$rablock)), row.names = NULL
  )
}))
sites <- sites[order(sites$subject, sites$school), ]
cat(sprintf(
  "J = %d sites (%d Algebra, %d English); n_j: min %d, median %.0f, max %d\n",
  nrow(sites), sum(sites$subject == "Algebra 1"),
  sum(sites$subject == "English 9"),
  min(sites$n1 + sites$n0), median(sites$n1 + sites$n0),
  max(sites$n1 + sites$n0)
))
cat(sprintf(
  "block adjustment: max |tau_blk - tau_pooled| = %.4f, median = %.4f\n",
  max(abs(sites$tau_hat - sites$tau_hat_pool)),
  median(abs(sites$tau_hat - sites$tau_hat_pool))
))

# Grid fits and posterior summaries
analyze <- function(y, v, label, M = 20000L, seed = 42) {
  fit <- fit_rubin_grid(y, v, M = M, seed = seed)
  pdg <- posterior_diagnostics(fit, v)
  gm <- exp(mean(log(v)))
  I_plug <- mean(fit$sigma^2) / (mean(fit$sigma^2) + gm)
  list(
    label = label, fit = fit, pdg = pdg, v = v, y = y,
    gm = gm, I_plug = I_plug,
    p_I_gt_half = mean(pdg$draws$I > 0.5),
    sigma_summ = c(mean = mean(fit$sigma), quantile(fit$sigma, c(0.025, 0.5, 0.975)))
  )
}

fits <- list(
  pooled = analyze(sites$tau_hat, sites$se^2, "Pooled (J = 32)"),
  algebra = with(
    subset(sites, subject == "Algebra 1"),
    analyze(tau_hat, se^2, "Algebra 1 (J = 13)")
  ),
  english = with(
    subset(sites, subject == "English 9"),
    analyze(tau_hat, se^2, "English 9 (J = 19)")
  ),
  pooled_adj = analyze(
    sites$tau_hat_adj, sites$se_adj^2,
    "Pooled, covariate-adjusted (J = 32)"
  ),
  credits = analyze(
    sites$tau_hat_cr, sites$se_cr^2,
    "Credits by year 4, subject-matched (J = 32)"
  )
)
summarize_fit <- function(a) {
  s <- function(x) sprintf("%.3f [%.3f, %.3f]", x["mean"], x["2.5%"], x["97.5%"])
  data.frame(
    analysis = a$label, J = length(a$y),
    tau = sprintf("%.3f", mean(a$fit$tau)),
    sigma = sprintf(
      "%.3f [%.3f, %.3f]", a$sigma_summ["mean"],
      a$sigma_summ["2.5%"], a$sigma_summ["97.5%"]
    ),
    I = s(a$pdg$I), I_plug = sprintf("%.3f", a$I_plug),
    P_I_gt_half = sprintf("%.2f", a$p_I_gt_half),
    H = s(a$pdg$H),
    V_logit = sprintf("%.3f", a$pdg$V["mean"]),
    J_eff = sprintf(
      "%.1f [%.1f, %.1f]", a$pdg$J_eff["mean"],
      a$pdg$J_eff["2.5%"], a$pdg$J_eff["97.5%"]
    )
  )
}

diag_tab <- do.call(rbind, lapply(fits, summarize_fit))
print(diag_tab, row.names = FALSE)
sites$S_post <- fits$pooled$pdg$S_mean
sites$pm <- fits$pooled$fit$pm

# Sensitivity analyses
sens <- list()

# Eight schools contribute two subject sites; keep one site from each.
paired <- names(which(table(sites$school) == 2))
cat("\npaired schools:", length(paired), "\n")
grid8 <- expand.grid(rep(list(c("Algebra 1", "English 9")), length(paired)),
  stringsAsFactors = FALSE
)
uni <- t(vapply(seq_len(nrow(grid8)), function(k) {
  keep <- rep(TRUE, nrow(sites))
  for (pi in seq_along(paired)) {
    drop_subj <- setdiff(c("Algebra 1", "English 9"), grid8[k, pi])
    keep[sites$school == paired[pi] & sites$subject == drop_subj] <- FALSE
  }
  s2 <- sites[keep, ]
  fitp <- fit_rubin_grid(s2$tau_hat, s2$se^2, M = 4000L, seed = 1000 + k)
  fitc <- fit_rubin_grid(s2$tau_hat_cr, s2$se_cr^2, M = 4000L, seed = 5000 + k)
  gm_p <- exp(mean(log(s2$se^2)))
  gm_c <- exp(mean(log(s2$se_cr^2)))
  c(
    I_pass = mean(fitp$sigma^2 / (fitp$sigma^2 + gm_p)),
    I_cred = mean(fitc$sigma^2 / (fitc$sigma^2 + gm_c))
  )
}, numeric(2)))
sens$unique_school <- data.frame(
  outcome = c("passing", "credits"),
  I_min = c(min(uni[, 1]), min(uni[, 2])),
  I_med = c(median(uni[, 1]), median(uni[, 2])),
  I_max = c(max(uni[, 1]), max(uni[, 2]))
)
cat("unique-school sensitivity (2^8 = 256 selections):\n")
print(round(sens$unique_school[, -1], 3))

# Log odds with Woolf variances; add 0.5 only when a cell count is zero.
lo <- do.call(rbind, lapply(split(d, d$site), function(dd) {
  a <- sum(dd$scm_passcrs[dd$ragrp == 1])
  b <- sum(dd$ragrp == 1) - a
  c0 <- sum(dd$scm_passcrs[dd$ragrp == 0])
  e <- sum(dd$ragrp == 0) - c0
  cc <- if (min(a, b, c0, e) == 0) 0.5 else 0
  data.frame(
    site = dd$site[1],
    lor = log((a + cc) / (b + cc)) - log((c0 + cc) / (e + cc)),
    vlor = 1 / (a + cc) + 1 / (b + cc) + 1 / (c0 + cc) + 1 / (e + cc),
    zero_cell = cc > 0, row.names = NULL
  )
}))
cat("\nlog-odds scale: sites needing continuity correction:", sum(lo$zero_cell), "\n")
flo <- analyze(lo$lor, lo$vlor, "Passing, log-odds scale (J = 32)", seed = 7)
sens$logodds <- data.frame(
  scale = c("risk difference", "log odds"),
  I = c(
    sprintf(
      "%.3f [%.3f, %.3f]", fits$pooled$pdg$I["mean"],
      fits$pooled$pdg$I["2.5%"], fits$pooled$pdg$I["97.5%"]
    ),
    sprintf("%.3f [%.3f, %.3f]", flo$pdg$I["mean"], flo$pdg$I["2.5%"], flo$pdg$I["97.5%"])
  ),
  P_I_gt_half = c(fits$pooled$p_I_gt_half, flo$p_I_gt_half),
  H = c(sprintf("%.3f", fits$pooled$pdg$H["mean"]), sprintf("%.3f", flo$pdg$H["mean"])),
  V = c(sprintf("%.3f", fits$pooled$pdg$V["mean"]), sprintf("%.3f", flo$pdg$V["mean"]))
)
print(sens$logodds, row.names = FALSE)

# Resample students within site, arm, and randomization block. These intervals
# measure first-stage resampling uncertainty, not conditional posterior uncertainty.
# A resampled arm may be constant, giving zero HC2 variance and undefined log(v).
# For course passing, use a minimum variance based on p = 0.5 / (n + 1).
# For credits, use 5% of the original site sampling variance.
boot_one <- function(outcome, B = 500L, M = 6000L) {
  split_cells <- split(seq_len(nrow(d)), list(d$site, d$ragrp, d$rablock), drop = TRUE)
  if (outcome == "scm_passcrs") {
    v_floor <- vapply(split(d, d$site), function(dd) {
      n1 <- sum(dd$ragrp == 1)
      n0 <- sum(dd$ragrp == 0)
      pm1 <- 0.5 / (n1 + 1)
      pm0 <- 0.5 / (n0 + 1)
      pm1 * (1 - pm1) / n1 + pm0 * (1 - pm0) / n0
    }, numeric(1))
  } else {
    v_floor <- 0.05 * vapply(
      split(d, d$site),
      function(dd) site_effect_block(dd, outcome)["se"]^2,
      numeric(1)
    )
  }
  v_floor <- v_floor[sort(names(v_floor))]
  res <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("I", "H", "V")))
  for (b in seq_len(B)) {
    idx <- unlist(lapply(split_cells, function(ii) ii[sample.int(length(ii),
      length(ii),
      replace = TRUE
    )]))
    db <- d[idx, ]
    db <- db[order(db$site), ]
    sb <- do.call(rbind, lapply(
      split(db, db$site),
      function(dd) site_effect_block(dd, outcome)
    ))
    tau_b <- sb[, "est"]
    v_b <- pmax(sb[, "se"]^2, v_floor)
    fit <- fit_rubin_grid(tau_b, v_b, M = M, seed = 20000 + b)
    gm <- exp(mean(log(v_b)))
    Ib <- mean(fit$sigma^2 / (fit$sigma^2 + gm))
    Smat <- outer(fit$sigma^2, v_b, function(s2, vv) s2 / (s2 + vv))
    Lmat <- log(Smat) - log1p(-Smat)
    Ipd <- plogis(rowMeans(Lmat))
    Hb <- mean(vapply(
      seq_len(nrow(Smat)),
      function(m) mean(kl_bern(Ipd[m], Smat[m, ])), numeric(1)
    ))
    Lv <- log(v_b)
    res[b, ] <- c(Ib, Hb, mean((Lv - mean(Lv))^2))
  }
  res
}

cat("\nfirst-stage bootstrap (B = 500) ...\n")
bs_pass <- boot_one("scm_passcrs")
bs_cred <- boot_one("credits_y4")
qs <- function(x) sprintf("[%.3f, %.3f]", quantile(x, .025), quantile(x, .975))
sens$firststage <- data.frame(
  outcome = c("passing", "credits"),
  I_boot95 = c(qs(bs_pass[, "I"]), qs(bs_cred[, "I"])),
  H_boot95 = c(qs(bs_pass[, "H"]), qs(bs_cred[, "H"])),
  V_boot95 = c(qs(bs_pass[, "V"]), qs(bs_cred[, "V"]))
)
print(sens$firststage, row.names = FALSE)

# Save grid-based results separately from the quadrature-based analysis.
dir.create("outputs/historical-ocrs/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(sites, "outputs/historical-ocrs/tables/tbl_ocrs_sites.csv", row.names = FALSE)
write.csv(diag_tab, "outputs/historical-ocrs/tables/tbl_ocrs_diagnostics.csv", row.names = FALSE)
sens_out <- rbind(
  data.frame(
    block = "unique_school", metric = sens$unique_school$outcome,
    value = sprintf(
      "I in [%.3f, %.3f], median %.3f",
      sens$unique_school$I_min, sens$unique_school$I_max,
      sens$unique_school$I_med
    )
  ),
  data.frame(
    block = "logodds", metric = sens$logodds$scale,
    value = sprintf(
      "I = %s, Pr(I > .5) = %.2f, H = %s, V = %s",
      sens$logodds$I, as.numeric(sens$logodds$P_I_gt_half),
      sens$logodds$H, sens$logodds$V
    )
  ),
  data.frame(
    block = "firststage_boot", metric = sens$firststage$outcome,
    value = sprintf(
      "I %s, H %s, V %s", sens$firststage$I_boot95,
      sens$firststage$H_boot95, sens$firststage$V_boot95
    )
  )
)
write.csv(sens_out, "outputs/historical-ocrs/tables/tbl_ocrs_sensitivity.csv", row.names = FALSE)
saveRDS(
  list(fits = fits, sens = sens, logodds_fit = flo),
  "outputs/historical-ocrs/ocrs_fits.rds"
)
cat("\nwritten: tbl_ocrs_sites.csv, tbl_ocrs_diagnostics.csv, tbl_ocrs_sensitivity.csv, ocrs_fits.rds\n")
writeLines(capture.output(sessionInfo()), "outputs/historical-ocrs/sessionInfo_ocrs.txt")
