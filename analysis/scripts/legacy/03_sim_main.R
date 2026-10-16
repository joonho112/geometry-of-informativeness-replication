# Simulation with specified informativeness I and logit variance V.
# Within each (J, I, G) family, profiles share latent effects and sampling errors
# so that loss comparisons are paired. M = 2,000 posterior draws, REPS = 500.
#
# Rao-Blackwellized moments preserve the PM/raw rank equality when v_j is
# constant. The two concentrated profiles cover opposite skew directions.
# The semi-oracle estimator fixes sigma2 and estimates the overall mean.
#
# Factors: J = 25, 50, 100; I = .05, .15, .30, .45, .60, .75;
# V = 0, .5, 1, 2; four shapes when V > 0; Gaussian, mixture, and ALD effects.
# sigma_tau = .20. Estimators: ml, opm, spm, pm, cb, gr.
# opm is the known-hyperparameter Gaussian linear benchmark.
# Run through analysis/scripts/run_historical.R --arm=stress.

source("00_utils.R")
suppressPackageStartupMessages({
  library(parallel)
  source("historical-dgp.R")
})

SIGMA_T <- 0.20
REPS    <- 500L
M_DRAWS <- 2000L
SEED0 <- 20260710L

families <- expand.grid(
  J = c(25L, 50L, 100L),
  I = c(0.05, 0.15, 0.30, 0.45, 0.60, 0.75),
  G = c("Gaussian", "Mixture", "ALD"),
  stringsAsFactors = FALSE
)
combos <- rbind(
  data.frame(V = 0, shape = "even", stringsAsFactors = FALSE),
  expand.grid(
    V = c(0.5, 1, 2),
    shape = c("even", "bimodal", "concentrated", "concentrated_noisy"),
    stringsAsFactors = FALSE
  )
)
cat("families:", nrow(families), " combos/family:", nrow(combos),
    " cells:", nrow(families) * nrow(combos), " reps/cell:", REPS, "
")

gen_tau_matrix <- function(G, J, reps) {
  n <- J * reps
  z <- switch(G,
    Gaussian = gen_effects_gaussian(J = n, tau = 0, sigma_tau = SIGMA_T)$tau_j,
    Mixture = gen_effects_mixture(
      J = n, tau = 0, sigma_tau = SIGMA_T,
      delta = 5, eps = 0.3, ups = 2
    )$tau_j,
    ALD = gen_effects_ald(J = n, tau = 0, sigma_tau = SIGMA_T, rho = 0.3)$tau_j
  )
  matrix(z, nrow = reps, ncol = J)
}

run_family <- function(f) {
  g <- families[f, ]
  gm <- SIGMA_T^2 * (1 - g$I) / g$I

  # common random numbers for the whole family
  set.seed(SEED0 + f)
  TAU <- gen_tau_matrix(g$G, g$J, REPS)
  EPS <- matrix(rnorm(REPS * g$J), REPS, g$J)

  fam_rows <- vector("list", nrow(combos))
  for (ci in seq_len(nrow(combos))) {
    cmb <- combos[ci, ]
    v <- make_profile(cmb$shape, g$J, gm, cmb$V)
    dg <- diagnostics_profile(SIGMA_T^2, v)
    S_true <- dg$S
    sv <- sqrt(v)

    rows <- vector("list", REPS)
    for (r in seq_len(REPS)) {
      tau_j <- TAU[r, ]
      y <- tau_j + EPS[r, ] * sv
      set.seed(SEED0 + 1000L * f + 61L * ci + 100000L + r) # Posterior-draw seed.
      fit <- fit_rubin_grid(y, v, M = M_DRAWS)

      ests <- list(
        ml  = y,
        opm = S_true * y, # known-(tau, sigma2) Gaussian linear benchmark
        spm = est_spm(y, v, SIGMA_T^2), # sigma2 known, tau estimated
        pm  = est_pm(fit), # Rao-Blackwellized
        cb  = est_cb(fit),
        gr  = est_gr(fit)
      )
      L <- do.call(rbind, lapply(ests, all_losses, tau = tau_j))

      gm_v <- exp(mean(log(v)))
      I_post <- mean(fit$sigma^2 / (fit$sigma^2 + gm_v))

      rows[[r]] <- data.frame(
        family = f, combo = ci, rep = r, est = rownames(L),
        rmse = L[, "rmse"], mselp = L[, "mselp"], isel = L[, "isel"],
        I_hat = I_post, row.names = NULL
      )
    }
    out <- do.call(rbind, rows)
    out$J <- g$J
    out$I <- g$I
    out$G <- g$G
    out$V <- cmb$V
    out$shape <- cmb$shape
    out$H <- dg$H
    out$V_logit <- dg$V
    out$J_eff <- dg$J_eff
    out$S_bar <- dg$S_bar
    fam_rows[[ci]] <- out
  }
  do.call(rbind, fam_rows)
}

t0 <- Sys.time()
res_list <- mclapply(seq_len(nrow(families)), run_family,
  mc.cores = max(1L, detectCores() - 2L),
  mc.preschedule = FALSE
)
bad <- vapply(res_list, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) stop("failed families: ", paste(which(bad), collapse = ", "))
res <- do.call(rbind, res_list)
res$cell <- interaction(res$family, res$combo, drop = TRUE)
res$cell <- as.integer(res$cell)
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

dir.create("outputs", showWarnings = FALSE)
saveRDS(res, "outputs/sim_main_reps.rds", compress = TRUE)

cells <- aggregate(
  cbind(rmse, mselp, isel, I_hat) ~
    cell + est + J + I + V + shape + G + H + V_logit + J_eff + S_bar,
  data = res, FUN = mean
)
write.csv(cells, "outputs/sim_main_cells.csv", row.names = FALSE)
cat("written: outputs/sim_main_reps.rds, outputs/sim_main_cells.csv\n")
cat("rows (reps):", nrow(res), " rows (cells):", nrow(cells), "\n")

writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_sim_main.txt")
