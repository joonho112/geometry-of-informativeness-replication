# Simulation at fixed total enrollment N = J * nbar, with v_j = kappa / n_j.
# Unlike the fixed-I arm, this arm compares site sizes under a common budget.
# The balanced-design theorem concerns informativeness under its variance
# model; observed loss differences remain finite-simulation comparisons.
#
# Factors: J = 25, 50, 100; nbar = 20, 50, 100; parent Gamma CV = 0, .25, .5, .75;
# Gaussian, mixture, and ALD effects; kappa = 4; sigma_tau = .20; REPS = 500.
# Each cell uses one site-size vector for all replications. Latent effects and
# sampling errors are shared across CV levels within each (J, nbar, G) family.
# Sizes are clipped at 5 and rescaled twice to sum to N. They are continuous;
# final sizes can be below 5, and their CV need not equal the parent Gamma CV.
# Run through analysis/scripts/run_historical.R --arm=fixedN.

source("00_utils.R")
suppressPackageStartupMessages({
  library(parallel)
  source("historical-dgp.R")
})

SIGMA_T <- 0.20
KAPPA <- 4
REPS    <- 500L
M_DRAWS <- 2000L
SEED0 <- 20260731L

families <- expand.grid(
  J = c(25L, 50L, 100L),
  nbar = c(20, 50, 100),
  G = c("Gaussian", "Mixture", "ALD"),
  stringsAsFactors = FALSE
)
CVs <- c(0, 0.25, 0.5, 0.75)
cat("families:", nrow(families), " x CV levels:", length(CVs),
    " = cells:", nrow(families) * length(CVs), " reps/cell:", REPS, "
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

# Clip and rescale Gamma sizes twice; the final rescaling can give sizes below 5.
draw_sizes <- function(J, nbar, cv) {
  N <- J * nbar
  if (cv <= 0) return(rep(nbar, J))
  n <- pmax(5, rgamma(J, shape = 1 / cv^2, scale = nbar * cv^2))
  n <- n * (N / sum(n))
  n <- pmax(5, n)
  n * (N / sum(n))
}

run_family <- function(f) {
  g <- families[f, ]
  set.seed(SEED0 + f)
  TAU <- gen_tau_matrix(g$G, g$J, REPS)
  EPS <- matrix(rnorm(REPS * g$J), REPS, g$J)
  size_list <- lapply(CVs, function(cv) draw_sizes(g$J, g$nbar, cv))

  fam_rows <- vector("list", length(CVs))
  for (ci in seq_along(CVs)) {
    n_j <- size_list[[ci]]
    v <- KAPPA / n_j
    dg <- diagnostics_profile(SIGMA_T^2, v)
    S_true <- dg$S
    sv <- sqrt(v)

    rows <- vector("list", REPS)
    for (r in seq_len(REPS)) {
      tau_j <- TAU[r, ]
      y <- tau_j + EPS[r, ] * sv
      set.seed(SEED0 + 1000L * f + 61L * ci + 200000L + r)
      fit <- fit_rubin_grid(y, v, M = M_DRAWS)
      ests <- list(
        ml = y, opm = S_true * y,
        pm = est_pm(fit), cb = est_cb(fit), gr = est_gr(fit)
      )
      L <- do.call(rbind, lapply(ests, all_losses, tau = tau_j))
      rows[[r]] <- data.frame(
        family = f, cv_i = ci, rep = r, est = rownames(L),
        rmse = L[, "rmse"], mselp = L[, "mselp"], isel = L[, "isel"],
        row.names = NULL
      )
    }
    out <- do.call(rbind, rows)
    out$J <- g$J
    out$nbar <- g$nbar
    out$G <- g$G
    out$CV <- CVs[ci]
    out$I <- dg$I
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
res$cell <- as.integer(interaction(res$family, res$cv_i, drop = TRUE))
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

dir.create("outputs", showWarnings = FALSE)
saveRDS(res, "outputs/sim_fixedN_reps.rds", compress = TRUE)
cells <- aggregate(
  cbind(rmse, mselp, isel) ~
    cell + est + J + nbar + G + CV + I + H + V_logit + J_eff + S_bar,
  data = res, FUN = mean
)
write.csv(cells, "outputs/sim_fixedN_cells.csv", row.names = FALSE)
cat("written: outputs/sim_fixedN_reps.rds, outputs/sim_fixedN_cells.csv\n")
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_sim_fixedN.txt")
