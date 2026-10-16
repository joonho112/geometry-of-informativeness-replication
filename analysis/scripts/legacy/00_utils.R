# Reliability diagnostics, simulation profiles, Rubin-model fitting, and losses.
# Sampling variances v are treated as known; sigma2 is the cross-site variance.
# All diagnostics use natural logarithms.
#
# Posterior summaries follow Ghosh (1992), Shen and Louis (1998), and
# Lee et al. (2025). The model uses a grid on sigma, with a uniform prior
# on the overall mean and sigma; it does not use MCMC.
# ggplot2 is optional and used only by the plotting helpers.

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
}

# Reliability diagnostics

logit <- function(p) log(p) - log1p(-p)
logistic <- function(x) plogis(x)

# Site reliability, also the weight on the site estimate.
site_reliability <- function(sigma2, v) sigma2 / (sigma2 + v)

# KL-centroid informativeness (Lee et al. 2025, eq. 4): geometric-mean form
informativeness_I <- function(sigma2, v) {
  gm <- exp(mean(log(v)))
  sigma2 / (sigma2 + gm)
}

# Equivalent logit-mean form (Theorem 1).
kl_centroid <- function(S) logistic(mean(logit(S)))

# Bernoulli KL divergence KL(Ber(p) || Ber(q)), vectorized in q
kl_bern <- function(p, q) {
  p * (log(p) - log(q)) + (1 - p) * (log1p(-p) - log1p(-q))
}

# Bregman (KL) dispersion: H_KL = (1/J) sum_j KL(Ber(I) || Ber(S_j))
H_KL <- function(S, I = kl_centroid(S)) mean(kl_bern(I, S))

# Logit-scale variance (population form, 1/J)
V_logit <- function(S) {
  L <- logit(S)
  mean((L - mean(L))^2)
}

# Effective site count for the reliability profile (Kish concentration formula).
J_eff <- function(S) sum(S)^2 / sum(S^2)

# Calculate all profile summaries from sigma2 and the sampling variances.
diagnostics_profile <- function(sigma2, v) {
  S <- site_reliability(sigma2, v)
  I <- kl_centroid(S)
  list(
    S = S, I = I, H = H_KL(S, I), V = V_logit(S),
    J_eff = J_eff(S), S_bar = mean(S)
  )
}

# Sampling-variance profiles with specified geometric mean and logit variance
#
# Construction: log v_j = log(gm) + x_j with mean(x) = 0 and mean(x^2) = V.
# Then logit(S_j) = log(sigma2/gm) - x_j, giving the specified
# target informativeness I = sigma2/(sigma2 + gm) and logit variance V.

.standardize_exact <- function(x) {
  x <- x - mean(x)
  s <- sqrt(mean(x^2))
  if (s < .Machine$double.eps) return(rep(0, length(x)))
  x / s
}

# Even spread: equally spaced log-variances (log-uniform grid)
profile_even <- function(J, gm, V) {
  if (V <= 0) return(rep(gm, J))
  x <- .standardize_exact(seq(0, 1, length.out = J)) * sqrt(V)
  gm * exp(x)
}

# Symmetric bimodal: half precise, half noisy
profile_bimodal <- function(J, gm, V) {
  if (V <= 0) return(rep(gm, J))
  half <- rep(c(-1, 1), length.out = J)
  x <- .standardize_exact(sort(half)) * sqrt(V)
  gm * exp(x)
}

# Concentrated-precise: a small fraction of very precise sites
# (k sites far below gm, the rest slightly above; exact mean 0, variance V)
profile_concentrated <- function(J, gm, V, frac = 0.2) {
  if (V <= 0) return(rep(gm, J))
  k <- max(1L, round(frac * J))
  a <- -sqrt(V * (J - k) / k) # precise block (low v, high S)
  b <- sqrt(V * k / (J - k)) # diffuse block
  gm * exp(c(rep(a, k), rep(b, J - k)))
}

# Concentrated-noisy: the mirror profile, with a small fraction of very noisy
# sites (k sites far above gm, the rest slightly below). Included because the
# sign of S_bar - I at fixed (I, V) depends on the skew orientation of the
# logit profile (the O(rho_3) term of Proposition 2).
profile_concentrated_noisy <- function(J, gm, V, frac = 0.2) {
  if (V <= 0) return(rep(gm, J))
  k <- max(1L, round(frac * J))
  a <- sqrt(V * (J - k) / k) # noisy block (high v, low S)
  b <- -sqrt(V * k / (J - k)) # remaining, more precise sites
  gm * exp(c(rep(a, k), rep(b, J - k)))
}

make_profile <- function(shape = c(
                           "even", "bimodal", "concentrated",
                           "concentrated_noisy"
                         ),
                         J, gm, V) {
  shape <- match.arg(shape)
  switch(shape,
    even               = profile_even(J, gm, V),
    bimodal            = profile_bimodal(J, gm, V),
    concentrated       = profile_concentrated(J, gm, V),
    concentrated_noisy = profile_concentrated_noisy(J, gm, V)
  )
}

# Rubin model: grid approximation over sigma
#
# Model: y_j | tau_j ~ N(tau_j, v_j), tau_j | tau, sigma ~ N(tau, sigma^2),
# p(tau, sigma) ∝ 1 on (−Inf, Inf) × [0, sigma_max]. The marginal posterior
# p(sigma | y) is available in closed form (Rubin 1981; Gelman et al. 2013,
# sec. 5.4); we sample sigma on a fine grid, then tau | sigma, then tau_j.

log_post_sigma <- function(sigma, y, v) {
  # vectorized over sigma
  vapply(sigma, function(s) {
    w <- 1 / (v + s^2)
    mu <- sum(w * y) / sum(w)
    0.5 * log(1 / sum(w)) + sum(-0.5 * log(v + s^2) - 0.5 * w * (y - mu)^2)
  }, numeric(1))
}

fit_rubin_grid <- function(y, v, M = 2000L, n_grid = 250L, sigma_max = NULL,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  J <- length(y)
  stopifnot(length(v) == J, all(v > 0))

  if (is.null(sigma_max)) {
    mom <- sqrt(max(0, var(y) - mean(v)))
    sigma_max <- max(3 * mom, 2 * sd(y), 3 * sqrt(median(v)), 1e-3)
  }

  # grid on sigma (uniform prior on sigma)
  repeat {
    sg <- seq(0, sigma_max, length.out = n_grid)
    sg[1] <- sigma_max / (10 * n_grid) # avoid sigma = 0 degeneracies
    lp <- log_post_sigma(sg, y, v)
    p <- exp(lp - max(lp))
    p <- p / sum(p)
    # extend the grid if the upper tail still carries mass
    if (sum(p[sg > 0.95 * sigma_max]) < 1e-4 || sigma_max > 100 * sd(y)) break
    sigma_max <- 2 * sigma_max
  }

  idx <- sample.int(n_grid, M, replace = TRUE, prob = p)
  step <- sg[2] - sg[1]
  sigma <- pmax(sg[idx] + runif(M, -step / 2, step / 2), step / 20)

  # tau | sigma, y
  Vmat <- outer(sigma^2, v, `+`) # M x J: sigma^2 + v_j
  W <- 1 / Vmat
  sw <- rowSums(W)
  mu_h <- rowSums(W * matrix(y, M, J, byrow = TRUE)) / sw
  tau <- rnorm(M, mu_h, sqrt(1 / sw))

  # tau_j | tau, sigma, y  (conditional posterior, eq. 3 of Lee et al. 2025)
  Smat <- sigma^2 / Vmat # M x J shrinkage weights
  mean_ <- matrix(tau, M, J) + Smat * (matrix(y, M, J, byrow = TRUE) - matrix(tau, M, J))
  sd_ <- sqrt(sigma^2 * (1 - Smat))
  theta <- mean_ + matrix(rnorm(M * J), M, J) * sd_

  # Rao-Blackwellized posterior means and variances: average the conditional
  # moments over the hyperparameter draws instead of the sampled site effects.
  # pm_j = M^-1 sum_m E(tau_j | tau^m, sigma^m, y) is affine in y_j with a
  # common positive slope whenever the v_j are equal, so posterior-mean ranks
  # tie exactly with the raw ranks in homogeneous profiles (no Monte Carlo
  # rank noise); pv_j uses the law of total variance.
  pm_rb <- colMeans(mean_)
  pv_rb <- colMeans(sd_^2) + colMeans(mean_^2) - pm_rb^2

  list(
    theta = theta, # M x J joint posterior draws of tau_j
    sigma = sigma, tau = tau, # hyperparameter draws
    pm = pm_rb, # Rao-Blackwellized posterior means
    pv = pv_rb, # Rao-Blackwellized posterior variances
    grid = list(sigma = sg, prob = p)
  )
}

# Semi-oracle posterior mean: sigma^2 known at its true value, tau estimated
# (flat prior on tau). E[tau_j | y, sigma2] = S_j y_j + (1 - S_j) mu_hat.
# Used to separate the sigma^2-learning channel from the tau-learning channel.
est_spm <- function(y, v, sigma2) {
  w <- 1 / (sigma2 + v)
  mu <- sum(w * y) / sum(w)
  S <- sigma2 / (sigma2 + v)
  S * y + (1 - S) * mu
}

# Posterior summaries

est_pm <- function(fit) fit$pm

# Constrained Bayes (Ghosh 1992; Lee et al. 2025, eq. 9-10)
est_cb <- function(fit) {
  pm <- fit$pm
  J <- length(pm)
  tbar <- mean(pm)
  vfin <- sum((pm - tbar)^2) / (J - 1)
  sig2h <- mean(fit$pv) + vfin
  if (vfin < .Machine$double.eps) return(rep(tbar, J))
  tbar + (pm - tbar) * sqrt(sig2h / vfin)
}

# Triple-goal (Shen & Louis 1998): quantiles of Gbar at estimated ranks
est_gr <- function(fit) {
  theta <- fit$theta
  M <- nrow(theta)
  J <- ncol(theta)
  rmat <- t(apply(theta, 1L, rank, ties.method = "average"))
  Rbar <- colMeans(rmat) # posterior expected ranks
  Rhat <- rank(Rbar, ties.method = "first")
  Uk <- quantile(theta,
    probs = (2 * seq_len(J) - 1) / (2 * J),
    type = 1, names = FALSE
  ) # Gbar^{-1} at mid-quantiles
  Uk[Rhat]
}

# Finite-population losses (Lee et al. 2025, eq. 5-8)

loss_msel <- function(a, tau) mean((a - tau)^2)
loss_rmse <- function(a, tau) sqrt(loss_msel(a, tau))

loss_mselp <- function(a, tau) {
  J <- length(tau)
  A <- rank(a, ties.method = "average") / J
  R <- rank(tau, ties.method = "average") / J
  mean((A - R)^2)
}

# Exact integrated squared difference between the EDFs of a and tau
loss_isel <- function(a, tau) {
  z <- sort(c(a, tau))
  Fa <- findInterval(z, sort(a)) / length(a)
  Ft <- findInterval(z, sort(tau)) / length(tau)
  n <- length(z)
  sum((Fa[-n] - Ft[-n])^2 * diff(z))
}

all_losses <- function(a, tau) {
  c(rmse = loss_rmse(a, tau), mselp = loss_mselp(a, tau), isel = loss_isel(a, tau))
}

# Plotting helpers for exploratory figures

okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")

theme_ms <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold")
    )
}

save_fig <- function(plot, path, width, height) {
  ggsave(path, plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  if (!file.exists(path)) stop("figure not written: ", path)
  invisible(path)
}

# Posterior summaries of the reliability diagnostics
# For each posterior draw sigma^(m), compute S_j, I, H, J_eff; report
# posterior mean and central credible intervals (Proposition 3 favors the
# posterior functional over the plug-in).

posterior_diagnostics <- function(fit, v, probs = c(0.025, 0.5, 0.975)) {
  Smat <- outer(fit$sigma^2, v, function(s2, vv) s2 / (s2 + vv)) # M x J
  Lmat <- logit(Smat)
  Imat <- logistic(rowMeans(Lmat))
  Hvec <- vapply(
    seq_len(nrow(Smat)),
    function(m) mean(kl_bern(Imat[m], Smat[m, ])), numeric(1)
  )
  Jeff <- rowSums(Smat)^2 / rowSums(Smat^2)
  Vlog <- apply(Lmat, 1, function(l) mean((l - mean(l))^2))
  summarize <- function(x) c(mean = mean(x), quantile(x, probs, names = TRUE))
  list(
    draws  = list(I = Imat, H = Hvec, J_eff = Jeff, V = Vlog, S = Smat),
    I      = summarize(Imat),
    H      = summarize(Hvec),
    J_eff  = summarize(Jeff),
    V      = summarize(Vlog),
    S_mean = colMeans(Smat)
  )
}
