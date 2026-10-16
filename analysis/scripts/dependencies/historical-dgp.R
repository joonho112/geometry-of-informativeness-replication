# Effect generators from multisiteDGP 0.1.0 used in the simulation study.
# These functions support only the no-covariate calls used in this study.
# See LICENSE-multisiteDGP.txt and verification/simulation-dgp-provenance.json.
# They return data.frames rather than tibbles. For the manuscript parameters,
# effect values and random-number draws match that package version.

.historical_effect_frame <- function(z, J, tau, sigma_tau) {
  data.frame(
    site_index = seq_len(J), z_j = z,
    tau_j = tau + 0 + sigma_tau * z
  )
}

gen_effects_gaussian <- function(J, tau = 0, sigma_tau = 0.20) {
  .historical_effect_frame(stats::rnorm(J), J, tau, sigma_tau)
}

gen_effects_mixture <- function(J, tau = 0, sigma_tau = 0.20,
                                delta, eps, ups) {
  component <- as.integer(ifelse(stats::runif(J) < (1 - eps), 1L, 2L))
  scale <- sqrt((1 - eps) + eps * ups^2 + eps * (1 - eps) * delta^2)
  component_1 <- stats::rnorm(J,
    mean = -eps * delta / scale,
    sd = sqrt(1 / scale^2)
  )
  component_2 <- stats::rnorm(J,
    mean = (1 - eps) * delta / scale,
    sd = sqrt(ups^2 / scale^2)
  )
  z <- ifelse(component == 1L, component_1, component_2)
  out <- .historical_effect_frame(z, J, tau, sigma_tau)
  out$latent_component <- component
  out
}

gen_effects_ald <- function(J, tau = 0, sigma_tau = 0.20, rho) {
  kappa <- sqrt(rho / (1 - rho))
  lambda <- sqrt(2 / (rho * (1 - rho)))
  # Use the two-uniform formula from the MIT-licensed multisiteDGP source,
  # preserving the operation order of the LaplacesDemon call.
  u1 <- stats::runif(J)
  u2 <- stats::runif(J)
  x <- lambda * log(u1^kappa / u2^(1 / kappa)) / sqrt(2)
  mu_x <- (1 - 2 * rho) / (rho * (1 - rho))
  var_x <- (1 - 2 * rho + 2 * rho^2) / (rho^2 * (1 - rho)^2)
  z <- (x - mu_x) / sqrt(var_x)
  .historical_effect_frame(z, J, tau, sigma_tau)
}
