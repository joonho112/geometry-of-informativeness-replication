# Check effective-site-count bounds and monotonicity numerically (OSM D.3).
# These finite-grid checks supplement the analytic proof.
# run_historical.R --arm=theory runs this in a separate output directory.

source("00_utils.R")
variances <- list(
  homogeneous = rep(.04, 25),
  even = profile_even(25, .04, 2),
  bimodal = profile_bimodal(50, .04, 2),
  concentrated = profile_concentrated(100, .04, 2),
  concentrated_noisy = profile_concentrated_noisy(100, .04, 2)
)
rows <- lapply(names(variances), function(nm) {
  v <- variances[[nm]]
  t <- exp(seq(log(1e-8), log(1e6), length.out = 501))
  count <- vapply(t, function(tt) J_eff(site_reliability(tt, v)), numeric(1))
  lower <- sum(1 / v)^2 / sum((1 / v)^2)
  stopifnot(
    all(diff(count) >= -1e-10), all(count >= lower - 1e-10),
    all(count <= length(v) + 1e-10)
  )
  # Independent finite differences check the derivative from the OSM proof.
  t0 <- .04
  h <- 1e-5 * t0
  q <- 1 / (t0 + v)
  analytic <- 2 * (sum(q)^2 / sum(q^2)) *
    (sum(q^3) / sum(q^2) - sum(q^2) / sum(q))
  finite <- (J_eff(site_reliability(t0 + h, v)) -
    J_eff(site_reliability(t0 - h, v))) / (2 * h)
  stopifnot(abs(analytic - finite) <= 1e-5 * max(1, abs(analytic)))
  data.frame(
    profile = nm, J = length(v), lower_limit = lower,
    minimum_increment = min(diff(count)),
    maximum_count = max(count), derivative = analytic,
    derivative_finite_difference = finite
  )
})
out <- do.call(rbind, rows)
stopifnot(
  max(abs(out$lower_limit[out$profile == "homogeneous"] - 25)) < 1e-10,
  max(abs(out$maximum_count[out$profile == "homogeneous"] - 25)) < 1e-10
)
write.csv(out, "outputs/effective-count-checks.csv", row.names = FALSE)
cat("Effective-count checks: 5 profiles x 501 positive-variance levels; PASS\n")
