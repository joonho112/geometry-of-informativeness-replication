#!/usr/bin/env Rscript
# Posterior summaries for independent normal site estimates, conditional on v.
# The prior is flat in the overall mean tau and between-site SD sigma, not sigma^2.
# Source this file and use posterior_reference(y, v)$summary for posterior means.
# Running the script also compares adaptive quadrature with grid-based fits.

options(stringsAsFactors = FALSE)
softplus_stable <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# Compute the reliability profile at each sigma. The entries of v are
# sampling variances; sigma is an SD in the same units as the site effects.

profile_values <- function(sigma, v) {
  sigma <- as.numeric(sigma)
  J <- length(v)
  lv <- log(v)
  lgm <- mean(lv)
  a <- sigma^2
  m <- 2 * log(sigma) - lgm
  S <- outer(a, v, function(aa, vv) aa / (aa + vv))
  I <- plogis(m)
  H <- vapply(m, function(z) {
    if (!is.finite(z)) return(0)

    # Cancel the common linear term before subtracting at large positive logits.
    if (z >= 0) mean(softplus_stable(-z + lv - lgm)) - softplus_stable(-z)
    else mean(softplus_stable(z - lv + lgm)) - softplus_stable(z)
  }, numeric(1))
  H <- pmax(H, 0)
  r <- outer(a, v, function(aa, vv) 1 / (aa + vv))
  Jeff <- rowSums(r)^2 / rowSums(r^2)
  Jeff[is.infinite(a)] <- J
  Sbar <- rowMeans(S)
  Sbar[is.infinite(a)] <- 1
  cbind(I = I, H = H, Sbar = Sbar, Jeff = Jeff, a = a, sigma = sigma)
}

# Integrate out tau analytically, then integrate over sigma from zero to infinity.
# This routine requires J > 6 to also check the fourth posterior moment of sigma.
make_reference_kernel <- function(y, v, rel.tol = 1e-10, abs.tol = 1e-12) {
  stopifnot(
    length(y) == length(v), length(y) > 6, all(is.finite(y)),
    all(is.finite(v)), all(v > 0)
  )
  J <- length(y)
  scale <- sqrt(exp(mean(log(v))))
  logg <- function(s) vapply(s, function(z) {
    if (!is.finite(z)) return(-Inf)
    if (z > 1e100) return((1 - J) * log(z) - .5 * log(J))
    w <- 1 / (z * z + v)
    mu <- sum(w * y) / sum(w)
    -.5 * log(sum(w)) - .5 * sum(log(z * z + v)) - .5 * sum(w * (y - mu)^2)
  }, numeric(1))
  upper <- max(2 * sd(y), 3 * sqrt(median(v)), 3 * sqrt(max(0, var(y) - mean(v))), scale)
  opt <- optimize(function(s) -logg(s), c(0, upper))
  shift <- max(logg(c(0, opt$minimum, upper)))
  density_un <- function(s) exp(logg(s) - shift)
  quad <- function(fun = function(s) rep(1, length(s)), lo = 0, hi = Inf,
                   relative = rel.tol, absolute = abs.tol) {
    ff <- function(s) {
      d <- density_un(s)
      ans <- d * fun(s)
      ans[d == 0] <- 0
      ans
    }
    z <- integrate(ff, lo, hi,
      subdivisions = 2000L, rel.tol = relative,
      abs.tol = absolute, stop.on.error = TRUE
    )
    c(value = z$value, abs.error = z$abs.error)
  }

  # Split at sqrt(GM(v)) so integration finds a posterior concentrated near zero.
  zz1 <- quad(hi = scale)
  zz2 <- quad(lo = scale)
  Z <- unname(zz1['value'] + zz2['value'])
  Zerr <- zz1['abs.error'] + zz2['abs.error']
  expect <- function(fun) unname((quad(fun, hi = scale)['value'] + quad(fun, lo = scale)['value']) / Z)
  cdf <- function(s) vapply(s, function(z) {
    if (z <= 0) return(0)
    if (is.infinite(z)) return(1)
    ans <- if (z <= scale) quad(hi = z)['value'] else zz1['value'] + quad(lo = scale, hi = z)['value']
    min(1, max(0, unname(ans / Z)))
  }, numeric(1))

  # Bound the posterior integral of sigma^q above B, using
  # g(sigma) <= J^(-1/2) * sqrt(1 + max(v)/B^2) * sigma^(1-J).
  # The bound requires J > q + 2; it does not truncate the reference integral.
  tail_bound <- function(B, q = 0) {
    stopifnot(J > q + 2)
    exp(-shift - .5 * log(J) + .5 * log1p(max(v) / B^2) +
      (q + 2 - J) * log(B) - log(J - q - 2) - log(Z))
  }
  B <- upper
  while (tail_bound(B) > 1e-12 || tail_bound(B, 2) > 1e-10 * scale^2 ||
    tail_bound(B, 4) > 1e-10 * scale^4) B <- B * 2
  qsigma <- function(p) vapply(p, function(pp) {
    if (pp <= 0) return(0)
    if (pp >= 1) return(Inf)
    ub <- B
    while (cdf(ub) < pp) ub <- 2 * ub
    uniroot(function(s) cdf(s) - pp, c(0, ub), tol = 1e-11 * scale)$root
  }, numeric(1))
  tau_cond <- function(s) vapply(s, function(z) {
    w <- 1 / (z * z + v)
    sum(w * y) / sum(w)
  }, numeric(1))
  list(
    y = y, v = v, J = J, scale = scale, logg = logg, shift = shift, quad = quad,
    Z = Z, Zerr = unname(Zerr), expect = expect, cdf = cdf, qsigma = qsigma, B = B,
    tail_bound = tail_bound, tau_cond = tau_cond, density = function(s) density_un(s) / Z
  )
}

# H can rise and fall as sigma increases. Compute its CDF over all branches,
# then invert the CDF to find H quantiles. Check turning points on two search grids.
make_H_distribution <- function(v, cdf_sigma, q_sigma = NULL) {
  lv <- log(v)
  lgm <- mean(lv)
  V <- mean((lv - lgm)^2)
  if (V < 1e-25) return(list(
    q = function(p) rep(0, length(p)), cdf = function(h) as.numeric(h >= 0),
    roots = numeric(), scan_agrees = TRUE, max = 0
  ))
  from_m <- function(m) exp((m + lgm) / 2)
  hv <- function(m) profile_values(from_m(m), v)[, 'H']
  deriv <- function(m) {
    z <- profile_values(from_m(m), v)
    z[, 'Sbar'] - z[, 'I']
  }
  get_roots <- function(n) {
    xx <- seq(-40, 40, length.out = n)
    yy <- deriv(xx)

    # Avoid cancellation-induced sign changes in the extreme saturated tail.
    ids <- which(yy[-n] * yy[-1] < 0 & pmax(abs(yy[-n]), abs(yy[-1])) > 1e-12)
    if (!length(ids)) stop('No interior H turning point found')
    sort(unique(round(
      vapply(ids, function(i) uniroot(deriv, xx[c(i, i + 1)], tol = 1e-12)$root, numeric(1)),
      10
    )))
  }
  roots <- get_roots(1601)
  roots2 <- get_roots(3201)
  scan_ok <- length(roots) == length(roots2) && max(abs(roots - roots2)) < 1e-7
  if (!scan_ok) stop('H turning-point refinement mismatch')
  mm <- c(-80, roots, 80)
  hh <- c(0, hv(roots), 0)
  ss <- c(0, from_m(roots), Inf)
  cdf_H <- function(h) vapply(h, function(t) {
    if (t < 0) return(0)
    if (t >= max(hh)) return(1)
    ans <- 0
    for (k in seq_len(length(mm) - 1)) {
      h1 <- hh[k]
      h2 <- hh[k + 1]
      s1 <- ss[k]
      s2 <- ss[k + 1]
      if (t >= max(h1, h2)) ans <- ans + cdf_sigma(s2) - cdf_sigma(s1)
      else if (t > min(h1, h2)) {
        root <- uniroot(function(z) hv(z) - t, mm[c(k, k + 1)], tol = 1e-11)$root
        sr <- from_m(root)
        if (h2 > h1) ans <- ans + cdf_sigma(sr) - cdf_sigma(s1)
        else ans <- ans + cdf_sigma(s2) - cdf_sigma(sr)
      }
    }
    min(1, max(0, ans))
  }, numeric(1))
  quant <- function(p) vapply(p, function(pp) {
    if (pp <= 0) return(0)
    if (pp >= 1) return(max(hh))
    uniroot(function(h) cdf_H(h) - pp, c(0, max(hh)), tol = 1e-12)$root
  }, numeric(1))
  list(q = quant, cdf = cdf_H, roots = roots, scan_agrees = scan_ok, max = max(hh))
}

# Return posterior means and, if requested, quantiles for each diagnostic.
# These condition on the supplied site estimates y and sampling variances v.
posterior_reference <- function(y, v, quantiles = FALSE, independent_check = FALSE) {
  k <- make_reference_kernel(y, v)
  met <- c('I', 'H', 'Sbar', 'Jeff', 'a', 'sigma')
  mn <- sapply(met, function(nm) k$expect(function(s) profile_values(s, v)[, nm]))
  names(mn) <- paste0(met, '.mean')
  mu <- k$expect(k$tau_cond)
  gm <- exp(mean(log(v)))
  summary <- c(mn,
    tau.mean = mu, I_plugin = mn['a.mean'] / (mn['a.mean'] + gm),
    P_I_gt_half = 1 - k$cdf(sqrt(gm)), V = mean((log(v) - log(gm))^2)
  )
  names(summary)[names(summary) == 'I_plugin.a.mean'] <- 'I_plugin'
  hdist <- NULL
  if (quantiles) {
    probs <- c(.025, .5, .975)
    suffix <- c('q025', 'q500', 'q975')
    sq <- k$qsigma(probs)
    vv <- profile_values(sq, v)
    for (nm in c('I', 'Sbar', 'Jeff', 'a', 'sigma'))
      summary <- c(summary, setNames(vv[, nm], paste0(nm, '.', suffix)))
    hdist <- make_H_distribution(v, k$cdf)
    summary <- c(summary, setNames(hdist$q(probs), paste0('H.', suffix)))
  }
  checks <- c(
    normalization_relative_error = k$Zerr / k$Z,
    tail_bound_mass = k$tail_bound(k$B), tail_bound_a = k$tail_bound(k$B, 2),
    tail_bound_a2 = k$tail_bound(k$B, 4), tail_cutoff = k$B
  )
  if (independent_check) {
    # Check the integral on the log-sigma scale, including the Jacobian exp(t).
    lf <- function(t, fun) {
      z <- exp(pmin(t, 700))
      dens <- exp(k$logg(z) - k$shift + t)
      ans <- dens * fun(z)
      ans[dens == 0 | !is.finite(z)] <- 0
      ans
    }
    lint <- function(fun) {
      f <- function(t) lf(t, fun)
      integrate(f, -Inf, log(k$scale), rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 2000)$value +
        integrate(f, log(k$scale), Inf, rel.tol = 1e-10, abs.tol = 1e-12, subdivisions = 2000)$value
    }
    zl <- lint(function(s) rep(1, length(s)))
    ml <- sapply(met, function(nm) lint(function(s) profile_values(s, v)[, nm]) / zl)
    checks <- c(checks,
      normalization_logsigma_relative_difference = abs(zl - k$Z) / k$Z,
      setNames(ml - unname(mn), paste0(met, '.logsigma_minus_sigma'))
    )
  }
  list(summary = summary, kernel = k, H_distribution = hdist, checks = checks)
}

# Nodes and weights for Gauss-Legendre integration on [-1, 1].
gauss_legendre <- function(n) {
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i * i - 1)
  A <- matrix(0, n, n)
  A[cbind(i, i + 1)] <- b
  A[cbind(i + 1, i)] <- b
  e <- eigen(A, symmetric = TRUE)
  oo <- order(e$values)
  list(x = e$values[oo], w = 2 * e$vectors[1, oo]^2)
}

# Reconstruct the grid sampler, including uniform jitter and its floor.
# Compare this approximation with adaptive quadrature.
legacy_law <- function(y, v, n_grid = 250L, stored_grid = NULL, gauss_n = 24L) {
  if (is.null(stored_grid)) {
    b <- max(3 * sqrt(max(0, var(y) - mean(v))), 2 * sd(y), 3 * sqrt(median(v)), 1e-3)
    logg <- make_reference_kernel(y, v)$logg
    repeat {
      sg <- seq(0, b, length.out = n_grid)
      sg[1] <- b / (10 * n_grid)
      lp <- logg(sg)
      pr <- exp(lp - max(lp))
      pr <- pr / sum(pr)
      mass <- sum(pr[sg > .95 * b])
      safeguard <- b > 100 * sd(y)
      if (mass < 1e-4 || safeguard) break
      b <- 2 * b
    }
  } else {
    sg <- stored_grid$sigma
    pr <- stored_grid$prob
    b <- max(sg)
    n_grid <- length(sg)
    mass <- sum(pr[sg > .95 * b])
    safeguard <- b > 100 * sd(y)
  }
  delta <- sg[2] - sg[1]
  floor <- delta / 20
  lo <- sg - delta / 2
  hi <- sg + delta / 2
  lower <- pmax(lo, floor)
  pa <- pmax(0, pmin(1, (floor - lo) / delta)) * pr
  gl <- gauss_legendre(gauss_n)
  points <- as.vector(outer((hi + lower) / 2, rep(1, gauss_n)) + outer((hi - lower) / 2, gl$x))
  weights <- as.vector(outer(pr * (hi - lower) / (2 * delta), gl$w))
  points <- c(points, floor)
  weights <- c(weights, sum(pa))
  cdf <- function(s) vapply(s, function(t) {
    if (t < floor) return(0)
    sum(pr * pmin(1, pmax(0, (t - lo) / delta)))
  }, numeric(1))
  q <- function(p) vapply(p, function(pp) {
    if (pp <= cdf(floor)) return(floor)
    if (pp >= 1) return(max(hi))
    uniroot(function(s) cdf(s) - pp, c(floor, max(hi)), tol = 1e-12 * sqrt(exp(mean(log(v)))))$root
  }, numeric(1))
  list(
    points = points, weights = weights, cdf = cdf, q = q, sg = sg, pr = pr,
    b = b, delta = delta, floor = floor, top_mass = mass, safeguard = safeguard, n_grid = n_grid, y = y, v = v
  )
}

# Integrate diagnostics under the reconstructed grid-and-jitter distribution.
summarize_law <- function(law, v, quantiles = TRUE) {
  val <- profile_values(law$points, v)
  mn <- colSums(val * law$weights)
  sm <- setNames(mn, paste0(names(mn), '.mean'))
  gm <- exp(mean(log(v)))
  sm <- c(sm, I_plugin = unname(mn['a'] / (mn['a'] + gm)), P_I_gt_half = 1 - law$cdf(sqrt(gm)))
  hd <- NULL
  if (quantiles) {
    probs <- c(.025, .5, .975)
    suffix <- c('q025', 'q500', 'q975')
    vv <- profile_values(law$q(probs), v)
    for (nm in c('I', 'Sbar', 'Jeff', 'a', 'sigma'))
      sm <- c(sm, setNames(vv[, nm], paste0(nm, '.', suffix)))
    hd <- make_H_distribution(v, law$cdf)
    sm <- c(sm, setNames(hd$q(probs), paste0('H.', suffix)))
  }
  vv <- pmax(0, colSums(val^2 * law$weights) - mn^2)
  names(vv) <- colnames(val)
  cc <- conditional_site(law$points, law$y, v)
  tau_mean <- sum(law$weights * cc$mu)
  sm <- c(sm, tau.mean = tau_mean)
  vv <- c(vv, tau = sum(law$weights * (1 / cc$W + cc$mu^2)) - tau_mean^2)
  list(summary = sm, variances = vv, H_distribution = hd)
}

stat_draws <- function(sig, v) {
  z <- profile_values(sig, v)
  sm <- setNames(colMeans(z), paste0(colnames(z), '.mean'))
  for (nm in colnames(z)) sm <- c(
    sm,
    setNames(
      quantile(z[, nm], c(.025, .5, .975), names = FALSE),
      paste0(nm, c('.q025', '.q500', '.q975'))
    )
  )
  c(sm, I_plugin = mean(sig^2) / (mean(sig^2) + exp(mean(log(v)))), P_I_gt_half = mean(z[, 'I'] > .5))
}

metric_tol <- function(nm, ref, tols, v) {
  if (nm %in% names(tols$deterministic_absolute_tolerance))
    return(unname(tols$deterministic_absolute_tolerance[[nm]]))
  if (nm == 'a.mean') return(max(tols$Ea_tolerance$absolute, tols$Ea_tolerance$relative * abs(ref)))
  if (grepl('^sigma\\.', nm)) return(sqrt(exp(mean(log(v)))) * if (grepl('mean', nm)) .001 else .002)
  if (nm == 'tau.mean') return(.001 * sqrt(exp(mean(log(v)))))
  NA_real_
}

# Site means and variances conditional on sigma, after integrating out tau.
# The variance includes uncertainty in tau through the (1 - S)^2 / W term.
conditional_site <- function(sigma, y, v) {
  a <- sigma^2
  den <- outer(a, v, `+`)
  Wmat <- 1 / den
  W <- rowSums(Wmat)
  mu <- rowSums(sweep(Wmat, 2, y, `*`)) / W
  S <- a / den
  mm <- sweep(S, 2, y, `*`) + (1 - S) * mu
  vv <- a * (1 - S) + (1 - S)^2 / W
  list(mean = mm, variance = vv, mu = mu, W = W, S = S)
}

# Average first and second site moments over sigma to obtain marginal means
# and variances. Either adaptive quadrature or the grid distribution can be used.
site_moments <- function(y, v, reference = NULL, law = NULL) {
  J <- length(v)
  if (!is.null(reference)) {
    pm <- vapply(
      seq_len(J),
      function(j) reference$kernel$expect(function(s) conditional_site(s, y, v)$mean[, j]),
      numeric(1)
    )
    second <- vapply(seq_len(J), function(j) reference$kernel$expect(function(s) {
      cc <- conditional_site(s, y, v)
      cc$variance[, j] + cc$mean[, j]^2
    }), numeric(1))
  } else {
    cc <- conditional_site(law$points, y, v)
    pm <- colSums(cc$mean * law$weights)
    second <- colSums((cc$variance + cc$mean^2) * law$weights)
  }
  list(pm = pm, pv = pmax(0, second - pm^2))
}

# Tabulate the sigma CDF for fast sampling, refining until its error is small.
# The squared grid gives more points near zero.
inverse_cdf_cache <- function(reference, n = 65537L, max_error = 2e-8) {
  k <- reference$kernel
  repeat {
    sg <- k$B * seq(0, 1, length.out = n)^2
    d <- k$density(sg)
    cell <- diff(sg) * (head(d, -1) + tail(d, -1)) / 2
    mass <- sum(cell)
    cdfg <- c(0, cumsum(cell)) / mass
    cdfg[length(cdfg)] <- 1
    keep <- c(TRUE, diff(cdfg) > 0)
    cc <- cdfg[keep]
    xx <- sg[keep]
    inv <- function(p) approx(cc, xx, xout = p, rule = 2, ties = 'ordered')$y
    pp <- sort(unique(c(1e-6, .0001, .001, seq(.005, .995, length.out = 101), .999, .9999, 1 - 1e-6)))
    error <- max(abs(k$cdf(inv(pp)) - pp))
    if (error <= max_error) break
    if (n >= 262145L) stop('Inverse-CDF cache failed accuracy criterion')
    n <- 2L * n - 1L
  }

  # Trapezoidal weights support plots and integration over the stored grid.
  h <- diff(sg)
  weights <- d * c(h[1] / 2, (head(h, -1) + tail(h, -1)) / 2, tail(h, 1) / 2) / mass
  list(
    inverse = inv, sigma = sg, density = d, cdf = cdfg, weights = weights,
    checks = c(n = n, max_CDF_error = error, trapezoid_mass_minus_one = mass - 1),
    checked_probabilities = pp
  )
}

# Use quadrature for summaries and site moments, and inverse-CDF samples for
# posterior draws. Midpoint samples give stable density plots; they are not iid.
corrected_fit <- function(z, reference, seed, M = 20000L) {
  y <- z$y
  v <- z$v
  J <- length(v)
  ic <- inverse_cdf_cache(reference)
  set.seed(seed)
  sig <- ic$inverse(runif(M))
  cc <- conditional_site(sig, y, v)
  tau <- rnorm(M, cc$mu, sqrt(1 / cc$W))
  theta_mean <- sweep(cc$S, 2, y, `*`) + (1 - cc$S) * tau
  theta_sd <- sqrt(sig^2 * (1 - cc$S))
  theta <- theta_mean + matrix(rnorm(M * J), M, J) * theta_sd
  moments <- site_moments(y, v, reference = reference)
  fit <- list(
    theta = theta, sigma = sig, tau = tau, pm = moments$pm, pv = moments$pv,
    grid = list(sigma = ic$sigma, prob = ic$weights)
  )
  sig_plot <- ic$inverse((seq_len(M) - .5) / M)
  make_draws <- function(s) {
    dd <- profile_values(s, v)
    Sm <- outer(s^2, v, function(a, vv) a / (a + vv))
    V <- mean((log(v) - mean(log(v)))^2)
    list(
      I = dd[, 'I'],
      H = dd[, 'H'],
      J_eff = dd[, 'Jeff'],
      S_bar = dd[, 'Sbar'],
      V = rep(V, length(s)),
      S = Sm
    )
  }
  rs <- reference$summary
  sumvec <- function(nm) setNames(
    unname(rs[paste0(nm, c('.mean', '.q025', '.q500', '.q975'))]),
    c('mean', '2.5%', '50%', '97.5%')
  )
  Smean <- vapply(
    seq_len(J),
    function(j) reference$kernel$expect(function(s) s * s / (s * s + v[j])),
    numeric(1)
  )
  pdg <- list(
    draws = make_draws(sig), I = sumvec('I'), H = sumvec('H'), J_eff = sumvec('Jeff'),
    S_bar = sumvec('Sbar'), V = setNames(rep(rs['V'], 4), c('mean', '2.5%', '50%', '97.5%')), S_mean = Smean,
    plot_draws = make_draws(sig_plot)
  )
  list(
    label = z$label,
    fit = fit,
    pdg = pdg,
    v = v,
    y = y,
    gm = exp(mean(log(v))),
    I_plug = unname(rs['I_plugin']),
    p_I_gt_half = unname(rs['P_I_gt_half']), sigma_summ = sumvec('sigma'),
    reference_summary = rs, numerics = list(
      method = 'untruncated flat-mean/SD posterior; deterministic adaptive quadrature summaries; cached inverse-CDF iid hyperparameter and conditional effect draws',
      seed = seed, M = M, inverse_cdf_checks = ic$checks,
      checked_probabilities = ic$checked_probabilities,
      plot_draws_method = 'deterministic midpoint inverse-CDF samples for stable density displays; not iid',
      pm_pv_method = 'conditional hypermean analytically integrated; marginal first/second site moments by adaptive quadrature',
      tau_mean = 'Use reference_summary[tau.mean] for tables; fit$tau contains Monte Carlo draws'
    )
  )
}

run_posterior_checks <- function() {
  script_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  script_path <- normalizePath(sub('^--file=', '', script_arg[1]))
  oldwd <- setwd(dirname(dirname(dirname(script_path))))
  on.exit(setwd(oldwd), add = TRUE)
  start <- Sys.time()
  out <- 'analysis/outputs/numerical-stability'
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  tolerance_file <- 'analysis/inputs/numerical-tolerances.json'
  stopifnot(file.exists(tolerance_file))
  file.copy(tolerance_file, file.path(out, 'tolerances.json'), overwrite = TRUE)
  tols <- jsonlite::fromJSON(tolerance_file, simplifyVector = TRUE)
  tols$deterministic_absolute_tolerance <- unlist(tols$deterministic_absolute_tolerance)
  x <- readRDS('analysis/inputs/ocrs_fits.rds')
  fits <- c(x$fits, list(logodds = x$logodds_fit))
  result <- list()
  comparison <- list()
  mc <- list()
  qa <- list()
  tails <- list()
  refined <- list()
  corrected <- list()
  cache_checks <- list()
  decomposition <- list()
  for (id in names(fits)) {
    cat('OCRS', id, 'start', format(Sys.time()), '\n')
    flush.console()
    z <- fits[[id]]
    r <- posterior_reference(z$y, z$v, TRUE, TRUE)
    law <- legacy_law(z$y, z$v, stored_grid = z$fit$grid)
    ls <- summarize_law(law, z$v)
    ss <- c(stat_draws(z$fit$sigma, z$v), tau.mean = mean(z$fit$tau))
    law48 <- legacy_law(z$y, z$v, stored_grid = z$fit$grid, gauss_n = 48)
    ls48 <- summarize_law(law48, z$v, FALSE)
    nms <- intersect(names(ls$summary), names(r$summary))
    tol <- sapply(nms, function(nm) metric_tol(nm, r$summary[nm], tols, z$v))
    comparison[[id]] <- data.frame(
      fit = id, metric = nms, reference = unname(r$summary[nms]),
      legacy_law = unname(ls$summary[nms]), saved_draws = unname(ss[nms]),
      grid_bias = unname(ls$summary[nms] - r$summary[nms]),
      saved_minus_legacy_law = unname(ss[nms] - ls$summary[nms]), tolerance = tol,
      deterministic_pass = ifelse(is.na(tol), NA, abs(ls$summary[nms] - r$summary[nms]) <= tol)
    )

    # Separate the effects of truncating sigma, discretizing it, and adding jitter.
    fields <- c('I', 'H', 'Sbar', 'Jeff', 'a', 'sigma')
    fb <- r$kernel$cdf(law$b)
    tr <- sapply(
      fields,
      function(nm) unname(r$kernel$quad(function(s) profile_values(s, z$v)[, nm], hi = law$b)['value'] / r$kernel$Z / fb)
    )
    pt <- colSums(profile_values(law$sg, z$v) * law$pr)
    tr <- c(
      setNames(tr, paste0(fields, '.mean')),
      I_plugin = tr['a'] / (tr['a'] + z$gm),
      P_I_gt_half = max(0, (fb - r$kernel$cdf(sqrt(z$gm))) / fb)
    )
    names(tr)[names(tr) == 'I_plugin.a'] <- 'I_plugin'
    pt <- c(
      setNames(pt, paste0(names(pt), '.mean')),
      I_plugin = pt['a'] / (pt['a'] + z$gm),
      P_I_gt_half = sum(law$pr[law$sg > sqrt(z$gm)])
    )
    names(pt)[names(pt) == 'I_plugin.a'] <- 'I_plugin'
    dn <- intersect(names(tr), names(r$summary))
    decomposition[[id]] <- data.frame(
      fit = id,
      metric = dn,
      reference = unname(r$summary[dn]),
      truncated_reference = unname(tr[dn]),
      node_only = unname(pt[dn]),
      jitter_law = unname(ls$summary[dn]),
      tail_component = unname(tr[dn] - r$summary[dn]),
      grid_component = unname(pt[dn] - tr[dn]),
      jitter_component = unname(ls$summary[dn] - pt[dn])
    )

    # Monte Carlo error here concerns draws from the grid sampler.
    # It does not measure participant-sampling uncertainty or grid bias.
    M <- length(z$fit$sigma)
    mcr <- list()
    for (nm in nms) {
      est <- ss[nm]
      target <- ls$summary[nm]
      se <- NA_real_
      low <- high <- NA_real_
      method <- ''
      if (grepl('\\.mean$', nm)) {
        field <- sub('\\.mean$', '', nm)
        se <- sqrt(ls$variances[field] / M)
        method <- 'iid mean MCSE under exact legacy jitter law'
        low <- target - 3.5 * se
        high <- target + 3.5 * se
      } else if (nm == 'P_I_gt_half') {
        se <- sqrt(target * (1 - target) / M)
        low <- target - 3.5 * se
        high <- target + 3.5 * se
        method <- 'Bernoulli MCSE'
      } else if (nm == 'I_plugin') {
        aa <- ls$summary['a.mean']
        gm <- exp(mean(log(z$v)))
        se <- gm / (aa + gm)^2 * sqrt(ls$variances['a'] / M)
        low <- target - 3.5 * se
        high <- target + 3.5 * se
        method <- 'delta-method MCSE for plugin of sample mean a'
      } else if (grepl('\\.q', nm)) {
        pp <- as.numeric(sub('.*\\.q', '', nm)) / 1000
        field <- sub('\\.q.*', '', nm)
        k <- floor((M - 1) * pp + 1)
        pl <- qbeta(.0005, k, M + 1 - k)
        pu <- qbeta(.9995, k + 1, M - k)
        getq <- function(p) if (field == 'H') ls$H_distribution$q(p) else profile_values(
          law$q(p),
          z$v
        )[, field]
        low <- getq(pl)
        high <- getq(pu)
        pl1 <- qbeta(pnorm(-1), k, M + 1 - k)
        pu1 <- qbeta(pnorm(1), k + 1, M - k)
        se <- (getq(pu1) - getq(pl1)) / 2
        method <- 'exact order-statistic/binomial inversion bounds, 99.9%; interpolated-quantile envelope'
      }
      mcr[[nm]] <- data.frame(
        fit = id, metric = nm, M = M, saved = unname(est), legacy_target = unname(target),
        mcse = unname(se),
        mc_error = unname(est - target),
        z = ifelse(se > 0, unname((est - target) / se), NA),
        sampling_low = unname(low),
        sampling_high = unname(high),
        within_sampling_interval = unname(est >= low && est <= high),
        method = method
      )
    }
    mc[[id]] <- do.call(rbind, mcr)
    qa[[id]] <- data.frame(fit = id, check = names(r$checks), value = unname(r$checks))
    qa[[paste0(id, 'gauss')]] <- data.frame(
      fit = id,
      check = 'legacy_quadrature_24_vs_48_max_abs',
      value = max(abs(ls48$summary - ls$summary[names(ls48$summary)]))
    )

    # Compare omitted tail integrals with their analytic upper bounds.
    for (q in c(0, 2, 4)) {
      fun <- function(s) s^q
      actual <- unname(r$kernel$quad(fun, lo = law$b)['value'] / r$kernel$Z)
      full <- r$kernel$expect(fun)
      tails[[paste(id, q)]] <- data.frame(
        fit = id, power = q, legacy_cutoff = law$b,
        legacy_top5pct_mass = law$top_mass, safeguard_stopped = law$safeguard,
        omitted_tail_integral = actual, tail_analytic_upper_bound = r$kernel$tail_bound(law$b, q),
        full_moment = full, omitted_relative = actual / full,
        reference_cutoff = r$kernel$B, reference_tail_bound = r$kernel$tail_bound(r$kernel$B, q)
      )
    }
    for (ng in c(1000L, 4000L)) {
      ll <- legacy_law(z$y, z$v, n_grid = ng)
      zz <- summarize_law(ll, z$v, FALSE)
      nn <- intersect(names(zz$summary), names(r$summary))
      refined[[paste(id, ng)]] <- data.frame(
        fit = id, n_grid = ng, metric = nn, value = unname(zz$summary[nn]),
        reference = unname(r$summary[nn]), difference = unname(zz$summary[nn] - r$summary[nn])
      )
    }

    # Save numeric results; integration functions stay in memory.
    sg <- seq(0, r$kernel$B, length.out = 4001)
    result[[id]] <- list(
      label = z$label, y = z$y, v = z$v, summary = r$summary, checks = r$checks,
      posterior_curve = data.frame(sigma = sg, density = r$kernel$density(sg)),
      H_critical_logits = r$H_distribution$roots, legacy_grid = z$fit$grid,
      legacy_law_summary = ls$summary, saved_draw_summary = ss
    )
    corrected[[id]] <- corrected_fit(z, r, seed = 2026101400L + match(id, names(fits)))
    cache_checks[[id]] <- data.frame(
      fit = id,
      check = names(corrected[[id]]$numerics$inverse_cdf_checks),
      value = unname(corrected[[id]]$numerics$inverse_cdf_checks)
    )
    cat(
      'OCRS', id, 'done; max abs I/H mean difference',
      max(abs(ls$summary[c('I.mean', 'H.mean')] - r$summary[c('I.mean', 'H.mean')])), '\n'
    )
    flush.console()
  }
  cmp <- do.call(rbind, comparison)
  mct <- do.call(rbind, mc)
  qt <- do.call(rbind, qa)
  tt <- do.call(rbind, tails)
  rt <- do.call(rbind, refined)
  write.csv(cmp, file.path(out, 'ocrs-comparison.csv'), row.names = FALSE)
  stopifnot(all(is.finite(mct$mcse)), all(!is.na(mct$within_sampling_interval)))
  write.csv(mct, file.path(out, 'ocrs-monte-carlo.csv'), row.names = FALSE)
  write.csv(
    do.call(rbind, decomposition),
    file.path(out, 'grid-jitter-tail-decomposition.csv'),
    row.names = FALSE
  )
  write.csv(qt, file.path(out, 'reference-quadrature-checks.csv'), row.names = FALSE)
  write.csv(tt, file.path(out, 'tail-contributions.csv'), row.names = FALSE)
  write.csv(rt, file.path(out, 'refined-grid-comparison.csv'), row.names = FALSE)
  saveRDS(result, file.path(out, 'ocrs-reference-summaries.rds'))

  # Retain the site summaries and plotting values used by the figures.
  # corrected_fit() can recreate the draws. This object has no participant records.
  compact <- function(z) {
    z$fit <- z$fit[c('pm', 'pv')]
    z$pdg$draws <- NULL
    z$pdg$plot_draws <- z$pdg$plot_draws['I']
    z
  }
  saveRDS(
    list(
      fits = lapply(corrected[names(x$fits)], compact),
      logodds_fit = compact(corrected$logodds),
      sensitivity_status = 'Quadrature-based sensitivity results are in ../ocrs-audit/*reference*.csv.'
    ),
    file.path(out, 'ocrs_fits_corrected.rds'),
    compress = 'xz'
  )
  write.csv(do.call(rbind, cache_checks), file.path(out, 'inverse-cdf-cache-checks.csv'), row.names = FALSE)
  write.csv(
    do.call(
      rbind,
      lapply(
        names(result),
        function(id) data.frame(
          fit = id,
          metric = names(result[[id]]$summary),
          value = unname(result[[id]]$summary)
        )
      )
    ),
    file.path(out, 'ocrs-corrected-summaries.csv'),
    row.names = FALSE
  )

  # Compare the two numerical methods in six simulated examples.
  # These are sensitivity checks, not a rerun of the full simulation study.
  legacy <- new.env(parent = globalenv())
  sys.source('analysis/scripts/legacy/00_utils.R', envir = legacy)
  cases <- tols$selected_simulation_checks$cases
  set.seed(tols$selected_simulation_checks$seed)
  sim <- list()
  sim_inputs <- list()
  sim_losses <- list()
  sim_estimators <- list()
  for (i in seq_len(nrow(cases))) {
    cc <- cases[i, ]
    a0 <- .2^2
    v <- legacy$make_profile(cc$shape, cc$J, a0 * (1 - cc$I) / cc$I, cc$V)
    theta <- rnorm(cc$J, 0, .2)
    y <- theta + rnorm(cc$J, 0, sqrt(v))
    case_rng <- .Random.seed
    cat('Boundary case', i, '\n')
    flush.console()
    rr <- posterior_reference(y, v, TRUE, TRUE)
    ll <- legacy_law(y, v)
    zz <- summarize_law(ll, v)
    nn <- intersect(names(rr$summary), names(zz$summary))
    tol <- sapply(nn, function(nm) metric_tol(nm, rr$summary[nm], tols, v))
    sim[[i]] <- data.frame(
      case = i, J = cc$J, I_true = cc$I, V = cc$V, shape = cc$shape, metric = nn,
      reference = unname(rr$summary[nn]), legacy_law = unname(zz$summary[nn]),
      difference = unname(zz$summary[nn] - rr$summary[nn]), tolerance = tol,
      pass = ifelse(is.na(tol), NA, abs(zz$summary[nn] - rr$summary[nn]) <= tol)
    )
    mr <- site_moments(y, v, reference = rr)
    ml <- site_moments(y, v, law = ll)
    ff <- legacy$fit_rubin_grid(y, v, M = 2000L, seed = 2026101500L + i)

    # Restore the seed below so fitting does not change the next simulated dataset.
    aa <- list(
      reference_PM = mr$pm, reference_CB = legacy$est_cb(mr),
      legacy_limit_PM = ml$pm, legacy_limit_CB = legacy$est_cb(ml),
      legacy_M2000_PM = ff$pm, legacy_M2000_CB = legacy$est_cb(ff)
    )
    sim_losses[[i]] <- do.call(rbind, lapply(names(aa), function(nm) {
      ls0 <- legacy$all_losses(aa[[nm]], theta)
      data.frame(
        case = i,
        J = cc$J,
        I_true = cc$I,
        V = cc$V,
        shape = cc$shape,
        estimator = nm,
        rmse = ls0['rmse'],
        rank_mse = ls0['mselp'],
        isel = ls0['isel']
      )
    }))
    sim_estimators[[i]] <- list(actions = aa, reference_moments = mr, legacy_limit_moments = ml)
    sim_inputs[[i]] <- list(
      spec = cc,
      y = y,
      v = v,
      theta_true = theta,
      reference = rr$summary,
      checks = rr$checks
    )
    assign('.Random.seed', case_rng, envir = .GlobalEnv)
  }
  write.csv(do.call(rbind, sim), file.path(out, 'selected-boundary-cases.csv'), row.names = FALSE)
  saveRDS(sim_inputs, file.path(out, 'selected-boundary-inputs.rds'))
  loss_table <- do.call(rbind, sim_losses)
  write.csv(loss_table, file.path(out, 'selected-boundary-losses.csv'), row.names = FALSE)
  lr <- subset(loss_table, grepl('^reference_', estimator))
  lr$kind <- sub('reference_', '', lr$estimator)
  ll <- subset(loss_table, grepl('^legacy_limit_', estimator))
  ll$kind <- sub('legacy_limit_', '', ll$estimator)
  ld <- merge(ll, lr, by = c('case', 'kind'), suffixes = c('_legacy', '_reference'))
  for (nm in c('rmse', 'rank_mse', 'isel')) {
    ld[[paste0(nm, '_difference')]] <- ld[[paste0(nm, '_legacy')]] - ld[[paste0(nm, '_reference')]]
    ld[[paste0(nm, '_relative')]] <- ld[[paste0(nm, '_difference')]] / ld[[paste0(nm, '_reference')]]
  }
  write.csv(ld, file.path(out, 'selected-boundary-loss-differences.csv'), row.names = FALSE)
  saveRDS(sim_estimators, file.path(out, 'selected-boundary-estimators.rds'))
  main_metrics <- names(tols$deterministic_absolute_tolerance)
  fail <- subset(cmp, !is.na(deterministic_pass) & !deterministic_pass)
  receipt <- list(
    started = format(start, tz = 'UTC', usetz = TRUE),
    finished = format(Sys.time(), tz = 'UTC', usetz = TRUE),
    elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = 'secs')), R = R.version.string,
    input_md5 = as.character(tools::md5sum('analysis/inputs/ocrs_fits.rds')),
    script_md5 = as.character(tools::md5sum('analysis/scripts/05_posterior_stability.R')),
    tolerance_md5 = as.character(tools::md5sum(tolerance_file)),
    grid_tolerance_exceeded = nrow(fail) > 0, deterministic_failures = fail,
    monte_carlo_validation = list(
      rows = nrow(mct),
      nonfinite_mcse = sum(!is.finite(mct$mcse)),
      missing_interval_pass = sum(is.na(mct$within_sampling_interval)),
      outside_sampling_intervals = sum(!mct$within_sampling_interval)
    ),
    boundary_scope = 'six fresh fixed-seed Gaussian examples; diagnostic and PM/CB action/loss sensitivity only, not all original 702 cells or 500 replications',
    source_scope = 'saved fixed y/v; no rawdata first-stage fit; six fresh fixed-seed Gaussian scope checks',
    reference_method = 'physical sigma integration split at sqrt(GMv), independent log-sigma integration, analytic entire-tail moment bound',
    legacy_method = 'exact normalized-grid plus pmax-uniform-jitter mixture, Gaussian quadrature 24 vs48 checks',
    corrected_output = 'ocrs-reference-summaries.rds, ocrs-corrected-summaries.csv, and ocrs_fits_corrected.rds; same likelihood/prior/y/v; cached corrected iid draws plus deterministic plot samples; legacy input fit unchanged'
  )
  jsonlite::write_json(receipt,
    file.path(out, 'run-receipt.json'),
    pretty = TRUE,
    auto_unbox = TRUE,
    digits = 16,
    na = 'null'
  )
  capture.output(sessionInfo(), file = file.path(out, 'session-info.txt'))
  cat('Posterior checks finished; grid comparison exceeds tolerance:',
      receipt$grid_tolerance_exceeded, '; runtime', receipt$elapsed_seconds, 'seconds\n')
}

if (sys.nframe() == 0L) run_posterior_checks()
