# Numerical checks of selected identities and approximations.
# These supplement the proofs; they do not establish every manuscript claim.
# Outputs: verification_report.md and key_numbers.csv in outputs/.
# Run through analysis/scripts/run_historical.R --arm=theory.

source("00_utils.R")
set.seed(20260710)

report <- character(0)
keynum <- list()
note <- function(...) {
  line <- paste0(...)
  report <<- c(report, line)
  cat(line, "\n")
}
ANY_FAIL <- FALSE
check <- function(label, cond, detail = "") {
  status <- if (isTRUE(cond)) "PASS" else "FAIL"
  if (!isTRUE(cond)) ANY_FAIL <<- TRUE
  note(sprintf(
    "- [%s] %s%s", status, label,
    if (nzchar(detail)) paste0(" -- ", detail) else ""
  ))
  invisible(cond)
}

note("# Verification report -- Geometry of Informativeness (", format(Sys.Date()), ")\n")

note("## 1. Theorem 1: KL-centroid characterization of I")

sigma2 <- 0.04
v <- exp(rnorm(40, log(0.02), 0.9))
S <- site_reliability(sigma2, v)
I_formula <- informativeness_I(sigma2, v)
I_centroid <- kl_centroid(S)
check(
  "logit-mean form equals geometric-mean form",
  abs(I_formula - I_centroid) < 1e-12,
  sprintf("I = %.10f vs %.10f", I_formula, I_centroid)
)

# argmin verification over a fine grid
pgrid <- seq(1e-4, 1 - 1e-4, length.out = 200001)
obj <- vapply(pgrid, function(p) mean(kl_bern(p, S)), numeric(1))
p_star <- pgrid[which.min(obj)]
check(
  "I minimizes mean KL(Ber(p) || Ber(S_j)) over p-grid",
  abs(p_star - I_centroid) < 1e-4,
  sprintf("grid argmin %.6f vs I %.6f", p_star, I_centroid)
)

# Check convexity along the probability grid.
d2 <- diff(diff(obj[seq(1, length(obj), by = 100)]))
check(
  "objective convex along p-grid (second differences >= 0 up to fp noise)",
  min(d2) > -1e-12
)

note("\n## 2. Theorem 2: constant-manifold Pythagorean decomposition")

I <- I_centroid
H <- H_KL(S, I)
dec_err <- vapply(c(0.05, 0.2, 0.5, 0.77, 0.95), function(p) {
  lhs <- mean(kl_bern(p, S))
  rhs <- kl_bern(p, I) + H
  abs(lhs - rhs)
}, numeric(1))
check(
  "decomposition holds for all tested p", max(dec_err) < 1e-12,
  sprintf("max |lhs - rhs| = %.2e", max(dec_err))
)

# Check failure of the decomposition at a non-centroid center.
Ibad <- plogis(qlogis(I) + 0.3)
Hbad <- mean(kl_bern(Ibad, S))
bad_err <- abs(mean(kl_bern(0.3, S)) - (kl_bern(0.3, Ibad) + Hbad))
check(
  "decomposition fails for a non-centroid center (iff direction)",
  bad_err > 1e-4, sprintf("violation %.2e", bad_err)
)

check(
  "H_KL >= 0 and H_KL = 0 iff constant profile",
  H > 0 && abs(H_KL(rep(0.4, 10))) < 1e-15
)

note("\n## 3. Proposition 2: local structure and the risk link")

# (a) H ~ (1/2) I(1-I) V ; (b) S_bar - I ~ (1-2I) H
res <- do.call(rbind, lapply(c(0.05, 0.1, 0.2, 0.4), function(Vtar) {
  do.call(rbind, lapply(c(0.15, 0.3, 0.5, 0.65), function(Itar) {
    gm <- sigma2 * (1 - Itar) / Itar
    vv <- profile_even(60, gm, Vtar)
    SS <- site_reliability(sigma2, vv)
    II <- kl_centroid(SS)
    HH <- H_KL(SS, II)
    VV <- V_logit(SS)
    data.frame(Itar, Vtar,
      I = II, H = HH, V = VV,
      H_approx = 0.5 * II * (1 - II) * VV,
      gap = mean(SS) - II,
      gap_approx = (1 - 2 * II) * HH
    )
  }))
}))
rel_err_H <- with(res, abs(H_approx - H) / H)
abs_err_g <- with(res, abs(gap_approx - gap))
check(
  "H_KL = (1/2) I(1-I) V to second order (rel. err < 12% at V <= 0.4)",
  max(rel_err_H) < 0.12,
  sprintf("max rel err %.3f at V = %.2f", max(rel_err_H), res$Vtar[which.max(rel_err_H)])
)
check(
  "S_bar - I = (1-2I) H to second order (abs err < 6e-3 at V <= 0.4)",
  max(abs_err_g) < 6e-3, sprintf("max abs err %.1e", max(abs_err_g))
)
keynum$prop2_max_relerr_H <- max(rel_err_H)

# shrinking-V convergence: ratios must -> 1
conv <- vapply(c(0.4, 0.2, 0.1, 0.05, 0.025), function(Vtar) {
  gm <- sigma2 * (1 - 0.3) / 0.3
  vv <- profile_even(60, gm, Vtar)
  SS <- site_reliability(sigma2, vv)
  II <- kl_centroid(SS)
  H_KL(SS, II) / (0.5 * II * (1 - II) * V_logit(SS))
}, numeric(1))
check(
  "H / [(1/2) I(1-I) V] -> 1 as V -> 0 (second-order exactness)",
  abs(conv[length(conv)] - 1) < 0.01,
  paste0("ratios: ", paste(sprintf("%.4f", conv), collapse = ", "))
)

# Average PM risk with known hyperparameters equals sigma2(1 - S_bar).
J <- 50
gm <- sigma2 * (1 - 0.3) / 0.3
vv <- profile_bimodal(J, gm, 1.0)
SS <- site_reliability(sigma2, vv)
mc <- replicate(4000, {
  tau_j <- rnorm(J, 0, sqrt(sigma2))
  y <- rnorm(J, tau_j, sqrt(vv))
  pm <- SS * y # known tau = 0, sigma2
  mean((pm - tau_j)^2)
})
risk_mc <- mean(mc)
risk_th <- sigma2 * (1 - mean(SS))
check(
  "oracle-PM average risk equals sigma2 (1 - S_bar)",
  abs(risk_mc - risk_th) / risk_th < 0.02,
  sprintf("MC %.6f vs theory %.6f", risk_mc, risk_th)
)

# sign flip of the dispersion effect at I = 1/2 (oracle PM risk, fixed I)
flip <- vapply(c(0.3, 0.7), function(Itar) {
  gm <- sigma2 * (1 - Itar) / Itar
  r0 <- sigma2 * (1 - mean(site_reliability(sigma2, profile_even(J, gm, 0))))
  r1 <- sigma2 * (1 - mean(site_reliability(sigma2, profile_even(J, gm, 1.5))))
  r1 - r0
}, numeric(1))
check(
  "dispersion lowers oracle-PM risk at I = .3 and raises it at I = .7",
  flip[1] < 0 && flip[2] > 0,
  sprintf("delta risk: %.5f (I=.3), %.5f (I=.7)", flip[1], flip[2])
)
keynum$risk_flip_low <- flip[1]
keynum$risk_flip_high <- flip[2]

note("\n## 4. Proposition 4: design corollary")

kappa <- 4
N <- 2000
J <- 40
I_of_n <- function(n) informativeness_I(sigma2, kappa / n)
n_eq <- rep(N / J, J)
worse <- replicate(2000, {
  w <- as.numeric(rmultinom(1, N - 5 * J, rep(1, J))) + 5 # random allocations, n_j >= 5
  I_of_n(w)
})
check(
  "equal allocation beats 2000 random allocations",
  all(worse <= I_of_n(n_eq) + 1e-12),
  sprintf("I_eq = %.5f, max random = %.5f", I_of_n(n_eq), max(worse))
)

# efficiency ratio GM/AM ~ exp(-CV^2/2) for Gamma site sizes
eff <- vapply(c(0.25, 0.5, 0.75), function(cv) {
  n <- rgamma(2e5, shape = 1 / cv^2, scale = 60 * cv^2)
  exp(mean(log(n))) / mean(n)
}, numeric(1))
eff_th <- exp(-c(0.25, 0.5, 0.75)^2 / 2)
check(
  "GM/AM efficiency ratio approx exp(-CV^2/2) (Gamma sizes)",
  max(abs(eff - eff_th)) < 0.03,
  paste0(sprintf(
    "emp %.3f/%.3f/%.3f vs th %.3f/%.3f/%.3f",
    eff[1], eff[2], eff[3], eff_th[1], eff_th[2], eff_th[3]
  ))
)

note("\n## 5. Proposition 3: plug-in overstates informativeness (Jensen)")

J <- 25
sigma_true <- 0.15
vv <- profile_even(J, 0.04, 0.5)
jensen <- replicate(300, {
  tau_j <- rnorm(J, 0, sigma_true)
  y <- rnorm(J, tau_j, sqrt(vv))
  fit <- fit_rubin_grid(y, vv, M = 1200)
  gm <- exp(mean(log(vv)))
  I_plug <- mean(fit$sigma^2) / (mean(fit$sigma^2) + gm) # I(E[sigma^2 | y])
  I_post <- mean(fit$sigma^2 / (fit$sigma^2 + gm)) # E[I(sigma^2) | y]
  c(plug = I_plug, post = I_post)
})
check(
  "I(E[sigma^2|y]) >= E[I(sigma^2)|y] in every replication",
  all(jensen["plug", ] >= jensen["post", ] - 1e-10),
  sprintf("mean plug %.4f vs mean post %.4f", mean(jensen["plug", ]), mean(jensen["post", ]))
)
keynum$jensen_gap <- mean(jensen["plug", ] - jensen["post", ])

note("\n## 6. J_eff bounds and behavior")

check(
  "1 <= J_eff <= J with equality iff constant",
  abs(J_eff(rep(0.5, 30)) - 30) < 1e-12 &&
    J_eff(c(0.9, rep(1e-6, 29))) < 1.01
)

note("\n## 7. Fitting machinery sanity (8 schools; oracle recovery)")

# Rubin (1981) 8-schools data
y8 <- c(28, 8, -3, 7, -1, 1, 18, 12)
v8 <- c(15, 10, 16, 11, 9, 11, 10, 18)^2
fit8 <- fit_rubin_grid(y8, v8, M = 8000, seed = 1)
# published reference values: posterior means of tau_j all pulled toward ~8,
# school A largest with PM around 10-12 (BDA3, table 5.3)
check(
  "8-schools: PM ordering preserved and school A in [9, 13]",
  which.max(fit8$pm) == 1 && fit8$pm[1] > 9 && fit8$pm[1] < 13,
  sprintf("PM_A = %.2f, median sigma = %.2f", fit8$pm[1], median(fit8$sigma))
)
check(
  "8-schools: posterior median of sigma in [4, 10] (BDA ~ 6-8)",
  median(fit8$sigma) > 4 && median(fit8$sigma) < 10
)

# oracle recovery at large J, small v: sigma posterior concentrates near truth
J <- 200
sigma_true <- 0.2
tau_j <- rnorm(J, 0.1, sigma_true)
y <- rnorm(J, tau_j, sqrt(0.002))
fitL <- fit_rubin_grid(y, rep(0.002, J), M = 3000)
check(
  "large-J recovery: posterior mean of sigma within 10% of truth",
  abs(mean(fitL$sigma) - sigma_true) / sigma_true < 0.10,
  sprintf("E[sigma|y] = %.4f (truth %.2f)", mean(fitL$sigma), sigma_true)
)

# CB restores finite-population spread: var(CB) matches sigma2_hat
fitC <- fit_rubin_grid(y8, v8, M = 8000, seed = 2)
cb <- est_cb(fitC)
vfin <- var(fitC$pm)
sig2h <- mean(fitC$pv) + vfin
check(
  "CB variance equals sigma2_hat (Ghosh 1992 calibration)",
  abs(var(cb) - sig2h) / sig2h < 1e-10,
  sprintf("var(CB) = %.2f vs sigma2_hat = %.2f (var(PM) = %.2f)", var(cb), sig2h, vfin)
)

# GR: estimates are a permutation of the Gbar quantiles; ranks follow Rbar
gr <- est_gr(fitC)
check(
  "GR output is J mid-quantiles of pooled posterior in Rbar order",
  length(unique(gr)) == length(gr) &&
    all(sort(gr) == quantile(fitC$theta, probs = (2 * 1:8 - 1) / 16, type = 1, names = FALSE))
)

note("\n## 8. ISEL implementation vs brute-force numerical integral")

a <- rnorm(30)
b <- rnorm(30, 0.3)
zz <- seq(min(c(a, b)) - 1, max(c(a, b)) + 1, length.out = 2e5)
num <- sum((ecdf(a)(zz) - ecdf(b)(zz))^2) * diff(zz)[1]
check(
  "exact ISEL matches numeric integral",
  abs(loss_isel(a, b) - num) < 1e-3,
  sprintf("exact %.6f vs numeric %.6f", loss_isel(a, b), num)
)

note("\n## 9. Exact identity dH/dm = S_bar - I (fixed centered profile)")

set.seed(99)
for (Itar in c(0.2, 0.5, 0.8)) {
  m0 <- qlogis(Itar)
  D <- rnorm(40)
  D <- (D - mean(D)) / sqrt(mean((D - mean(D))^2)) * sqrt(1.3)
  Hm <- function(mm) mean(log1p(exp(mm + D))) - log1p(exp(mm))
  dH <- (Hm(m0 + 1e-6) - Hm(m0 - 1e-6)) / 2e-6
  gap <- mean(plogis(m0 + D)) - Itar
  check(
    sprintf("dH/dm equals S_bar - I at I = %.1f", Itar),
    abs(dH - gap) < 1e-8, sprintf("dH/dm %.7f vs gap %.7f", dH, gap)
  )
}

note("\n## 10. Sign of S_bar - I is profile-dependent at fixed (I, V)")

# Asymmetric counterexample (Appendix C, remark (b)): J = 25,
# I = .30, V = 2 exactly, one large positive logit deviation.
D <- c(sqrt(48), rep(-sqrt(1 / 12), 24))
S <- plogis(qlogis(0.30) + D)
check(
  "asymmetric profile at I = .30, V = 2 has S_bar < I",
  abs(mean(D)) < 1e-12 && abs(mean(D^2) - 2) < 1e-12 && mean(S) < 0.30,
  sprintf("S_bar = %.6f", mean(S))
)
# The concentrated-noisy mirror at the same (I, V) has S_bar > I.
vv <- profile_concentrated_noisy(25, 0.04 * (1 - 0.3) / 0.3, 2)
SS <- site_reliability(0.04, vv)
check(
  "concentrated-noisy profile at I = .30, V = 2 has S_bar > I",
  mean(SS) > 0.30, sprintf("S_bar = %.6f", mean(SS))
)

note("\n## 11. Mixed-endpoint limits need not reach a boundary")

for (eps in c(1e-3, 1e-6, 1e-9)) {
  Ie <- kl_centroid(c(eps, 1 - eps))
  # tolerance reflects double-precision representation error of 1 - eps
  # (the exact value is 1/2 for every eps; fp error is O(eps) near eps = 1e-9)
  check(
    sprintf("S = (%.0e, 1 - %.0e) gives I = 1/2 (logit cancellation)", eps, eps),
    abs(Ie - 0.5) < 1e-6, sprintf("I = %.10f", Ie)
  )
}

note("\n## 12. Rao-Blackwellized PM: exact homogeneous rank tie")

set.seed(7)
J <- 30
vhom <- rep(0.03, J)
tau_j <- rnorm(J, 0, 0.2)
y <- rnorm(J, tau_j, sqrt(vhom))
fitH <- fit_rubin_grid(y, vhom, M = 1000)
check(
  "rank(pm) == rank(y) under equal v_j (affine invariance)",
  all(rank(fitH$pm) == rank(y))
)
check(
  "rank(CB) == rank(y) under equal v_j",
  all(rank(est_cb(fitH)) == rank(y))
)

note("\n## 13. Proposition 1: Monte Carlo variance-ratio check")

set.seed(11)
sig2 <- 0.05
vj <- 0.02
n <- 4e5
tj <- rnorm(n, 0, sqrt(sig2))
xj <- rnorm(n, tj, sqrt(vj))
rel_mc <- var(tj) / var(xj)
check(
  "Var(tau)/Var(tau_hat) matches S_j = sigma2/(sigma2 + v_j)",
  abs(rel_mc - sig2 / (sig2 + vj)) < 0.005,
  sprintf("MC %.4f vs S %.4f", rel_mc, sig2 / (sig2 + vj))
)

dir.create("outputs", showWarnings = FALSE)
writeLines(report, "outputs/verification_report.md")
write.csv(data.frame(key = names(keynum), value = unlist(keynum)),
  "outputs/key_numbers.csv",
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_verification.txt")
note("\nReport written to outputs/verification_report.md")
if (ANY_FAIL) {
  note("\nAT LEAST ONE CHECK FAILED")
  quit(status = 1)
}
