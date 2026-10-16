# Exploratory pretest with Gaussian effects and J = 50.
# This checks whether dispersion and shape help explain performance at fixed I.
# It supplies no manuscript results; the main arms add effect laws and J values.
#
# The 5% criterion uses the raw increase in R-squared. Table 4,
# Panel A instead reports the fraction of remaining variation explained.
# Those quantities should not be compared directly.
# Run through analysis/scripts/run_historical.R --arm=pretest.

source("00_utils.R")
suppressPackageStartupMessages({
  library(parallel)
})

SIGMA_T <- 0.20 # cross-site SD, effect-size units
J <- 50
REPS <- 300
M_DRAWS <- 1500
SEED0 <- 20260710

grid <- expand.grid(
  I = c(0.15, 0.30, 0.45, 0.60),
  V = c(0, 0.5, 1.0, 2.0),
  shape = c("even", "bimodal", "concentrated"),
  stringsAsFactors = FALSE
)
grid <- grid[!(grid$V == 0 & grid$shape != "even"), ] # V = 0 collapses shapes
rownames(grid) <- NULL

run_cell <- function(g) {
  gm <- SIGMA_T^2 * (1 - g$I) / g$I
  v <- make_profile(g$shape, J, gm, g$V)
  dg <- diagnostics_profile(SIGMA_T^2, v)
  set.seed(SEED0 + round(1e5 * g$I) + round(17 * g$V) + nchar(g$shape))

  out <- vector("list", REPS)
  for (r in seq_len(REPS)) {
    tau_j <- rnorm(J, 0, SIGMA_T)
    y <- rnorm(J, tau_j, sqrt(v))
    fit <- fit_rubin_grid(y, v, M = M_DRAWS)
    ests <- list(ml = y, pm = est_pm(fit), cb = est_cb(fit), gr = est_gr(fit))
    out[[r]] <- do.call(rbind, lapply(names(ests), function(nm) {
      L <- all_losses(ests[[nm]], tau_j)
      data.frame(
        est = nm, rmse = L["rmse"], mselp = L["mselp"], isel = L["isel"],
        row.names = NULL
      )
    }))
  }
  res <- do.call(rbind, out)
  agg <- aggregate(cbind(rmse, mselp, isel) ~ est, data = res, FUN = mean)
  cbind(
    I = g$I, V = g$V, shape = g$shape,
    H = dg$H, J_eff = dg$J_eff, S_bar = dg$S_bar, agg
  )
}

cells <- mclapply(seq_len(nrow(grid)), function(i) run_cell(grid[i, ]),
  mc.cores = max(1L, detectCores() - 2L)
)
res <- do.call(rbind, cells)
dir.create("outputs", showWarnings = FALSE)
write.csv(res, "outputs/pretest_cellmeans.csv", row.names = FALSE)

# Increase in R-squared from adding dispersion and shape, by loss and estimator.

cat("\n== Incremental R^2 of (H, shape) beyond factor(I) on log mean loss ==\n")
gate <- expand.grid(
  goal = c("rmse", "mselp", "isel"),
  est = c("pm", "cb", "gr"), stringsAsFactors = FALSE
)
gate$R2_I <- gate$R2_IH <- gate$inc <- NA_real_
for (k in seq_len(nrow(gate))) {
  d <- res[res$est == gate$est[k], ]
  yy <- log(d[[gate$goal[k]]])
  m0 <- lm(yy ~ factor(I), data = d)
  m1 <- lm(yy ~ factor(I) + H + factor(shape), data = d)
  gate$R2_I[k] <- summary(m0)$r.squared
  gate$R2_IH[k] <- summary(m1)$r.squared
  gate$inc[k] <- gate$R2_IH[k] - gate$R2_I[k]
}
print(transform(gate,
  R2_I = round(R2_I, 3), R2_IH = round(R2_IH, 3),
  inc = round(inc, 3)
), row.names = FALSE)
write.csv(gate, "outputs/pretest_gate.csv", row.names = FALSE)

cat("\n== Does the BEST method change with dispersion at fixed I? ==\n")
for (goal in c("rmse", "mselp", "isel")) {
  b <- do.call(rbind, lapply(
    split(res, list(res$I, res$V, res$shape), drop = TRUE),
    function(d) data.frame(
      I = d$I[1], V = d$V[1], shape = d$shape[1],
      best = d$est[which.min(d[[goal]])]
    )
  ))
  tab <- with(b, table(I, best))
  cat("\n--", goal, "-- best-method counts by I (across V x shape cells):\n")
  print(tab)
}

cat("\n== Direction check at I = .30 (theory: higher V lowers PM rmse) ==\n")
d <- subset(res, est == "pm" & I == 0.30 & shape %in% c("even"))
print(d[order(d$V), c("V", "H", "rmse", "mselp", "isel")], row.names = FALSE)
cat(
  "\nGATE:", if (max(gate$inc, na.rm = TRUE) >= 0.05) "PASS" else "FAIL",
  sprintf("(max incremental R^2 = %.3f)\n", max(gate$inc, na.rm = TRUE))
)
