# Monte Carlo SEs for the fixed-grid geometric loss ratios (OSM F.3).
# Rows sharing a fitted-posterior seed are grouped before estimating variance.
# The calculation assumes iid latent rows and independent distinct seed streams.
# Run last in the saved-loss analysis; run_all.sh --simulation sets the order.

library(jsonlite)
r <- readRDS('inputs/sim_main_reps.rds')
b <- subset(r, V == 0, c(family, rep, est, rmse, mselp, isel))
names(b)[4:6] <- paste0(names(b)[4:6], '0')
# Join paired rows that share a numerator or balanced-profile fit seed.
components <- function(family, rep, combo) {
  n <- length(family)
  parent <- seq_len(n)
  findroot <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  seeds <- c(
    20260710L + 1000L * family + 61L * combo + 100000L + rep,
    20260710L + 1000L * family + 61L + 100000L + rep
  )
  rows <- rep(seq_len(n), 2)
  first <- match(seeds, seeds)
  edge <- which(first != seq_along(seeds))
  for (j in edge) {
    a <- findroot(rows[j])
    b <- findroot(rows[first[j]])
    if (a != b) parent[b] <- a
  }
  comp <- vapply(seq_len(n), findroot, integer(1))
  parts <- split(seq_len(n), comp)
  types <- vapply(parts, function(ii) paste(sort(family[ii]), collapse = ','), character(1))
  list(
    comp = comp,
    parts = parts,
    types = types,
    seed_edges = length(edge),
    ncomponents = length(parts),
    maxsize = max(lengths(parts))
  )
}
out <- list()
details <- list()
kk <- 0
for (shp in c('even', 'bimodal', 'concentrated', 'concentrated_noisy')) {
  m <- merge(subset(r, V == 2 & shape == shp), b, by = c('family', 'rep', 'est'))
  for (z in split(m, list(m$est, m$I), drop = TRUE)) {
    # The reported target averages these nine families, not a sampled population.
    K <- length(unique(z$family))
    R <- unique(as.integer(table(z$family)))
    stopifnot(length(R) == 1, R == 500, K == 9, !anyDuplicated(z[c('family', 'rep')]))
    ell <- cbind(
      rmse = log(z$rmse / z$rmse0),
      mselp = log(z$mselp / z$mselp0),
      isel = log(z$isel / z$isel0)
    )
    cc <- components(z$family, z$rep, z$combo)
    sums <- rowsum(ell, cc$comp)
    types <- cc$types[rownames(sums)]
    B <- table(types)
    stopifnot(min(B) >= 2, sum(lengths(cc$parts)) == nrow(z))
    # Estimate variance separately for groups with the same family membership.
    vg <- Reduce(
      `+`,
      lapply(
        split(seq_len(nrow(sums)), types),
        function(ii) length(ii) * apply(sums[ii, , drop = FALSE], 2, var)
      )
    ) / (K * R)^2
    fam <- rowsum(ell, z$family) / R
    point <- exp(colMeans(fam))
    naive <- Reduce(
      `+`,
      lapply(split(as.data.frame(ell), z$family), function(q) apply(q, 2, var) / nrow(q))
    ) / K^2
    # Transform log-ratio uncertainty to the ratio scale by the delta method.
    se <- point * sqrt(vg)
    kk <- kk + 1
    out[[kk]] <- data.frame(
      est = z$est[1],
      I = z$I[1],
      shape = shp,
      r_rmse = point[1],
      se_rmse = se[1],
      r_mselp = point[2],
      se_mselp = se[2],
      r_isel = point[3],
      se_isel = se[3],
      log_mcse_rmse = sqrt(vg[1]),
      log_mcse_mselp = sqrt(vg[2]),
      log_mcse_isel = sqrt(vg[3]),
      naive_se_rmse = point[1] * sqrt(naive[1]),
      naive_se_mselp = point[2] * sqrt(naive[2]),
      naive_se_isel = point[3] * sqrt(naive[3]),
      logratio_family_sd_rmse = sd(fam[, 1]),
      logratio_family_sd_mselp = sd(fam[, 2]),
      logratio_family_sd_isel = sd(fam[, 3]),
      uncertainty_status = 'component-type delta MCSE; recorded seed sharing; independent unique-seed-stream assumption',
      row.names = NULL
    )
    details[[kk]] <- data.frame(
      est = z$est[1],
      I = z$I[1],
      shape = shp,
      K,
      R,
      n = nrow(z),
      components = cc$ncomponents,
      types = length(B),
      min_components_per_type = min(B),
      max_component_size = cc$maxsize,
      shared_seed_edges = cc$seed_edges
    )
  }
}
out <- do.call(rbind, out)
out <- out[order(out$shape, out$est, out$I), ]
detail <- do.call(rbind, details)
old <- read.csv('inputs/tables/tbl_ratios_shape.csv')
mm <- merge(out, old, by = c('est', 'I', 'shape'), suffixes = c('', '.old'))
diff <- sapply(
  c('r_rmse', 'r_mselp', 'r_isel'),
  function(g) max(abs(mm[[g]] - mm[[paste0(g, '.old')]]))
)
stopifnot(max(diff) < 1e-12)
write.csv(out, 'outputs/tables/tbl_ratios_shape.csv', row.names = FALSE)
ev <- subset(out, shape == 'even')
write.csv(ev, 'outputs/tables/tbl_ratios.csv', row.names = FALSE)
write.csv(
  subset(ev, est %in% c('opm', 'spm', 'pm'), select = c(est, I, r_rmse, se_rmse)),
  'outputs/tables/tbl_oracle.csv',
  row.names = FALSE
)
write.csv(detail, 'outputs/evidence-audit/component-mcse-structure.csv', row.names = FALSE)
write.csv(mm, 'outputs/evidence-audit/paired-mcse-comparison.csv', row.names = FALSE)
summary <- list(
  point_maxabs = diff,
  max_component_size = max(detail$max_component_size),
  min_components_per_type = min(detail$min_components_per_type),
  family_se_to_component_mcse_range = range(unlist(mm[c('se_rmse.old', 'se_mselp.old', 'se_isel.old')]) / unlist(mm[c('se_rmse', 'se_mselp', 'se_isel')])),
  all_shapes_max_se = sapply(out[c('se_rmse', 'se_mselp', 'se_isel')], max),
  scope = 'Fixed K=9 family geometric point target. Components join rows sharing either actual numerator or baseline fit seed; type is sorted exact family-membership multiset. Sum B_t Var(component sums)/(KR)^2, delta transform. Depends on intended independent unique RNG streams; no claim of formal independence for actual PRNG bits. Raw estimators use same valid coarsened partition.'
)
write_json(
  summary,
  'outputs/evidence-audit/component-mcse-summary.json',
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = 16
)
print(summary)
