#!/usr/bin/env Rscript
# Check OCRS site estimates and run the empirical sensitivity analyses.
# With --aggregate-reference, use saved site summaries without participant data.
# Other modes read the participant CSV supplied through OCRS_DATA.
# Outputs contain site and block summaries, not participant records or identifiers.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]))
analysis_dir <- dirname(dirname(script_path))
in_dir <- file.path(analysis_dir, "inputs")
aggregate_only <- "--aggregate-reference" %in% args
if (aggregate_only && any(c("--rerun-sensitivity", "--pilot", "--final-checks") %in% args))
  stop("--aggregate-reference cannot be combined with raw-data audit options")
saved <- readRDS(file.path(in_dir, "ocrs_fits.rds"))
raw_sha256 <- "d20cbbaaf07ee8e7aad7e5d0e516b74f0a1114617fb3c7db9a4c7c2afcb681c7"
raw_path <- Sys.getenv(
  "OCRS_DATA",
  file.path(dirname(analysis_dir), "data", "OCRS_publicuse_impactanalysisdata.csv")
)
out_dir <- file.path(analysis_dir, "outputs", "ocrs-audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (aggregate_only) {
  # Start from the supplied summaries. The participant-data mode recreates them.
  aggregate_files <- c(
    "site-audit.csv", "bootstrap-replication.rds",
    "subset-replication.rds", "bootstrap-summary.csv"
  )
  for (name in aggregate_files) {
    baseline <- file.path(in_dir, "ocrs-aggregate", name)
    if (!file.exists(baseline)) stop("Missing aggregate input: ", name)
    if (!file.copy(baseline, file.path(out_dir, name), overwrite = TRUE))
      stop("Cannot copy aggregate input: ", name)
  }
}

source(file.path(analysis_dir, "scripts", "legacy", "00_utils.R"))

# Load the four helper functions without running the grid analysis.
legacy <- parse(file.path(analysis_dir, "scripts", "legacy", "04_ocrs.R"))
for (ex in legacy) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) &&
    is.symbol(ex[[2]]) && as.character(ex[[2]]) %in%
    c("hc2_ols", "site_effect_block", "analyze", "summarize_fit")) eval(ex)
}

# Check how flooring 1 - leverage at 1e-8 affects the HC2 treatment variance.
# Use the same independent design columns as hc2_ols().
hc2_clamp_diagnostics <- function(X, y) {
  qx <- qr(X)
  keep <- qx$pivot[seq_len(qx$rank)]
  Xk <- X[, keep, drop = FALSE]
  inv <- chol2inv(chol(crossprod(Xk)))
  e <- as.numeric(y - Xk %*% (inv %*% crossprod(Xk, y)))
  h <- rowSums((Xk %*% inv) * Xk)
  clamp <- 1 - h < 1e-8
  g <- as.numeric(inv[match("ragrp", colnames(Xk)), , drop = FALSE] %*% t(Xk))
  contrib <- g^2 * e^2 / pmax(1 - h, 1e-8)
  c(
    clamp_rows = sum(clamp), max_h = max(h),
    clamped_treatment_score_max = if (any(clamp)) max(abs(g[clamp])) else 0,
    clamped_variance_contribution = sum(contrib[clamp]),
    clamped_variance_fraction = sum(contrib[clamp]) / sum(contrib)
  )
}

# Fit the block-adjusted effect and report any effect of the leverage floor.
site_effect_clamp <- function(dd, outcome) {
  X <- cbind(`(Intercept)` = 1, ragrp = dd$ragrp)
  if (length(unique(dd$rablock)) > 1L)
    X <- cbind(X, model.matrix(~ factor(rablock), dd)[, -1, drop = FALSE])
  c(hc2_ols(X, dd[[outcome]]), hc2_clamp_diagnostics(X, dd[[outcome]]))
}

write_csv <- function(x, name) write.csv(x, file.path(out_dir, name), row.names = FALSE)
write_json <- function(x, name) jsonlite::write_json(x, file.path(out_dir, name),
  pretty = TRUE, auto_unbox = TRUE, digits = 16, null = "null", na = "null"
)
qsum <- function(x) c(min = min(x), median = median(x), max = max(x))
if (!aggregate_only) {
  if (!file.exists(raw_path)) stop("OCRS participant file not supplied. Set OCRS_DATA to the external CSV; see data/DATA_ACCESS.md. Use --aggregate-reference for saved aggregate sensitivities.")
  actual_raw_sha256 <- digest::digest(raw_path, file = TRUE, algo = "sha256")
  if (!identical(actual_raw_sha256, raw_sha256))
    stop("OCRS_DATA SHA256 does not match the manuscript input; refusing a different data release. See data/DATA_ACCESS.md.")
  d <- read.csv(raw_path, stringsAsFactors = FALSE)
  d$subject <- ifelse(d$crssub == "A1", "Algebra 1", "English 9")
  d$site <- paste(d$studyschids, d$subject, sep = " x ")
  d$credits_y4 <- ifelse(d$subject == "Algebra 1", d$cby_y4_math, d$cby_y4_eng)
  required <- c(
    "stdids", "studyschids", "rablock", "ragrp", "term", "crssub",
    "scm_passcrs", "credits_y4", "state_math_imp", "state_ela_imp", "state_miss",
    "gpa9th_imp", "gpa9th_miss"
  )
  missingness <- data.frame(
    field = required,
    missing = vapply(d[required], function(x) sum(is.na(x)), integer(1))
  )
  write_csv(missingness, "required-input-missingness.csv")
  stopifnot(!any(missingness$missing), all(d$ragrp %in% 0:1), all(d$scm_passcrs %in% 0:1))
  site_groups <- split(d, d$site)

  # Original order retained internally only; exported identifiers are S01,... .
  site_order <- order(
    vapply(site_groups, function(x) x$subject[1], character(1)),
    vapply(site_groups, function(x) x$studyschids[1], numeric(1))
  )
  site_groups <- site_groups[site_order]
  site_names <- names(site_groups)
  site_codes <- setNames(sprintf("S%02d", seq_along(site_groups)), site_names)

  # Independently check the treatment coefficient and HC2 variance by
  # residualizing treatment on the other regressors (Frisch-Waugh-Lovell).
  matrix_audit <- function(dd, outcome, adjust = FALSE) {
    X <- cbind(`(Intercept)` = 1, ragrp = dd$ragrp)
    if (length(unique(dd$rablock)) > 1L)
      X <- cbind(X, model.matrix(~ factor(rablock), dd)[, -1, drop = FALSE])
    if (adjust) {
      pre <- if (dd$subject[1] == "Algebra 1") dd$state_math_imp else dd$state_ela_imp
      X <- cbind(X, pre = pre, as.matrix(dd[, c("state_miss", "gpa9th_imp", "gpa9th_miss")]))
    }
    qx <- qr(X)
    keep <- qx$pivot[seq_len(qx$rank)]
    Xk <- X[, keep, drop = FALSE]
    fit <- lm.fit(Xk, dd[[outcome]])
    Q <- qr.Q(qr(Xk))
    h <- rowSums(Q^2)
    C <- Xk[, colnames(Xk) != "ragrp", drop = FALSE]
    rz <- qr.resid(qr(C), dd$ragrp)
    denom <- sum(rz^2)
    b_fwl <- sum(rz * dd[[outcome]]) / denom
    v_fwl <- sum(rz^2 * fit$residuals^2 / pmax(1 - h, 1e-8)) / denom^2
    legacy_value <- hc2_ols(X, dd[[outcome]])
    c(
      n = nrow(X), columns = ncol(X), rank = qx$rank, alias_columns = ncol(X) - qx$rank,
      residual_df = nrow(X) - qx$rank, max_leverage = max(h),
      leverage_floor_rows = sum(1 - h < 1e-8), treatment_residual_ss = denom,
      est = unname(legacy_value["est"]), v = unname(legacy_value["se"]^2),
      independent_coef_absdiff = abs(b_fwl - legacy_value["est"]),
      independent_v_absdiff = abs(v_fwl - legacy_value["se"]^2),
      hc2_clamp_diagnostics(X, dd[[outcome]])
    )
  }
  site_rows <- list()
  matrix_rows <- list()
  block_rows <- list()
  for (j in seq_along(site_groups)) {
    dd <- site_groups[[j]]
    code <- site_codes[names(site_groups)[j]]
    blocks <- split(dd, dd$rablock)
    bl <- do.call(rbind, lapply(blocks, function(b) {
      n <- nrow(b)
      n1 <- sum(b$ragrp == 1)
      n0 <- n - n1
      p <- n1 / n
      contrast <- if (n1 > 0 && n0 > 0) mean(b$scm_passcrs[b$ragrp == 1]) - mean(b$scm_passcrs[b$ragrp == 0]) else NA_real_
      c(n = n, n0 = n0, n1 = n1, p = p, overlap = n * p * (1 - p), rd = contrast)
    }))
    informative <- bl[, "overlap"] > 0
    wfe <- bl[, "overlap"] / sum(bl[, "overlap"])
    wperson <- bl[, "n"] / sum(bl[, "n"])
    wb <- ifelse(informative, wfe * ifelse(is.na(bl[, "rd"]), 0, bl[, "rd"]), 0)
    un <- matrix_audit(dd, "scm_passcrs")
    ad <- matrix_audit(dd, "scm_passcrs", TRUE)
    cr <- matrix_audit(dd, "credits_y4")
    for (spec in c("block_rd", "adjusted_rd", "block_credits")) {
      vals <- switch(spec,
        block_rd = un,
        adjusted_rd = ad,
        block_credits = cr
      )
      matrix_rows[[length(matrix_rows) + 1L]] <- data.frame(
        site = unname(code),
        subject = dd$subject[1], specification = spec, t(vals), check.names = FALSE
      )
    }
    n1 <- sum(dd$ragrp == 1)
    n0 <- sum(dd$ragrp == 0)
    p1 <- mean(dd$scm_passcrs[dd$ragrp == 1])
    p0 <- mean(dd$scm_passcrs[dd$ragrp == 0])
    var_unadjusted <- var(dd$scm_passcrs[dd$ragrp == 1]) / n1 +
      var(dd$scm_passcrs[dd$ragrp == 0]) / n0
    ap <- sum(dd$scm_passcrs[dd$ragrp == 1])
    bf <- n1 - ap
    cp <- sum(dd$scm_passcrs[dd$ragrp == 0])
    df <- n0 - cp
    cc <- if (min(ap, bf, cp, df) == 0) .5 else 0
    logor <- log((ap + cc) / (bf + cc)) - log((cp + cc) / (df + cc))
    vlogor <- sum(1 / (c(ap, bf, cp, df) + cc))
    floor_rd <- sum((.5 / (c(n0, n1) + 1)) * (1 - .5 / (c(n0, n1) + 1)) / c(n0, n1))

    # Size-weighted effects cover only blocks with both treatment arms.
    # With unsupported blocks, they do not estimate the complete site's SATE.
    winf <- bl[informative, "n"] / sum(bl[informative, "n"])
    rd_informative_sizeweighted <- sum(winf * bl[informative, "rd"])

    # Each arm needs at least two students to estimate its sample variance.
    all_big <- all(bl[informative, "n0"] >= 2 & bl[informative, "n1"] >= 2)
    v_informative_sizeweighted <- if (all_big) {
      vv <- vapply(blocks[informative], function(b)
        var(b$scm_passcrs[b$ragrp == 1]) / sum(b$ragrp == 1) +
          var(b$scm_passcrs[b$ragrp == 0]) / sum(b$ragrp == 0), numeric(1))
      sum(winf^2 * vv)
    } else NA_real_
    termp <- vapply(split(dd, dd$term), function(x) mean(x$ragrp), numeric(1))
    site_rows[[j]] <- data.frame(
      site = unname(code), subject = dd$subject[1], n = nrow(dd),
      n0 = n0, n1 = n1, blocks = nrow(bl), informative_blocks = sum(informative),
      unsupported_blocks = sum(!informative), unsupported_records = sum(bl[!informative, "n"]),
      blocks_small_arm = sum(informative & (bl[, "n0"] < 2 | bl[, "n1"] < 2)),
      n_terms = length(termp), term_allocation_range = diff(range(termp)),
      weights_l1 = sum(abs(wfe - wperson)), weights_maxabs = max(abs(wfe - wperson)),
      overlap_product_range = diff(range(bl[, "p"] * (1 - bl[, "p"]))),
      rd = un["est"], v = un["v"], rd_adj = ad["est"], v_adj = ad["v"],
      credits = cr["est"], v_credits = cr["v"], rd_unadjusted = p1 - p0,
      v_unadjusted = var_unadjusted, log_or = logor, v_log_or = vlogor,
      correction = cc, binary_floor = floor_rd, binary_floor_fraction = floor_rd / un["v"],
      credits_floor = .05 * cr["v"], contrast_identity_error = abs(sum(wb) - un["est"]),
      rd_informative_personweighted = rd_informative_sizeweighted,
      v_informative_personweighted = v_informative_sizeweighted, check.names = FALSE
    )
    block_rows[[j]] <- data.frame(
      n = bl[, "n"], n0 = bl[, "n0"], n1 = bl[, "n1"],
      informative = informative
    )
  }
  sites <- do.call(rbind, site_rows)
  rownames(sites) <- NULL
  mats <- do.call(rbind, matrix_rows)
  rownames(mats) <- NULL
  blall <- do.call(rbind, block_rows)
  write_csv(sites, "site-audit.csv")
  write_csv(mats, "design-matrix-audit.csv")
  block_dist <- aggregate(rep(1, nrow(blall)), by = blall, FUN = sum)
  names(block_dist)[ncol(block_dist)] <- "block_count"
  write_csv(block_dist, "block-support-distribution.csv")
  fits_inputs <- list(
    pooled = list(y = sites$rd, v = sites$v),
    algebra = list(y = sites$rd[sites$subject == "Algebra 1"], v = sites$v[sites$subject == "Algebra 1"]),
    english = list(y = sites$rd[sites$subject == "English 9"], v = sites$v[sites$subject == "English 9"]),
    pooled_adj = list(y = sites$rd_adj, v = sites$v_adj),
    credits = list(y = sites$credits, v = sites$v_credits),
    logodds = list(y = sites$log_or, v = sites$v_log_or)
  )

  # Saved log-odds rows use lexical order rather than subject/school order.
  ord_lex <- order(site_names)
  fits_inputs$logodds <- list(y = sites$log_or[ord_lex], v = sites$v_log_or[ord_lex])
  input_match <- do.call(rbind, lapply(names(fits_inputs), function(nm) {
    old <- if (nm == "logodds") saved$logodds_fit else saved$fits[[nm]]
    new <- fits_inputs[[nm]]
    data.frame(
      fit = nm, J = length(new$y), max_y_diff = max(abs(new$y - old$y)),
      max_v_diff = max(abs(new$v - old$v)), all_v_positive = all(new$v > 0)
    )
  }))
  write_csv(input_match, "saved-input-comparison.csv")
  stopifnot(
    max(input_match$max_y_diff) < 1e-10, max(input_match$max_v_diff) < 1e-10,
    max(sites$contrast_identity_error) < 1e-10,
    max(mats$independent_coef_absdiff) < 1e-9,
    max(mats$independent_v_absdiff) < 1e-9
  )

  # Count repeated student records without writing student identifiers.
  student_groups <- split(d, d$stdids)
  student_n <- vapply(student_groups, nrow, integer(1))
  school_sites <- vapply(split(d, d$studyschids), function(x) length(unique(x$site)), integer(1))
  shared <- student_groups[student_n > 1]
  overview <- list(
    status = "audited_aggregate_inputs", rows = nrow(d), columns_raw = ncol(d) - 3L,
    unique_students = length(student_groups), students_multiple_records = sum(student_n > 1),
    students_multiple_subjects = sum(vapply(
      student_groups,
      function(x) length(unique(x$subject)) > 1,
      logical(1)
    )),
    students_multiple_schools = sum(vapply(
      student_groups,
      function(x) length(unique(x$studyschids)) > 1,
      logical(1)
    )),
    shared_students_same_arm = sum(vapply(shared, function(x) length(unique(x$ragrp)) == 1, logical(1))),
    duplicate_student_subject_keys = sum(duplicated(paste(d$stdids, d$subject))),
    exact_duplicate_rows = sum(duplicated(d)), schools = length(unique(d$studyschids)),
    sites = nrow(sites), algebra_sites = sum(sites$subject == "Algebra 1"),
    english_sites = sum(sites$subject == "English 9"), paired_schools = sum(school_sites == 2),
    total_blocks = nrow(blall), informative_blocks = sum(blall$informative),
    unsupported_blocks = sum(!blall$informative), unsupported_records = sum(blall$n[!blall$informative]),
    sites_with_unsupported_blocks = sum(sites$unsupported_blocks > 0),
    informative_blocks_arm_below2 = sum(blall$informative & (blall$n0 < 2 | blall$n1 < 2)),
    sites_person_weights_differ = sum(sites$weights_maxabs > 1e-12),
    site_n = qsum(sites$n), site_arm_n = qsum(c(sites$n0, sites$n1)),
    blocks_per_site = qsum(sites$blocks), sites_multiple_blocks = sum(sites$blocks > 1),
    sites_multiple_terms = sum(sites$n_terms > 1),
    term_allocation_range_multiterm = qsum(sites$term_allocation_range[sites$n_terms > 1]),
    weights_l1 = qsum(sites$weights_l1), weights_maxabs = qsum(sites$weights_maxabs),
    absolute_block_adjustment = qsum(abs(sites$rd - sites$rd_unadjusted)),
    v_positive = all(c(sites$v, sites$v_adj, sites$v_credits, sites$v_log_or, sites$v_unadjusted) > 0),
    logodds_zero_cells = sum(sites$correction > 0),
    credits_range = c(
      Algebra_min = min(d$cby_y4_math[d$subject == "Algebra 1"]),
      Algebra_max = max(d$cby_y4_math[d$subject == "Algebra 1"]),
      English_min = min(d$cby_y4_eng[d$subject == "English 9"]),
      English_max = max(d$cby_y4_eng[d$subject == "English 9"])
    ),
    test_missing_fraction = c(
      sts_pcorrect = mean(is.na(d$sts_pcorrect)),
      sts_irtscore = mean(is.na(d$sts_irtscore))
    ),
    manuscript_values_status = "Counts and first-stage y/v independently reproduced; posterior numerical convergence audited separately.",
    source_mode = "external participant CSV verified by SHA256; no raw participant rows output"
  )
  write_json(overview, "overview.json")
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  cat(sprintf(
    "Input checks: %d records, %d students, %d sites, %d blocks; %d unsupported blocks (%d records); all six saved y/v vectors matched.\n",
    nrow(d),
    length(student_groups),
    nrow(sites),
    nrow(blall),
    sum(!blall$informative),
    sum(blall$n[!blall$informative])
  ))
  cat(sprintf(
    "Person-weight mismatch: %d/%d sites; HC2 variance checks passed.\n",
    sum(sites$weights_maxabs > 1e-12),
    nrow(sites)
  ))
} else {
  sites <- read.csv(file.path(out_dir, "site-audit.csv"), check.names = FALSE)
  stopifnot(nrow(sites) == 32L, all(sites$v > 0), all(sites$v_credits > 0))
}

# Summarize existing posterior draws without drawing new samples.
# Check the faster dual-Jensen calculation of H against the direct KL formula.
summary_draws <- function(fit, v) {
  a <- fit$sigma^2
  gm <- exp(mean(log(v)))
  I <- a / (a + gm)
  S <- outer(a, v, function(aa, vv) aa / (aa + vv))
  eta <- outer(log(a), log(v), `-`)
  m <- log(a) - mean(log(v))
  A <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  H <- rowMeans(A(eta)) - A(m)
  c(
    I = mean(I), H = mean(H), V = mean((log(v) - mean(log(v)))^2),
    I_lo = unname(quantile(I, .025)), I_hi = unname(quantile(I, .975)),
    Pr_I_gt_half = mean(I > .5), Sbar = mean(S),
    Jeff = mean(rowSums(S)^2 / rowSums(S^2)),
    Ea = mean(a), Iplug = mean(a) / (mean(a) + gm)
  )
}

if (any(c("--rerun-sensitivity", "--pilot") %in% args)) {
  pilot <- "--pilot" %in% args
  run_B <- if (pilot) 5L else 500L
  prefix <- if (pilot) "pilot-" else ""
  t0 <- proc.time()[[3]]

  # Use the same site estimates and seeds for all one-site-per-school subsets.
  if (!pilot) {
    school_ids <- vapply(site_groups, function(x) x$studyschids[1], numeric(1))
    paired <- names(which(table(school_ids) == 2))
    grid8 <- expand.grid(rep(list(c("Algebra 1", "English 9")), length(paired)),
      stringsAsFactors = FALSE
    )
    subset_out <- vector("list", nrow(grid8))
    subset_keep <- matrix(FALSE, nrow(grid8), nrow(sites), dimnames = list(NULL, sites$site))
    for (k in seq_len(nrow(grid8))) {
      keep <- rep(TRUE, nrow(sites))
      for (pi in seq_along(paired)) {
        drop_subj <- setdiff(c("Algebra 1", "English 9"), grid8[k, pi])
        keep[as.character(school_ids) == paired[pi] & sites$subject == drop_subj] <- FALSE
      }
      fp <- fit_rubin_grid(sites$rd[keep], sites$v[keep], M = 4000L, seed = 1000 + k)
      fc <- fit_rubin_grid(sites$credits[keep], sites$v_credits[keep], M = 4000L, seed = 5000 + k)
      gp <- exp(mean(log(sites$v[keep])))
      gc <- exp(mean(log(sites$v_credits[keep])))
      subset_out[[k]] <- data.frame(
        subset = k, J = sum(keep),
        I_pass = mean(fp$sigma^2 / (fp$sigma^2 + gp)),
        I_credits = mean(fc$sigma^2 / (fc$sigma^2 + gc))
      )
      subset_keep[k, ] <- keep
    }
    subset_out <- do.call(rbind, subset_out)
    write_csv(subset_out, "subset-replication.csv")
    saveRDS(
      list(
        summary = subset_out, keep = subset_keep,
        method = "legacy250grid_M4000_seeds1000k_5000k"
      ),
      file.path(out_dir, "subset-replication.rds")
    )
    subset_summary <- data.frame(
      outcome = c("passing", "credits"),
      I_min = c(min(subset_out$I_pass), min(subset_out$I_credits)),
      I_med = c(median(subset_out$I_pass), median(subset_out$I_credits)),
      I_max = c(max(subset_out$I_pass), max(subset_out$I_credits))
    )
    write_csv(subset_summary, "subset-summary.csv")
    old <- saved$sens$unique_school
    subset_comparison <- data.frame(
      outcome = old$outcome,
      max_abs_diff = vapply(seq_len(nrow(old)), function(i)
        max(abs(as.numeric(old[i, -1]) - as.numeric(subset_summary[i, -1]))), numeric(1))
    )
    write_csv(subset_comparison, "subset-saved-comparison.csv")
    cat(sprintf(
      "Subset enumeration complete in %.1fs; maximum difference from saved results %.3g.\n",
      proc.time()[[3]] - t0, max(subset_comparison$max_abs_diff)
    ))
  }

  # Refit the grid-based log-odds model to reproduce the bootstrap's starting seed.
  flo_input <- fits_inputs$logodds
  flo_fit <- fit_rubin_grid(flo_input$y, flo_input$v, M = 20000L, seed = 7)
  boot_initial_state <- .Random.seed
  split_cells <- split(seq_len(nrow(d)), list(d$site, d$ragrp, d$rablock), drop = TRUE)

  # Resample students within site, arm, and randomization block, then refit.
  # These intervals describe first-stage resampling, not the conditional posterior.
  boot_run <- function(outcome, B) {
    # Keep resampled variances positive: use the continuity-corrected minimum
    # for course passing and 5% of the original variance for credits.
    v_floor <- vapply(split(d, d$site), function(dd) {
      if (outcome == "scm_passcrs") {
        n1 <- sum(dd$ragrp == 1)
        n0 <- sum(dd$ragrp == 0)
        pm1 <- .5 / (n1 + 1)
        pm0 <- .5 / (n0 + 1)
        pm1 * (1 - pm1) / n1 + pm0 * (1 - pm0) / n0
      } else .05 * site_effect_block(dd, outcome)["se"]^2
    }, numeric(1))
    v_floor <- v_floor[sort(names(v_floor))]
    output_codes <- unname(site_codes[names(v_floor)])
    yb <- vb <- rawv <- matrix(NA_real_, B, length(v_floor),
      dimnames = list(NULL, output_codes)
    )
    sm <- matrix(NA_real_, B, 12, dimnames = list(
      NULL,
      c(
        "I", "H", "V", "floor_sites", "zero_variance_sites", "min_raw_v",
        "max_floor_inflation", "direct_kl_check", "clamp_rows", "max_h",
        "clamped_treatment_score_max", "clamped_variance_fraction_max"
      )
    ))
    for (b in seq_len(B)) {
      idx <- unlist(lapply(split_cells, function(ii)
        ii[sample.int(length(ii), length(ii), replace = TRUE)]))
      db <- d[idx, ]
      db <- db[order(db$site), ]
      sb <- do.call(rbind, lapply(split(db, db$site), function(dd) site_effect_clamp(dd, outcome)))
      tau_b <- sb[, "est"]
      raw_v <- sb[, "se"]^2
      v_b <- pmax(raw_v, v_floor)
      fit <- fit_rubin_grid(tau_b, v_b, M = 6000L, seed = 20000 + b)
      ds <- summary_draws(fit, v_b)

      # Compare the two H formulas in the first replicate of each outcome.
      direct_error <- NA_real_
      if (b == 1L) {
        S <- outer(fit$sigma^2, v_b, function(aa, vv) aa / (aa + vv))
        ipd <- plogis(rowMeans(log(S) - log1p(-S)))
        Hdirect <- mean(vapply(seq_len(nrow(S)), function(m)
          mean(kl_bern(ipd[m], S[m, ])), numeric(1)))
        direct_error <- abs(Hdirect - ds["H"])
        stopifnot(direct_error < 1e-12)
      }
      sm[b, ] <- c(
        ds[c("I", "H", "V")], sum(raw_v < v_floor), sum(raw_v == 0),
        min(raw_v), max(v_b / pmax(raw_v, 1e-300)), direct_error,
        sum(sb[, "clamp_rows"]), max(sb[, "max_h"]),
        max(sb[, "clamped_treatment_score_max"]), max(sb[, "clamped_variance_fraction"])
      )
      yb[b, ] <- tau_b
      vb[b, ] <- v_b
      rawv[b, ] <- raw_v
      if (b %% 50L == 0L) cat(sprintf(
        "Bootstrap %s: %d/%d, elapsed %.1fs.\n",
        outcome, b, B, proc.time()[[3]] - t0
      ))
    }
    list(
      summary = data.frame(replicate = seq_len(B), sm),
      y = yb, v = vb, raw_v = rawv, floor = setNames(unname(v_floor), output_codes),
      final_rng = .Random.seed, method = "legacy250grid_M6000_originalRNGsequence",
      B = B, outcome = outcome
    )
  }
  bpass <- boot_run("scm_passcrs", run_B)
  bcred <- boot_run("credits_y4", run_B)
  saveRDS(
    list(
      passing = bpass, credits = bcred, initial_rng = boot_initial_state,
      participant_data_included = FALSE
    ),
    file.path(out_dir, paste0(prefix, "bootstrap-replication.rds"))
  )
  write_csv(bpass$summary, paste0(prefix, "bootstrap-passing-replicates.csv"))
  write_csv(bcred$summary, paste0(prefix, "bootstrap-credits-replicates.csv"))
  summ <- do.call(rbind, lapply(c("passing", "credits"), function(nm) {
    z <- if (nm == "passing") bpass$summary else bcred$summary
    do.call(rbind, lapply(c("I", "H", "V"), function(metric)
      data.frame(
        outcome = nm, metric = metric, B = nrow(z), mean = mean(z[[metric]]),
        median = median(z[[metric]]), lo = unname(quantile(z[[metric]], .025)),
        hi = unname(quantile(z[[metric]], .975))
      )))
  }))
  write_csv(summ, paste0(prefix, "bootstrap-summary.csv"))
  floor_summary <- do.call(rbind, lapply(c("passing", "credits"), function(nm) {
    z <- if (nm == "passing") bpass else bcred
    hit <- sweep(z$raw_v, 2, z$floor, `<`)
    data.frame(
      outcome = nm, B = nrow(hit), site_fits = length(hit), floor_hits = sum(hit),
      fraction_site_fits = sum(hit) / length(hit),
      replicates_any_floor = sum(rowSums(hit) > 0),
      max_sites_per_replicate = max(rowSums(hit)),
      exact_zero_variances = sum(z$raw_v == 0)
    )
  }))
  write_csv(floor_summary, paste0(prefix, "bootstrap-floor-summary.csv"))
  if (!pilot) {
    qs <- function(x) sprintf("[%.3f, %.3f]", quantile(x, .025), quantile(x, .975))
    now <- data.frame(
      outcome = c("passing", "credits"),
      I_boot95 = c(qs(bpass$summary$I), qs(bcred$summary$I)),
      H_boot95 = c(qs(bpass$summary$H), qs(bcred$summary$H)),
      V_boot95 = c(qs(bpass$summary$V), qs(bcred$summary$V))
    )
    old <- saved$sens$firststage
    matchrows <- data.frame(
      outcome = now$outcome,
      I_matches = now$I_boot95 == old$I_boot95,
      H_matches = now$H_boot95 == old$H_boot95,
      V_matches = now$V_boot95 == old$V_boot95
    )
    write_csv(now, "bootstrap-formatted.csv")
    write_csv(matchrows, "bootstrap-saved-comparison.csv")
    cat("Bootstrap intervals match the saved results: ", all(as.matrix(matchrows[, -1])), ".\n", sep = "")
  }
  write_json(list(
    mode = if (pilot) "pilot" else "full_sensitivity_replication", B = run_B,
    seconds = proc.time()[[3]] - t0, legacy_grid = 250L, bootstrap_M = 6000L,
    participant_rows_saved = FALSE
  ), paste0(prefix, "sensitivity-runtime.json"))
  cat(sprintf("Sensitivity mode finished in %.1fs.\n", proc.time()[[3]] - t0))
}

if (any(c("--reference-comparisons", "--aggregate-reference") %in% args)) {
  # Integrate over sigma in [0, Inf), with a flat prior on tau and sigma.
  # The supplied y and v stay fixed within each posterior calculation.
  source(file.path(analysis_dir, "scripts", "05_posterior_stability.R"))
  t0 <- proc.time()[[3]]
  specs <- list(
    block_adjusted_RD = list(y = sites$rd, v = sites$v),
    unadjusted_RD = list(y = sites$rd_unadjusted, v = sites$v_unadjusted),
    unadjusted_log_OR = list(y = sites$log_or, v = sites$v_log_or)
  )
  comp <- lapply(names(specs), function(nm) {
    x <- specs[[nm]]
    ref <- posterior_reference(x$y, x$v, TRUE, TRUE)
    data.frame(
      specification = nm, J = length(x$y), GM_v = exp(mean(log(x$v))),
      observed_y_variance = var(x$y), t(ref$summary), check.names = FALSE
    )
  })
  write_csv(do.call(rbind, comp), "firststage-scale-reference.csv")

  # Reuse the same bootstrap samples to isolate the change in numerical integration.
  boot <- readRDS(file.path(out_dir, "bootstrap-replication.rds"))
  boot_ref <- list()
  floor_ref <- list()
  floor_effects <- list()
  metrics <- c("I.mean", "H.mean", "V")
  for (nm in c("passing", "credits")) {
    x <- boot[[nm]]
    ref <- t(vapply(seq_len(nrow(x$y)), function(b) {
      posterior_reference(x$y[b, ], x$v[b, ], FALSE)$summary
    }, posterior_reference(x$y[1, ], x$v[1, ], FALSE)$summary))
    boot_ref[[nm]] <- data.frame(replicate = seq_len(nrow(ref)), ref, check.names = FALSE)
    write_csv(boot_ref[[nm]], paste0("bootstrap-reference-", nm, ".csv"))
    base <- ref[, metrics, drop = FALSE]
    for (mult in c(.5, 1, 2)) {
      altered_v <- pmax(x$raw_v, matrix(mult * x$floor, nrow(x$raw_v), ncol(x$raw_v), byrow = TRUE))
      changed <- rowSums(abs(altered_v - x$v) > 1e-14) > 0
      alt <- base
      for (b in which(changed)) {
        alt[b, ] <- posterior_reference(x$y[b, ], altered_v[b, ], FALSE)$summary[metrics]
      }
      floor_ref[[paste(nm, mult, sep = "_")]] <- data.frame(
        outcome = nm,
        floor_multiplier = mult, replicate = seq_len(nrow(alt)), alt, changed = changed,
        check.names = FALSE
      )
      floor_effects[[paste(nm, mult, sep = "_")]] <- do.call(rbind, lapply(metrics, function(metric)
        data.frame(
          outcome = nm, floor_multiplier = mult, metric = metric,
          changed_replicates = sum(changed),
          mean_difference = mean(alt[, metric] - base[, metric]),
          max_absolute_difference = max(abs(alt[, metric] - base[, metric])),
          mean = mean(alt[, metric]), median = median(alt[, metric]),
          lo = unname(quantile(alt[, metric], .025)),
          hi = unname(quantile(alt[, metric], .975))
        )))
    }
    cat(sprintf(
      "Reference bootstrap/floor %s complete; elapsed %.1fs.\n", nm,
      proc.time()[[3]] - t0
    ))
  }
  write_csv(do.call(rbind, floor_ref), "bootstrap-floor-reference-replicates.csv")
  write_csv(do.call(rbind, floor_effects), "bootstrap-floor-reference-summary.csv")
  boot_summary <- do.call(rbind, lapply(names(boot_ref), function(nm) {
    x <- boot_ref[[nm]]
    do.call(rbind, lapply(metrics, function(metric) data.frame(
      outcome = nm, metric = metric,
      mean = mean(x[[metric]]), median = median(x[[metric]]),
      lo = unname(quantile(x[[metric]], .025)), hi = unname(quantile(x[[metric]], .975))
    )))
  }))
  write_csv(boot_summary, "bootstrap-reference-summary.csv")
  old_boot <- read.csv(file.path(out_dir, "bootstrap-summary.csv"))
  old_boot$metric <- ifelse(old_boot$metric == "V", "V", paste0(old_boot$metric, ".mean"))
  merged <- merge(boot_summary, old_boot, by = c("outcome", "metric"), suffixes = c("_reference", "_legacy"))
  for (col in c("mean", "median", "lo", "hi"))
    merged[[paste0(col, "_change")]] <- merged[[paste0(col, "_reference")]] - merged[[paste0(col, "_legacy")]]
  write_csv(merged, "bootstrap-reference-legacy-comparison.csv")
  sub <- readRDS(file.path(out_dir, "subset-replication.rds"))
  subref <- do.call(rbind, lapply(seq_len(nrow(sub$keep)), function(k) {
    keep <- sub$keep[k, ]
    pp <- posterior_reference(sites$rd[keep], sites$v[keep], FALSE)$summary
    cc <- posterior_reference(sites$credits[keep], sites$v_credits[keep], FALSE)$summary
    data.frame(subset = k, J = sum(keep), I_pass = pp["I.mean"], I_credits = cc["I.mean"])
  }))
  write_csv(subref, "subset-reference-replicates.csv")
  sub_summary <- data.frame(
    outcome = c("passing", "credits"),
    I_min = c(min(subref$I_pass), min(subref$I_credits)),
    I_med = c(median(subref$I_pass), median(subref$I_credits)),
    I_max = c(max(subref$I_pass), max(subref$I_credits))
  )
  write_csv(sub_summary, "subset-reference-summary.csv")
  saveRDS(list(
    bootstrap = boot_ref, subset = subref, method = "untruncated deterministic quadrature",
    empirical_resamples = "same saved bootstrap and subset samples",
    floor_multipliers = c(.5, 1, 2)
  ), file.path(out_dir, "reference-sensitivity.rds"))
  write_json(list(
    seconds = proc.time()[[3]] - t0, method = "shared posterior_reference",
    bootstrap_refits_per_outcome = 500L, subset_selections = 256L,
    floor_multipliers = c(.5, 1, 2), no_new_resampling = TRUE
  ), "reference-comparisons-runtime.json")
  cat(sprintf("Reference comparisons complete in %.1fs.\n", proc.time()[[3]] - t0))
}

if ("--final-checks" %in% args) {
  source(file.path(analysis_dir, "scripts", "05_posterior_stability.R"))

  # Check sigma and log-sigma integration for each effect-scale specification.
  specs <- list(
    block_adjusted_RD = list(y = sites$rd, v = sites$v),
    unadjusted_RD = list(y = sites$rd_unadjusted, v = sites$v_unadjusted),
    unadjusted_log_OR = list(y = sites$log_or, v = sites$v_log_or)
  )
  qchecks <- lapply(names(specs), function(nm) {
    x <- specs[[nm]]
    rr <- posterior_reference(x$y, x$v, FALSE, TRUE)
    ch <- rr$checks
    stopifnot(
      ch["normalization_relative_error"] < 1e-8,
      abs(ch["normalization_logsigma_relative_difference"]) < 1e-8,
      max(abs(ch[grepl("logsigma_minus_sigma", names(ch))])) < 1e-8
    )
    list(specification = nm, checks = as.list(ch))
  })
  write_json(qchecks, "firststage-scale-quadrature-checks.json")

  # Check which linearly dependent covariates were dropped by the QR decomposition.
  dropped <- character()
  for (dd in site_groups) {
    X <- cbind(`(Intercept)` = 1, ragrp = dd$ragrp)
    if (length(unique(dd$rablock)) > 1L)
      X <- cbind(X, model.matrix(~ factor(rablock), dd)[, -1, drop = FALSE])
    pre <- if (dd$subject[1] == "Algebra 1") dd$state_math_imp else dd$state_ela_imp
    X <- cbind(X, pre = pre, as.matrix(dd[, c("state_miss", "gpa9th_imp", "gpa9th_miss")]))
    qx <- qr(X)
    if (qx$rank < ncol(X)) dropped <- c(dropped, colnames(X)[qx$pivot[(qx$rank + 1):ncol(X)]])
  }
  write_csv(as.data.frame(table(dropped), stringsAsFactors = FALSE), "adjusted-aliased-columns.csv")
  clamp_summary <- do.call(rbind, lapply(unique(mats$specification), function(s) {
    z <- mats[mats$specification == s, ]
    data.frame(
      analysis = s, replicates = 1, clamp_rows_per_rep_min = sum(z$clamp_rows),
      clamp_rows_per_rep_max = sum(z$clamp_rows), max_h = max(z$max_h),
      max_abs_clamped_treatment_score = max(z$clamped_treatment_score_max),
      max_clamped_variance_fraction = max(z$clamped_variance_fraction)
    )
  }))
  boot <- readRDS(file.path(out_dir, "bootstrap-replication.rds"))
  for (nm in c("passing", "credits")) {
    z <- boot[[nm]]$summary
    clamp_summary <- rbind(clamp_summary, data.frame(
      analysis = paste0("bootstrap_", nm),
      replicates = nrow(z), clamp_rows_per_rep_min = min(z$clamp_rows),
      clamp_rows_per_rep_max = max(z$clamp_rows), max_h = max(z$max_h),
      max_abs_clamped_treatment_score = max(z$clamped_treatment_score_max),
      max_clamped_variance_fraction = max(z$clamped_variance_fraction_max)
    ))
  }
  write_csv(clamp_summary, "hc2-clamp-summary.csv")
  write_json(
    list(
      subject_matched_English_max = max(d$credits_y4[d$subject == "English 9"]),
      subject_matched_Algebra_max = max(d$credits_y4[d$subject == "Algebra 1"]),
      entire_file_English_field_max = max(d$cby_y4_eng),
      correction = "The analysis uses credits for the corresponding subject; the English maximum is 85."
    ),
    "credits-range-correction.json"
  )
  deps <- c(
    raw_path, file.path(in_dir, "ocrs_fits.rds"),
    file.path(analysis_dir, "scripts", "legacy", c("00_utils.R", "04_ocrs.R")),
    file.path(analysis_dir, "scripts", "05_posterior_stability.R"), script_path
  )
  dep_hash <- vapply(deps, digest::digest, character(1), file = TRUE, algo = "sha256")
  dep_paths <- c(
    "external/OCRS_publicuse_impactanalysisdata.csv",
    sub(paste0(dirname(analysis_dir), "/"), "", deps[-1], fixed = TRUE)
  )
  write_json(data.frame(path = dep_paths, sha256 = unname(dep_hash)), "dependencies-sha256.json")
  stopifnot(identical(actual_raw_sha256, raw_sha256))
  write_json(list(
    raw_input_sha256_verified = TRUE,
    raw_participant_rows_written = FALSE,
    exact_target_claim = "fixed supplied-input regression contrast; no all-student SATE or latent coverage claim",
    protected_sources_modified = FALSE
  ), "final-checks.json")
  cat("Quadrature, HC2 variance, input checksum, and aggregate-output checks passed.\n")
}
