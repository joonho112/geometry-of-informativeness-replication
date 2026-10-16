# Rerun theory checks or simulations in a new working directory.
# Rscript analysis/scripts/run_historical.R --arm=theory
# Rscript analysis/scripts/run_historical.R --arm=stress --smoke
# Rscript analysis/scripts/run_historical.R --arm=stress --spot-check
# Rscript analysis/scripts/run_historical.R --arm=fixedN --cores=2
# Full defaults preserve 500 replications and 2,000 posterior draws.
args <- commandArgs(trailingOnly = TRUE)
value <- function(prefix, default) {
  z <- args[startsWith(args, prefix)]
  if (length(z) > 1L) stop("Duplicate argument: ", prefix)
  if (length(z)) substring(z, nchar(prefix) + 1L) else default
}
known <- grepl("^--(arm|cores|workdir)=", args) | args %in% c("--smoke", "--spot-check")
if (any(!known)) stop("Unknown arguments: ", paste(args[!known], collapse = ", "))
arm <- value("--arm=", "theory")
stopifnot(arm %in% c("theory", "pretest", "stress", "fixedN"))
smoke <- "--smoke" %in% args
spot <- "--spot-check" %in% args
if ((smoke || spot) && !arm %in% c("stress", "fixedN")) stop("Smoke/spot modes are for stress/fixedN")
if (smoke && spot) stop("Choose either --smoke or --spot-check")
cores <- as.integer(value("--cores=", "1"))
stopifnot(length(cores) == 1L, !is.na(cores), cores >= 1L)
if (.Platform$OS.type == "windows" && cores > 1L) stop("Use --cores=1 on Windows")
filearg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", filearg[[1L]]), mustWork = TRUE)
analysis <- dirname(dirname(script))
tag <- paste0(arm, if (smoke) "-smoke" else if (spot) "-spot" else "", "-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid())
work <- value("--workdir=", file.path(analysis, "generated", tag))
if (file.exists(work)) stop("Choose a new, nonexistent work directory: ", work)
if (any(strsplit(work, "[/\\\\]")[[1L]] == "..")) stop("Use a workdir without '..' components")
work <- if (grepl("^/|^[A-Za-z]:", work)) work else file.path(getwd(), work)
frozen <- normalizePath(file.path(analysis, "inputs"), mustWork = TRUE)
if (startsWith(paste0(work, "/"), paste0(frozen, "/"))) stop("Workdir must not be within analysis/inputs")
dir.create(work, recursive = TRUE)
work <- normalizePath(work, mustWork = TRUE)
if (startsWith(paste0(work, "/"), paste0(frozen, "/"))) stop("Workdir must not be within analysis/inputs")
name <- switch(arm,
  theory = "01_verification.R",
  pretest = "02_sim_pretest.R",
  stress = "03_sim_main.R",
  fixedN = "08_sim_fixedN.R"
)
stopifnot(all(file.copy(file.path(analysis, "scripts", "legacy", c("00_utils.R", name)), work)))
stopifnot(file.copy(file.path(analysis, "scripts", "dependencies", "historical-dgp.R"), work))
code <- readLines(file.path(work, name), warn = FALSE)
# Each family has its own seed, so changing the worker count preserves its draws.
code <- sub("mc.cores = max(1L, detectCores() - 2L)", paste0("mc.cores = ", cores, "L"), code, fixed = TRUE)
if (smoke) {
  code <- sub("REPS    <- 500L", "REPS    <- 2L", code, fixed = TRUE)
  code <- sub("M_DRAWS <- 2000L", "M_DRAWS <- 200L", code, fixed = TRUE)
  # Keep the original family and profile IDs used to set random seeds.
  family_ids <- if (arm == "stress") "c(1L, 19L, 37L)" else "c(1L, 10L, 19L)"
  code <- sub("mclapply(seq_len(nrow(families)), run_family,",
    paste0("mclapply(", family_ids, ", run_family,"), code,
    fixed = TRUE
  )
  # Use two profiles in each selected family, covering all three effect laws.
  if (arm == "stress") {
    at <- which(startsWith(code, 'cat("families:'))[1L]
    code <- append(code, "combos <- combos[1:2, , drop = FALSE]", after = at - 1L)
  } else code <- sub("CVs <- c(0, 0.25, 0.5, 0.75)", "CVs <- c(0, 0.25)", code, fixed = TRUE)
  # Label the printed count as the full grid rather than the reduced run.
  code <- sub('"families:"', '"full-grid families (smoke executes 3):"', code, fixed = TRUE)
}
if (spot) {
  # Retain the original latent draws and fit settings, but fit only the
  # first/last profiles and replications 1/500 in three families.
  family_ids <- if (arm == "stress") "c(1L, 19L, 37L)" else "c(1L, 10L, 19L)"
  code <- sub("mclapply(seq_len(nrow(families)), run_family,",
    paste0("mclapply(", family_ids, ", run_family,"), code,
    fixed = TRUE
  )
  code <- sub("for (r in seq_len(REPS))", "for (r in c(1L, 500L))", code, fixed = TRUE)
  if (arm == "stress") code <- sub("for (ci in seq_len(nrow(combos)))",
    "for (ci in c(1L, nrow(combos)))", code,
    fixed = TRUE
  )
  else code <- sub("for (ci in seq_along(CVs))", "for (ci in c(1L, length(CVs)))", code, fixed = TRUE)
  code <- sub('"families:"', '"full-grid families (spot-check executes 3):"', code, fixed = TRUE)
}
writeLines(code, file.path(work, name))
dir.create(file.path(work, "outputs"))
writeLines(
  c(
    paste("arm:", arm), paste("smoke:", smoke), paste("spot_check:", spot), paste("cores:", cores),
    "supplied inputs are never overwritten",
    "Smoke uses 3 families x 2 profiles x 2 replications, M=200; not manuscript results.",
    "Spot-check uses original latent draws/M=2000; 12 fits per arm at original profile and replicate IDs.",
    "Full stress/fixedN defaults: 500 replications, M=2000, grid sampler with uniform jitter."
  ),
  file.path(work, "RUN-SCOPE.txt")
)
cat("Run directory:", work, "\n")
setwd(work)
runenv <- new.env(parent = globalenv())
source(name, local = runenv)
if (spot) {
  saved <- readRDS(file.path(frozen, if (arm == "stress") "sim_main_reps.rds" else "sim_fixedN_reps.rds"))
  key <- c("family", if (arm == "stress") "combo" else "cv_i", "rep", "est")
  m <- merge(runenv$res, saved, by = key, suffixes = c(".new", ".saved"))
  stopifnot(nrow(m) == nrow(runenv$res), !anyDuplicated(m[key]))
  cols <- c(
    "rmse", "mselp", "isel", "I", "H", "V_logit", "J_eff", "S_bar",
    if (arm == "stress") "I_hat"
  )
  err <- vapply(cols, function(nm) max(abs(m[[paste0(nm, ".new")]] -
    m[[paste0(nm, ".saved")]])), numeric(1))
  stopifnot(all(is.finite(err)), max(err) < 1e-12)
  write.csv(data.frame(quantity = names(err), maximum_absolute_difference = err),
    "outputs/spot-comparison.csv",
    row.names = FALSE
  )
  cat("Saved-result comparison:", nrow(m), "estimator rows; max abs difference", max(err), "\n")
}
if (arm == "theory") {
  source(file.path(analysis, "scripts", "01_effective_count_checks.R"), local = new.env(parent = globalenv()))
}
cat("Completed", arm, if (smoke) "smoke" else if (spot) "spot-check" else "run", "in", work, "\n")
