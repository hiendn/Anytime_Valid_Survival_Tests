# Source-safe SIM-2 comparators. No data generation, driver, or file writes.
# Primary reference: ter Schure et al. (2024), §§3.2--3.3, (3.7)--(3.8),
# https://nejsds.nestat.org/journal/NEJSDS/article/73/text
# Exact-event unrestricted prequential: the paper's two-virtual-event smoothed
# MLE, fitted only to earlier release times. Grouped likelihood fitting is an
# explicitly identified extension. Capping its predictor at theta=1 is a
# separate, direction-constrained adaptation, not the original procedure.
.sim2_script_file <- local({
  paths <- vapply(sys.frames(), function(e) if (is.null(e$ofile)) "" else as.character(e$ofile), "")
  paths <- paths[nzchar(paths)]
  if (!length(paths)) stop("Source scripts/sim2_comparators.R to load these functions")
  normalizePath(tail(paths, 1L), mustWork = TRUE)
})
.sim2_script_dir <- dirname(.sim2_script_file)
.sim2_kernels_ready <- FALSE

sim2_ensure_kernels <- function() {
  if (isTRUE(.sim2_kernels_ready) && exists("sim2_comparators_cpp", mode = "function")) return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required")
  cache <- Sys.getenv("SIMULATION_RCPP_CACHE_DIR", unset = file.path(tempdir(), "sim2_rcpp_cache"))
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  old <- Sys.getenv("PKG_CPPFLAGS", unset = "")
  on.exit(Sys.setenv(PKG_CPPFLAGS = old), add = TRUE)
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    sdk <- tryCatch(system2("xcrun", "--show-sdk-path", stdout = TRUE, stderr = FALSE), error = function(e) character())
    if (length(sdk)) {
      headers <- file.path(sdk[[1L]], "usr", "include", "c++", "v1")
      if (dir.exists(headers)) Sys.setenv(PKG_CPPFLAGS = paste(old, "-isystem", shQuote(headers)))
    }
  }
  Rcpp::sourceCpp(file.path(.sim2_script_dir, "sim2_kernels.cpp"), cacheDir = cache,
                  showOutput = FALSE, verbose = FALSE, env = .GlobalEnv)
  .sim2_kernels_ready <<- TRUE
  invisible(TRUE)
}

# One fit per release produces both unrestricted and protective variants.
# initial_n1/n0 are the ENROLLED counts, not the first (possibly depleted) risk
# set. Count restricts each packed row to its real cells; padded cells are ignored.
sim2_comparators <- function(updates, initial_n1, initial_n0, horizon,
                             eta = log(0.7), alpha = 0.05, keep_path = FALSE) {
  required <- c("n1", "n0", "d", "x", "time", "count")
  if (!all(required %in% names(updates))) stop("Missing packed update fields")
  shape <- dim(updates$n1)
  if (length(shape) != 2L || any(vapply(updates[required[1:5]], function(x) !identical(dim(x), shape), TRUE)))
    stop("Packed matrices must have identical dimensions")
  if (length(updates$count) != shape[[1L]] || anyNA(updates$count) ||
      any(!is.finite(updates$count)) || any(updates$count != as.integer(updates$count)) ||
      any(updates$count < 0 | updates$count > shape[[2L]]))
    stop("count must have one finite integer in [0,ncol] per dataset")
  expand <- function(x, name) {
    if (length(x) == 1L) x <- rep(x, shape[[1L]])
    if (length(x) != shape[[1L]] || anyNA(x) || any(!is.finite(x)) || any(x < 0) || any(x != as.integer(x)))
      stop(name, " must contain nonnegative integer enrolled counts")
    as.integer(x)
  }
  n1 <- expand(initial_n1, "initial_n1")
  n0 <- expand(initial_n0, "initial_n0")
  if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) stop("horizon must be positive")
  if (length(eta) != 1L || !is.finite(eta)) stop("eta must be finite")
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must lie in (0,1)")
  sim2_ensure_kernels()
  out <- sim2_comparators_cpp(updates$n1, updates$n0, updates$d, updates$x,
                              updates$time, as.integer(updates$count), n1, n0,
                              as.numeric(horizon), as.numeric(eta), as.numeric(alpha), isTRUE(keep_path))
  out$fixed$method <- "original_fixed_alternative_unit_weight"
  out$prequential$method <- "unrestricted_smoothed_prequential"
  out$prequential_protective$method <- "protective_smoothed_prequential_adaptation"
  grouped <- vapply(seq_len(shape[[1L]]), function(i)
    any(updates$d[i, seq_len(updates$count[[i]])] > 1L), TRUE)
  out$prequential$has_grouped_cells <- grouped
  out$prequential_protective$has_grouped_cells <- grouped
  out
}

sim2_fixed <- function(updates, eta = log(0.7), horizon, alpha = 0.05, keep_path = FALSE) {
  # Initial counts do not enter the fixed numerator. Any finite positive values
  # suffice for the unused simultaneous prequential calculation in this wrapper.
  sim2_comparators(updates, 1L, 1L, horizon, eta, alpha, keep_path)$fixed
}

sim2_prequential <- function(updates, initial_n1, initial_n0, horizon,
                             alpha = 0.05, direction = c("unrestricted", "protective"),
                             keep_path = FALSE) {
  direction <- match.arg(direction)
  out <- sim2_comparators(updates, initial_n1, initial_n0, horizon,
                          alpha = alpha, keep_path = keep_path)
  if (direction == "protective") out$prequential_protective else out$prequential
}
