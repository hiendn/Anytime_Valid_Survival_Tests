#!/usr/bin/env Rscript

# Covariate-stratification and delayed-entry simulation for the manuscript.
# Only the main-text LR comparison is run.
#
# This file is deliberately source-safe: sourcing it defines the panel-engine
# API but does not run the study.  The hot Fisher noncentral-hypergeometric
# likelihood-ratio recurrence is provided by simulation_core.R through Rcpp.

panel_script_path <- function() {
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      ofile <- frames[[i]]$ofile
      if (!is.null(ofile) && nzchar(ofile)) {
        return(normalizePath(ofile, mustWork = FALSE))
      }
    }
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE), mustWork = FALSE))
  }
  normalizePath("run_covariate_study.R", mustWork = FALSE)
}

PANEL_SCRIPT_PATH <- panel_script_path()
PANEL_CODE_DIR <- dirname(PANEL_SCRIPT_PATH)
source(file.path(PANEL_CODE_DIR, "simulation_config.R"), local = TRUE)
PANEL_PROJECT_ROOT <- dirname(PANEL_CODE_DIR)
PANEL_CORE_PATH <- file.path(PANEL_CODE_DIR, "simulation_core.R")

if (!file.exists(PANEL_CORE_PATH)) {
  stop("Shared simulation core not found: ", PANEL_CORE_PATH, call. = FALSE)
}

# Keep the core in a private namespace.  In particular, its mixed-resolution
# SCENARIOS registry must not be overwritten by the panel registry below.
.panel_core <- new.env(parent = globalenv())
sys.source(PANEL_CORE_PATH, envir = .panel_core)
# sys.source() adds both this driver and the sourced file to R's call stack.
# The core's intentionally defensive path discovery can therefore see this
# driver first on some R releases.  Pin the private core namespace to the
# directory that we have already resolved and verified above.
.panel_core$SIMULATION_CODE_DIR <- PANEL_CODE_DIR

panel_required_core_api <- c(
  "table_loge", "hypergeom_mean_var", "row_cumsums",
  "ensure_simulation_kernels", "configure_simulation_rng",
  "collect_r_metadata", "write_result_csv", "write_result_json"
)
panel_missing_core_api <- panel_required_core_api[
  !vapply(panel_required_core_api, exists, logical(1), envir = .panel_core,
          inherits = FALSE)
]
if (length(panel_missing_core_api)) {
  stop(
    "simulation_core.R is missing required API: ",
    paste(panel_missing_core_api, collapse = ", "),
    call. = FALSE
  )
}

ALPHA <- 0.05
LOG_BOUNDARY <- log(1 / ALPHA)
WEIGHT_NAMES <- "logrank"
DISPLAY_NAMES <- c(logrank = "Logrank / CMH")

new_scenario <- function(
    name,
    description,
    n,
    visits,
    base_rate,
    cure_control = 0,
    cure_treatment = 0,
    hr_early = 1,
    hr_late = 1,
    change_time = Inf,
    censor_rate = 0.04,
    covariate = FALSE,
    confounding = FALSE,
    left_truncation = FALSE,
    covariate_hr = 1) {
  structure(
    list(
      name = as.character(name),
      description = as.character(description),
      n = as.integer(n),
      visits = as.numeric(visits),
      base_rate = as.numeric(base_rate),
      cure_control = as.numeric(cure_control),
      cure_treatment = as.numeric(cure_treatment),
      hr_early = as.numeric(hr_early),
      hr_late = as.numeric(hr_late),
      change_time = as.numeric(change_time),
      censor_rate = as.numeric(censor_rate),
      covariate = isTRUE(covariate),
      confounding = isTRUE(confounding),
      left_truncation = isTRUE(left_truncation),
      covariate_hr = as.numeric(covariate_hr)
    ),
    class = "panel_scenario"
  )
}

PANEL_SCENARIOS <- list(
  covariate_null = new_scenario(
    "covariate_null",
    paste(
      "No treatment effect; binary prognostic covariate, confounding,",
      "censoring and delayed entry."
    ),
    n = REFERENCE_N, visits = 0:8, base_rate = 0.17,
    cure_control = 0.12, cure_treatment = 0.12,
    censor_rate = 0.05, covariate = TRUE, confounding = TRUE,
    left_truncation = TRUE, covariate_hr = 2.6
  ),
  covariate_alt = new_scenario(
    "covariate_alt",
    paste(
      "Hazard ratio 0.65 within binary prognostic strata, with",
      "confounding, censoring and delayed entry."
    ),
    n = REFERENCE_N, visits = 0:8, base_rate = 0.17,
    cure_control = 0.12, cure_treatment = 0.12,
    hr_early = 0.65, hr_late = 0.65,
    censor_rate = 0.05, covariate = TRUE, confounding = TRUE,
    left_truncation = TRUE, covariate_hr = 2.6
  )
)

# Consumers that also need the mixed engine should use a dedicated environment.
SCENARIOS <- PANEL_SCENARIOS

panel_rng_kind <- function() {
  .panel_core$configure_simulation_rng()
}

new_seed_registry <- function() new.env(parent = emptyenv())

derived_seed <- function(base_seed, offset, label = NULL, registry = NULL) {
  modulus <- 2147483646
  value <- as.integer((as.double(base_seed) + as.double(offset)) %% modulus)
  if (is.na(value) || value <= 0L) value <- value + modulus
  if (!is.null(registry) && !is.null(label)) {
    assign(as.character(label), value, envir = registry)
  }
  value
}

seed_registry_list <- function(registry) {
  keys <- sort(ls(registry, all.names = TRUE))
  stats::setNames(lapply(keys, function(k) get(k, envir = registry)), keys)
}

set_panel_seed <- function(seed) {
  panel_rng_kind()
  set.seed(as.integer(seed))
  invisible(seed)
}

.draw_treatment <- function(reps, n, z = NULL, confounding = FALSE) {
  if (isTRUE(confounding)) {
    if (is.null(z)) stop("confounding requires z", call. = FALSE)
    p <- ifelse(z == 1L, 0.28, 0.72)
    return(matrix(as.integer(runif(reps * n) < p), nrow = reps, ncol = n))
  }
  matrix(as.integer(runif(reps * n) < 0.5), nrow = reps, ncol = n)
}

.simulate_event_time <- function(scenario, entry, group, z, cured) {
  reps <- nrow(entry)
  n <- ncol(entry)
  e <- matrix(rexp(reps * n, rate = 1), nrow = reps, ncol = n)
  cov_mult <- ifelse(z == 1L, scenario$covariate_hr, 1)
  base <- scenario$base_rate * cov_mult

  if (is.infinite(scenario$change_time) ||
      abs(scenario$hr_early - scenario$hr_late) < 1e-14) {
    mult <- ifelse(group == 1L, scenario$hr_early, 1)
    event_time <- entry + e / (base * mult)
  } else {
    change <- scenario$change_time
    before <- entry < change
    rate1 <- base * ifelse(group == 1L, scenario$hr_early, 1)
    rate2 <- base * ifelse(group == 1L, scenario$hr_late, 1)
    available <- pmax(change - entry, 0)
    h_to_change <- rate1 * available
    fail_before <- before & (e <= h_to_change)
    event_time <- matrix(NA_real_, nrow = reps, ncol = n)
    event_time[fail_before] <-
      entry[fail_before] + e[fail_before] / rate1[fail_before]
    after <- !fail_before
    residual <- pmax(e[after] - h_to_change[after], 0)
    start <- ifelse(before[after], change, entry[after])
    event_time[after] <- start + residual / rate2[after]
  }
  event_time[cured] <- Inf
  event_time
}

.array_slice_matrix <- function(x, block, stratum = NULL) {
  d <- dim(x)
  if (length(d) == 3L) {
    if (is.null(stratum)) {
      return(matrix(x[, block, , drop = FALSE], nrow = d[[1L]], ncol = d[[3L]]))
    }
    return(matrix(x[, block, stratum, drop = FALSE], nrow = d[[1L]], ncol = 1L))
  }
  matrix(x[, block, drop = FALSE], nrow = d[[1L]], ncol = 1L)
}

.reps_blocks_matrix <- function(x, reps = dim(x)[[1L]], blocks = dim(x)[[2L]]) {
  matrix(x, nrow = reps, ncol = blocks)
}

.sum_over_strata <- function(x) {
  d <- dim(x)
  stopifnot(length(d) == 3L)
  flat <- matrix(x, nrow = d[[1L]] * d[[2L]], ncol = d[[3L]])
  matrix(rowSums(flat), nrow = d[[1L]], ncol = d[[2L]])
}

simulate_panel_tables <- function(scenario, reps, seed, stratified = NULL) {
  reps <- as.integer(reps)
  if (reps < 1L) stop("reps must be positive", call. = FALSE)
  set_panel_seed(seed)

  n <- scenario$n
  visits <- scenario$visits
  blocks <- length(visits) - 1L
  if (is.null(stratified)) stratified <- scenario$covariate
  strata <- if (isTRUE(stratified) && scenario$covariate) 2L else 1L

  z <- if (scenario$covariate) {
    matrix(as.integer(runif(reps * n) < 0.5), nrow = reps, ncol = n)
  } else {
    matrix(0L, nrow = reps, ncol = n)
  }
  group <- .draw_treatment(
    reps, n, if (scenario$covariate) z else NULL, scenario$confounding
  )

  entry <- if (scenario$left_truncation) {
    u <- matrix(runif(reps * n), nrow = reps, ncol = n)
    ifelse(u < 0.55, 0, ifelse(u < 0.85, 1, 2))
  } else {
    matrix(0, nrow = reps, ncol = n)
  }

  cure_prob <- ifelse(group == 1L, scenario$cure_treatment, scenario$cure_control)
  # Ordinary left truncation: enrolled subjects are sampled conditional on
  # survival to entry from one origin-time cure model. This keeps the null
  # event law identical within Z across entry cohorts.
  if (scenario$left_truncation) {
    if (scenario$hr_early != scenario$hr_late) stop("Delayed entry currently requires constant susceptible hazards")
    rate <- scenario$base_rate * ifelse(z == 1L, scenario$covariate_hr, 1) *
      ifelse(group == 1L, scenario$hr_early, 1)
    cure_prob <- cure_prob / (cure_prob + (1 - cure_prob) * exp(-rate * entry))
  }
  cured <- matrix(runif(reps * n), nrow = reps, ncol = n) < cure_prob
  event_time <- .simulate_event_time(scenario, entry, group, z, cured)
  censor_time <- entry + matrix(
    rexp(reps * n, rate = scenario$censor_rate),
    nrow = reps, ncol = n
  )

  dims <- c(reps, blocks, strata)
  n1 <- array(0L, dim = dims)
  n0 <- array(0L, dim = dims)
  d_events <- array(0L, dim = dims)
  x_events <- array(0L, dim = dims)
  raw_weights <- array(0, dim = c(reps, blocks, length(WEIGHT_NAMES)))

  active <- matrix(TRUE, nrow = reps, ncol = n)

  for (b in seq_len(blocks)) {
    t0 <- visits[[b]]
    t1 <- visits[[b + 1L]]
    current <- active & (entry <= t0 + 1e-12)
    complete <- current & (censor_time >= t1 - 1e-12)
    events <- complete & (event_time <= t1 + 1e-12)

    masks <- if (strata == 1L) {
      list(matrix(TRUE, nrow = reps, ncol = n))
    } else {
      list(z == 0L, z == 1L)
    }

    for (s in seq_len(strata)) {
      risk <- complete & masks[[s]]
      ev <- events & masks[[s]]
      n1[, b, s] <- as.integer(rowSums(risk & (group == 1L)))
      n0[, b, s] <- as.integer(rowSums(risk & (group == 0L)))
      x_events[, b, s] <- as.integer(rowSums(ev & (group == 1L)))
      d_events[, b, s] <- as.integer(rowSums(ev))
    }

    n_pool <- rowSums(.array_slice_matrix(n1 + n0, b))
    d_pool <- rowSums(.array_slice_matrix(d_events, b))
    raw_weights[, b, 1L] <- 1

    incomplete <- current & (censor_time < t1 - 1e-12)
    active <- active & !(events | incomplete)
  }

  list(
    n1 = n1,
    n0 = n0,
    d = d_events,
    x = x_events,
    raw_weights = raw_weights,
    n_subjects = as.integer(n), block_end_times = visits[-1L]
  )
}

collapse_strata <- function(tables) {
  collapse_one <- function(x) {
    summed <- .sum_over_strata(x)
    array(as.integer(summed), dim = c(nrow(summed), ncol(summed), 1L))
  }
  list(
    n1 = collapse_one(tables$n1),
    n0 = collapse_one(tables$n0),
    d = collapse_one(tables$d),
    x = collapse_one(tables$x),
    raw_weights = tables$raw_weights,
    n_subjects = tables$n_subjects, block_end_times = tables$block_end_times
  )
}

hypergeom_mean_var <- function(n1, n0, d) {
  dims <- dim(n1)
  if (is.null(dims) || length(dims) <= 2L) {
    return(.panel_core$hypergeom_mean_var(n1, n0, d))
  }
  if (length(dims) != 3L ||
      !identical(dim(n0), dims) || !identical(dim(d), dims)) {
    stop("Panel hypergeometric arrays must share reps x blocks x strata dimensions",
         call. = FALSE)
  }
  mean_array <- array(0, dim = dims)
  var_array <- array(0, dim = dims)
  for (s in seq_len(dims[[3L]])) {
    n1_s <- matrix(n1[, , s, drop = FALSE], nrow = dims[[1L]], ncol = dims[[2L]])
    n0_s <- matrix(n0[, , s, drop = FALSE], nrow = dims[[1L]], ncol = dims[[2L]])
    d_s <- matrix(d[, , s, drop = FALSE], nrow = dims[[1L]], ncol = dims[[2L]])
    moment <- .panel_core$hypergeom_mean_var(n1_s, n0_s, d_s)
    mean_array[, , s] <- moment$mean
    var_array[, , s] <- moment$var
  }
  list(mean = mean_array, var = var_array)
}

planned_scales_and_information <- function(null_tables) {
  hm <- hypergeom_mean_var(null_tables$n1, null_tables$n0, null_tables$d)
  var_block <- .sum_over_strata(hm$var)
  raw <- null_tables$raw_weights
  reps <- dim(raw)[[1L]]
  blocks <- dim(raw)[[2L]]
  scales <- rep(1, length(WEIGHT_NAMES))
  planned_info <- matrix(0, nrow = length(WEIGHT_NAMES), ncol = blocks)
  base_info <- mean(rowSums(var_block))

  for (j in seq_along(WEIGHT_NAMES)) {
    raw_j <- matrix(raw[, , j, drop = FALSE], nrow = reps, ncol = blocks)
    weighted_total <- mean(rowSums(raw_j^2 * var_block))
    scales[[j]] <- sqrt(max(weighted_total / max(base_info, 1e-12), 1e-12))
    w <- raw_j / scales[[j]]
    planned_info[j, ] <- cumsum(colMeans(w^2 * var_block))
  }
  names(scales) <- WEIGHT_NAMES
  rownames(planned_info) <- WEIGHT_NAMES
  list(scales = scales, planned_info = planned_info)
}

panel_table_loge <- function(n1, n0, d, x, u, tol = 1e-14) {
  as.numeric(.panel_core$table_loge(n1, n0, d, x, u, tol = tol))
}

loge_paths <- function(tables, scales, etas, stratified = TRUE) {
  n1 <- tables$n1
  n0 <- tables$n0
  d_events <- tables$d
  x_events <- tables$x
  raw <- tables$raw_weights
  reps <- dim(n1)[[1L]]
  blocks <- dim(n1)[[2L]]
  strata <- if (isTRUE(stratified)) dim(n1)[[3L]] else 1L
  paths <- list()

  for (j in seq_along(WEIGHT_NAMES)) {
    name <- WEIGHT_NAMES[[j]]
    w <- matrix(raw[, , j, drop = FALSE], nrow = reps, ncol = blocks) / scales[[j]]
    increments <- matrix(0, nrow = reps, ncol = blocks)
    eta <- etas[[name]]
    for (b in seq_len(blocks)) {
      u <- eta * w[, b]
      for (s in seq_len(strata)) {
        increments[, b] <- increments[, b] + panel_table_loge(
          n1[, b, s], n0[, b, s], d_events[, b, s], x_events[, b, s], u
        )
      }
    }
    paths[[name]] <- .panel_core$row_cumsums(increments)
  }
  paths
}

score_paths <- function(tables, scales, stratified = TRUE) {
  hm <- hypergeom_mean_var(tables$n1, tables$n0, tables$d)
  reps <- dim(tables$n1)[[1L]]
  blocks <- dim(tables$n1)[[2L]]
  strata <- if (isTRUE(stratified)) dim(tables$n1)[[3L]] else 1L
  U <- list()
  V <- list()
  Z <- list()

  for (j in seq_along(WEIGHT_NAMES)) {
    name <- WEIGHT_NAMES[[j]]
    w <- matrix(tables$raw_weights[, , j, drop = FALSE],
                nrow = reps, ncol = blocks) / scales[[j]]
    u_inc <- matrix(0, nrow = reps, ncol = blocks)
    v_inc <- matrix(0, nrow = reps, ncol = blocks)
    for (s in seq_len(strata)) {
      xs <- matrix(tables$x[, , s, drop = FALSE], nrow = reps, ncol = blocks)
      ms <- matrix(hm$mean[, , s, drop = FALSE], nrow = reps, ncol = blocks)
      vs <- matrix(hm$var[, , s, drop = FALSE], nrow = reps, ncol = blocks)
      u_inc <- u_inc + w * (xs - ms)
      v_inc <- v_inc + w^2 * vs
    }
    u_path <- .panel_core$row_cumsums(u_inc)
    v_path <- .panel_core$row_cumsums(v_inc)
    z_path <- matrix(0, nrow = reps, ncol = blocks)
    ok <- v_path > 1e-12
    z_path[ok] <- u_path[ok] / sqrt(v_path[ok])
    U[[name]] <- u_path
    V[[name]] <- v_path
    Z[[name]] <- z_path
  }
  list(U = U, V = V, Z = Z)
}

calibrate_eta_grid <- function(alt_tables, scales, eta_grid, stratified = TRUE) {
  n1 <- alt_tables$n1
  n0 <- alt_tables$n0
  d_events <- alt_tables$d
  x_events <- alt_tables$x
  raw <- alt_tables$raw_weights
  reps <- dim(n1)[[1L]]
  blocks <- dim(n1)[[2L]]
  strata <- if (isTRUE(stratified)) dim(n1)[[3L]] else 1L
  rows <- vector("list", length(WEIGHT_NAMES) * length(eta_grid))
  row_index <- 0L
  chosen <- stats::setNames(numeric(length(WEIGHT_NAMES)), WEIGHT_NAMES)

  for (j in seq_along(WEIGHT_NAMES)) {
    name <- WEIGHT_NAMES[[j]]
    w <- matrix(raw[, , j, drop = FALSE], nrow = reps, ncol = blocks) / scales[[j]]
    means <- numeric(length(eta_grid))
    for (k in seq_along(eta_grid)) {
      eta <- eta_grid[[k]]
      total <- numeric(reps)
      for (b in seq_len(blocks)) {
        u <- eta * w[, b]
        for (s in seq_len(strata)) {
          total <- total + panel_table_loge(
            n1[, b, s], n0[, b, s], d_events[, b, s], x_events[, b, s], u
          )
        }
      }
      means[[k]] <- mean(total)
      row_index <- row_index + 1L
      rows[[row_index]] <- list(
        weight = name,
        eta = eta,
        mean_terminal_loge = means[[k]]
      )
    }
    chosen[[name]] <- eta_grid[[which.max(means)]]
  }
  list(chosen = chosen, curve = data.table::rbindlist(rows))
}

fixed_obf_constant <- function(...) 1.96

.first_crossing <- function(condition, horizon = ncol(condition)) {
  condition <- as.matrix(condition)
  crossed <- rowSums(condition) > 0
  first <- max.col(condition, ties.method = "first")
  first[!crossed] <- horizon
  list(crossed = crossed, first = as.integer(first))
}

summarize_panel_methods <- function(
    scenario_name,
    truth,
    tables,
    scales,
    etas,
    planned_info,
    obf_constants,
    stratified) {
  loge <- loge_paths(tables, scales, etas, stratified = stratified)
  score <- score_paths(tables, scales, stratified = stratified)
  Z <- score$Z
  reps <- dim(tables$n1)[[1L]]
  blocks <- dim(tables$n1)[[2L]]
  rows <- list()
  alpha_z <- qnorm(ALPHA)
  n_subjects <- as.integer(tables$n_subjects)
  total_events <- rowSums(matrix(tables$d, nrow = reps))
  initial_complete <- rowSums(.array_slice_matrix(tables$n1 + tables$n0, 1L))
  terminal_complete <- rowSums(.array_slice_matrix(tables$n1 + tables$n0, blocks))
  common <- list(
    mean_total_detected_events = mean(total_events),
    mean_detected_event_fraction = mean(total_events / max(n_subjects, 1L)),
    mean_initial_complete_fraction = mean(initial_complete / max(n_subjects, 1L)),
    mean_terminal_complete_fraction = mean(terminal_complete / max(n_subjects, 1L)),
    n_blocks = as.integer(blocks)
  )

  add_row <- function(method, procedure, reject, stop = NULL, terminal_metric = NULL) {
    if (!is.null(stop)) stop <- tables$block_end_times[stop]
    p <- mean(reject)
    rows[[length(rows) + 1L]] <<- c(
      list(
        scenario = scenario_name,
        truth = truth,
        method = method,
        procedure = procedure,
        rejection_probability = p,
        mc_se = sqrt(max(p * (1 - p) / reps, 0)),
        n_reps = as.integer(reps),
        median_stop_given_reject = if (!is.null(stop) && any(reject)) {
          stats::median(stop[reject])
        } else {
          NA_real_
        },
        mean_stop_capped = if (!is.null(stop)) mean(stop) else NA_real_,
        mean_terminal_metric = if (!is.null(terminal_metric)) {
          mean(terminal_metric)
        } else {
          NA_real_
        }
      ),
      common
    )
  }

  for (j in seq_along(WEIGHT_NAMES)) {
    name <- WEIGHT_NAMES[[j]]
    lp <- loge[[name]]
    crossing <- .first_crossing(lp >= LOG_BOUNDARY, blocks)
    add_row(name, "anytime_e", crossing$crossed, crossing$first, lp[, blocks])

    z <- Z[[name]]
    fixed <- z[, blocks] <= alpha_z
    add_row(name, "fixed_classical", fixed, NULL, z[, blocks])

    naive_cross <- .first_crossing(z <= alpha_z, blocks)
    add_row(
      name, "naive_repeated", naive_cross$crossed, naive_cross$first, z[, blocks]
    )

    pi <- planned_info[j, ]
    t <- if (pi[[length(pi)]] <= 0) {
      seq_len(length(pi)) / length(pi)
    } else {
      pi / pi[[length(pi)]]
    }
    boundary <- -obf_constants[[name]] / sqrt(pmax(t, 1e-8))
    obf_condition <- sweep(z, 2L, boundary, `<=`)
    obf_cross <- .first_crossing(obf_condition, blocks)
    add_row(
      name, "obrien_fleming", obf_cross$crossed, obf_cross$first, z[, blocks]
    )
  }

  paths <- c(
    stats::setNames(loge, paste0("loge_", names(loge))),
    stats::setNames(Z, paste0("z_", names(Z)))
  )
  list(summary = data.table::rbindlist(rows, fill = TRUE), paths = paths)
}

validate_identities <- function(seed = 20260826L, iterations = 200L) {
  set_panel_seed(seed)
  expectation_errors <- numeric(iterations)
  derivative_errors <- numeric(iterations)
  for (i in seq_len(iterations)) {
    n1 <- sample.int(55L, 1L) + 4L
    n0 <- sample.int(55L, 1L) + 4L
    n <- n1 + n0
    d <- sample.int(n - 1L, 1L)
    lower <- max(0L, d - n0)
    upper <- min(n1, d)
    eta <- runif(1L, -1.2, 1.2)
    weight <- runif(1L, 0.1, 3)
    u <- eta * weight
    support <- lower:upper
    p0 <- dhyper(support, m = n1, n = n0, k = d)
    factors <- exp(panel_table_loge(
      rep(n1, length(support)), rep(n0, length(support)),
      rep(d, length(support)), support, rep(u, length(support))
    ))
    expectation_errors[[i]] <- abs(sum(p0 * factors) - 1)

    x_observed <- sample(support, 1L)
    h <- 1e-6
    plus <- panel_table_loge(n1, n0, d, x_observed, h * weight)
    minus <- panel_table_loge(n1, n0, d, x_observed, -h * weight)
    numerical <- (plus - minus) / (2 * h)
    analytic <- weight * (x_observed - d * n1 / n)
    derivative_errors[[i]] <- abs(numerical - analytic)
  }
  max_expectation <- max(expectation_errors)
  max_derivative <- max(derivative_errors)
  list(
    max_hypergeom_e_expectation_error = max_expectation,
    max_score_derivative_error = max_derivative,
    all_within_tolerance = isTRUE(
      max_expectation < 1e-10 && max_derivative < 1e-7
    )
  )
}

run_interval_self_tests <- function(seed = 20260826L, smoke = TRUE) {
  .panel_core$ensure_simulation_kernels()
  degenerate <- panel_table_loge(
    c(5L, 0L, 5L, 5L), c(5L, 5L, 0L, 5L),
    c(0L, 1L, 1L, 10L), c(0L, 0L, 0L, 5L), rep(-0.3, 4L)
  )
  if (any(degenerate != 0)) stop("Degenerate table factors must be one", call. = FALSE)

  checks <- validate_identities(
    seed = derived_seed(seed, 60000L),
    iterations = if (isTRUE(smoke)) 40L else 200L
  )
  if (!isTRUE(checks$all_within_tolerance)) {
    stop("Algebraic interval-engine validation failed", call. = FALSE)
  }

  tiny <- simulate_panel_tables(
    PANEL_SCENARIOS$covariate_null,
    reps = 3L,
    seed = derived_seed(seed, 60001L),
    stratified = FALSE
  )
  expected_dim <- c(3L, length(PANEL_SCENARIOS$covariate_null$visits) - 1L, 1L)
  if (!identical(dim(tiny$n1), expected_dim)) {
    stop("Panel simulator returned an unexpected table shape", call. = FALSE)
  }
  if (any(tiny$x > tiny$d) || any(tiny$d > tiny$n1 + tiny$n0)) {
    stop("Panel table margins are internally inconsistent", call. = FALSE)
  }
  invisible(checks)
}

.standard_summary_columns <- c(
  "scenario", "truth", "method", "procedure", "rejection_probability",
  "mc_se", "n_reps", "median_stop_given_reject", "mean_stop_capped",
  "mean_terminal_metric", "mean_total_detected_events",
  "mean_detected_event_fraction", "mean_initial_complete_fraction",
  "mean_terminal_complete_fraction", "n_blocks", "resolved_fraction"
)

.require_interval_packages <- function(make_plots = TRUE) {
  required <- c("Rcpp", "data.table", "jsonlite")
  if (isTRUE(make_plots)) required <- c(required, "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(required)
}

run_covariate_single_size <- function(args) {
  start_time <- proc.time()[["elapsed"]]
  .require_interval_packages(make_plots = FALSE)
  panel_rng_kind()
  .panel_core$ensure_simulation_kernels()
  seed_log <- new_seed_registry()
  eta_grid <- seq(-1.6, -0.015, length.out = args$eta_grid)
  all_summary <- list()
  all_eta <- list()
  all_design <- list()

  log_progress <- function(...) {
    if (!isTRUE(args$quiet)) message(...)
  }

  log_progress("Covariate adjustment and delayed entry")
  cov_alt <- PANEL_SCENARIOS$covariate_alt
  cov_null <- PANEL_SCENARIOS$covariate_null
  null_pilot_adjusted <- simulate_panel_tables(
    cov_null, args$pilot_null,
    derived_seed(args$seed, 7000L, "covariate_null_pilot", seed_log),
    stratified = TRUE
  )
  alt_pilot_adjusted <- simulate_panel_tables(
    cov_alt, args$pilot_alt,
    derived_seed(args$seed, 7001L, "covariate_alternative_pilot", seed_log),
    stratified = TRUE
  )
  planning_adjusted <- planned_scales_and_information(null_pilot_adjusted)
  calibration_adjusted <- calibrate_eta_grid(
    alt_pilot_adjusted, planning_adjusted$scales, eta_grid, stratified = TRUE
  )
  calibration_adjusted$curve[, scenario := "covariate_stratified"]
  all_eta[[length(all_eta) + 1L]] <- calibration_adjusted$curve
  obf_adjusted <- stats::setNames(numeric(length(WEIGHT_NAMES)), WEIGHT_NAMES)
  for (j in seq_along(WEIGHT_NAMES)) {
    times <- planning_adjusted$planned_info[j, ] /
      max(planning_adjusted$planned_info[j, ncol(planning_adjusted$planned_info)], 1e-12)
    obf_adjusted[[j]] <- fixed_obf_constant(
      times, ALPHA,
      derived_seed(args$seed, 7100L + (j - 1L),
                   paste0("covariate_stratified_obf_", WEIGHT_NAMES[[j]]), seed_log),
      args$obf_paths
    )
  }

  null_pilot_unadjusted <- collapse_strata(null_pilot_adjusted)
  alt_pilot_unadjusted <- collapse_strata(alt_pilot_adjusted)
  planning_unadjusted <- planned_scales_and_information(null_pilot_unadjusted)
  calibration_unadjusted <- calibrate_eta_grid(
    alt_pilot_unadjusted, planning_unadjusted$scales, eta_grid, stratified = FALSE
  )
  calibration_unadjusted$curve[, scenario := "covariate_unadjusted"]
  all_eta[[length(all_eta) + 1L]] <- calibration_unadjusted$curve
  obf_unadjusted <- stats::setNames(numeric(length(WEIGHT_NAMES)), WEIGHT_NAMES)
  for (j in seq_along(WEIGHT_NAMES)) {
    times <- planning_unadjusted$planned_info[j, ] /
      max(planning_unadjusted$planned_info[j, ncol(planning_unadjusted$planned_info)], 1e-12)
    obf_unadjusted[[j]] <- fixed_obf_constant(
      times, ALPHA,
      derived_seed(args$seed, 7150L + (j - 1L),
                   paste0("covariate_unadjusted_obf_", WEIGHT_NAMES[[j]]), seed_log),
      args$obf_paths
    )
  }

  obf_adjusted[] <- 1.96
  obf_unadjusted[] <- 1.96

  null_eval_adjusted <- simulate_panel_tables(
    cov_null, args$eval_reps,
    derived_seed(args$seed, 7200L, "covariate_null_evaluation", seed_log),
    stratified = TRUE
  )
  alt_eval_adjusted <- simulate_panel_tables(
    cov_alt, args$eval_reps,
    derived_seed(args$seed, 7300L, "covariate_alternative_evaluation", seed_log),
    stratified = TRUE
  )
  null_eval_unadjusted <- collapse_strata(null_eval_adjusted)
  alt_eval_unadjusted <- collapse_strata(alt_eval_adjusted)

  covariate_runs <- list(
    list(
      "covariate_stratified", "H0", null_eval_adjusted,
      planning_adjusted, calibration_adjusted, obf_adjusted, TRUE
    ),
    list(
      "covariate_stratified", "H1", alt_eval_adjusted,
      planning_adjusted, calibration_adjusted, obf_adjusted, TRUE
    ),
    list(
      "covariate_unadjusted", "H0", null_eval_unadjusted,
      planning_unadjusted, calibration_unadjusted, obf_unadjusted, FALSE
    ),
    list(
      "covariate_unadjusted", "H1", alt_eval_unadjusted,
      planning_unadjusted, calibration_unadjusted, obf_unadjusted, FALSE
    )
  )
  covariate_paths <- list()
  for (run_spec in covariate_runs) {
    result <- summarize_panel_methods(
      run_spec[[1L]], run_spec[[2L]], run_spec[[3L]],
      run_spec[[4L]]$scales, run_spec[[5L]]$chosen,
      run_spec[[4L]]$planned_info, run_spec[[6L]], run_spec[[7L]]
    )
    all_summary[[length(all_summary) + 1L]] <- result$summary
    covariate_paths[[paste(run_spec[[1L]], run_spec[[2L]], sep = "_")]] <- result$paths
  }

  for (design_spec in list(
    list("covariate_stratified", planning_adjusted, calibration_adjusted, obf_adjusted),
    list("covariate_unadjusted", planning_unadjusted, calibration_unadjusted, obf_unadjusted)
  )) {
    design <- list(
      scenario = design_spec[[1L]],
      description = PANEL_SCENARIOS$covariate_alt$description
    )
    for (j in seq_along(WEIGHT_NAMES)) {
      name <- WEIGHT_NAMES[[j]]
      design[[paste0("scale_", name)]] <- design_spec[[2L]]$scales[[j]]
      design[[paste0("eta_", name)]] <- design_spec[[3L]]$chosen[[name]]
      design[[paste0("obf_c_", name)]] <- design_spec[[4L]][[name]]
      design[[paste0("planned_information_", name)]] <-
        as.numeric(design_spec[[2L]]$planned_info[j, ])
    }
    all_design[[length(all_design) + 1L]] <- design
  }

  summary <- data.table::rbindlist(all_summary, fill = TRUE, use.names = TRUE)
  for (column in setdiff(.standard_summary_columns, names(summary))) summary[, (column) := NA]
  data.table::setcolorder(summary, .standard_summary_columns)
  eta_df <- data.table::rbindlist(all_eta, fill = TRUE)
  data.table::setcolorder(
    eta_df, c("weight", "eta", "mean_terminal_loge", "scenario")
  )
  results_dir <- file.path(args$output_root, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  .panel_core$write_result_csv(
    summary, file.path(results_dir, "covariate_stratification.csv")
  )
  .panel_core$write_result_csv(
    eta_df, file.path(results_dir, "covariate_eta_calibration.csv")
  )
  .panel_core$write_result_json(
    all_design, file.path(results_dir, "covariate_designs.json")
  )
  characteristic_columns <- c(
    "scenario", "truth", "n_reps", "n_blocks", "mean_total_detected_events",
    "mean_detected_event_fraction", "mean_initial_complete_fraction",
    "mean_terminal_complete_fraction", "resolved_fraction"
  )
  characteristics <- unique(summary[, ..characteristic_columns])
  .panel_core$write_result_csv(
    characteristics, file.path(results_dir, "covariate_characteristics.csv")
  )

  validation_seed <- derived_seed(
    args$seed, 60000L, "algebra_validation", seed_log
  )
  validation <- validate_identities(
    validation_seed, iterations = if (isTRUE(args$smoke)) 40L else 200L
  )
  elapsed <- proc.time()[["elapsed"]] - start_time
  packages <- c("Rcpp", "data.table", "jsonlite", "ggplot2")
  r_metadata <- .panel_core$collect_r_metadata(packages)
  validation <- c(
    validation,
    list(
      r_version = r_metadata$r,
      platform = r_metadata$platform,
      os = r_metadata$os,
      rng_kind = list(
        kind = r_metadata$rng_kind[[1L]],
        normal_kind = r_metadata$rng_kind[[2L]],
        sample_kind = r_metadata$rng_kind[[3L]],
        stream_strategy = paste(
          "L'Ecuyer-CMRG stream reset for every named design-cell",
          "seed offset; the complete seed map is recorded below."
        )
      ),
      package_versions = r_metadata$packages,
      seed = as.integer(args$seed),
      derived_seeds = seed_registry_list(seed_log),
      elapsed_seconds = elapsed,
      settings = list(
        seed = as.integer(args$seed),
        pilot_null = as.integer(args$pilot_null),
        pilot_alt = as.integer(args$pilot_alt),
        eval_reps = as.integer(args$eval_reps),
        eta_grid = as.integer(args$eta_grid),
        obf_constant = 1.96,
        obf_calibration = "Fixed manuscript constant; no Brownian calibration",
        smoke = isTRUE(args$smoke),
        figures = FALSE
      ),
      session_info = r_metadata$session_info
    )
  )
  .panel_core$write_result_json(
    validation, file.path(results_dir, "covariate_validation.json")
  )

  if (nrow(summary) != 16L) {
    stop("Expected 16 covariate-summary rows; obtained ", nrow(summary), call. = FALSE)
  }
  probability_columns <- summary$rejection_probability
  if (any(!is.finite(probability_columns)) || any(probability_columns < 0 | probability_columns > 1)) {
    stop("A rejection probability is outside [0,1]", call. = FALSE)
  }
  if (!isTRUE(validation$all_within_tolerance)) {
    stop("Algebraic validation failed after the study run", call. = FALSE)
  }
  log_progress(sprintf("Completed in %.2f seconds", elapsed))
  invisible(list(
    summary = summary,
    eta_calibration = eta_df,
    panel_designs = all_design,
    scenario_characteristics = characteristics,
    validation = validation
  ))
}

# Run every study at the same total-N grid. Per-size artifacts are retained
# as checkpoints; root CSVs are the complete, explicitly keyed study.
run_covariate_study <- function(args) {
  old_scenarios <- PANEL_SCENARIOS
  on.exit({ PANEL_SCENARIOS <<- old_scenarios; SCENARIOS <<- old_scenarios }, add=TRUE)
  outputs <- vector("list", length(SAMPLE_SIZES))
  for (k in seq_along(SAMPLE_SIZES)) {
    n <- SAMPLE_SIZES[[k]]
    PANEL_SCENARIOS <<- lapply(old_scenarios, function(sc) { sc$n <- n; sc })
    SCENARIOS <<- PANEL_SCENARIOS
    cell <- args
    cell$seed <- as.integer(args$seed + k*100000L)
    cell$output_root <- file.path(args$output_root, "work", "covariate_by_n", paste0("n",n))
    message("Covariate study: total n = ",n)
    out <- run_covariate_single_size(cell)
    for (field in c("summary", "eta_calibration", "scenario_characteristics")) {
      data.table::set(out[[field]], j="n", value=rep(n,nrow(out[[field]])))
    }
    out$panel_designs <- lapply(out$panel_designs, function(x) { x$n <- n; x })
    out$validation$n <- n
    outputs[[k]] <- out
    rm(out); gc()
  }
  rd <- file.path(args$output_root,"results")
  dir.create(rd, recursive=TRUE, showWarnings=FALSE)
  mapping <- c(summary="covariate_stratification.csv", eta_calibration="covariate_eta_calibration.csv",
    scenario_characteristics="covariate_characteristics.csv")
  combined <- list()
  for (field in names(mapping)) {
    combined[[field]] <- data.table::rbindlist(lapply(outputs, `[[`, field), fill=TRUE)
    .panel_core$write_result_csv(combined[[field]], file.path(rd,mapping[[field]]))
  }
  .panel_core$write_result_json(unlist(lapply(outputs, `[[`, "panel_designs"),recursive=FALSE), file.path(rd,"covariate_designs.json"))
  metadata <- .panel_core$collect_r_metadata(c("Rcpp","data.table","jsonlite","ggplot2"))
  validation <- list(sample_sizes=SAMPLE_SIZES, n_definition="Total enrolled subjects across both treatment groups",
    all_within_tolerance=all(vapply(outputs,function(x) isTRUE(x$validation$all_within_tolerance),logical(1))),
    r_version=metadata$r, rng_kind=metadata$rng_kind, package_versions=metadata$packages,
    sample_size_runs=lapply(outputs,`[[`,"validation"))
  .panel_core$write_result_json(validation,file.path(rd,"covariate_validation.json"))
  invisible(combined)
}

interval_usage <- function() {
  paste(
    "Usage: Rscript code/run_covariate_study.R [options]",
    "",
    "Options:",
    "  --seed N              Base seed (default 20260826)",
    "  --pilot-null N        Null pilot replications (default 1000)",
    "  --pilot-alt N         Alternative pilot replications (default 1000)",
    "  --eval-reps N         Evaluation replications (default 5000)",
    "  --eta-grid N          Eta-grid size (default 50)",
    "  --obf-paths N         Deprecated compatibility option; OBF uses fixed 1.96",
    "  --output-root PATH    Project/output root (default parent of code/)",
    "  --smoke               Reduced all-scenario integration run",
    "  --skip-figures        Compatibility option; no covariate figures are generated",
    "  --self-test-only      Run core/panel algebra and shape tests only",
    "  --quiet               Suppress progress messages",
    "  --help                Show this message",
    sep = "\n"
  )
}

parse_interval_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- list(
    seed = 20260826L,
    pilot_null = PILOT_REPS,
    pilot_alt = PILOT_REPS,
    eval_reps = EVALUATION_REPS,
    eta_grid = ETA_GRID_SIZE,
    obf_paths = 100000L,
    output_root = PANEL_PROJECT_ROOT,
    smoke = FALSE,
    skip_figures = FALSE,
    self_test_only = FALSE,
    quiet = FALSE,
    help = FALSE
  )
  value_options <- c(
    seed = "seed", `pilot-null` = "pilot_null", `pilot-alt` = "pilot_alt",
    `eval-reps` = "eval_reps", `eta-grid` = "eta_grid",
    `obf-paths` = "obf_paths",
    `output-root` = "output_root"
  )
  flag_options <- c(
    smoke = "smoke", `skip-figures` = "skip_figures",
    `self-test-only` = "self_test_only", quiet = "quiet", help = "help"
  )
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!startsWith(token, "--")) stop("Unexpected argument: ", token, call. = FALSE)
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      bits <- strsplit(token, "=", fixed = TRUE)[[1L]]
      option <- bits[[1L]]
      value <- paste(bits[-1L], collapse = "=")
    } else {
      option <- token
      value <- NULL
    }
    if (option %in% names(flag_options)) {
      args[[flag_options[[option]]]] <- TRUE
    } else if (option %in% names(value_options)) {
      if (is.null(value)) {
        i <- i + 1L
        if (i > length(argv)) stop("Missing value for --", option, call. = FALSE)
        value <- argv[[i]]
      }
      key <- value_options[[option]]
      if (identical(key, "output_root")) {
        args[[key]] <- normalizePath(value, mustWork = FALSE)
      } else {
        parsed <- suppressWarnings(as.integer(value))
        if (is.na(parsed) || parsed < 1L) {
          stop("--", option, " must be a positive integer", call. = FALSE)
        }
        args[[key]] <- parsed
      }
    } else {
      stop("Unknown option --", option, call. = FALSE)
    }
    i <- i + 1L
  }
  if (isTRUE(args$smoke)) {
    args$pilot_null <- min(args$pilot_null, 12L)
    args$pilot_alt <- min(args$pilot_alt, 12L)
    args$eval_reps <- min(args$eval_reps, 20L)
    args$eta_grid <- min(args$eta_grid, 5L)
    args$obf_paths <- min(args$obf_paths, 200L)
  }
  args
}

interval_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  args <- parse_interval_args(argv)
  if (isTRUE(args$help)) {
    cat(interval_usage(), "\n")
    return(invisible(NULL))
  }
  if (isTRUE(args$self_test_only)) {
    checks <- run_interval_self_tests(args$seed, smoke = TRUE)
    if (!isTRUE(args$quiet)) {
      cat(jsonlite::toJSON(checks, auto_unbox = TRUE, pretty = TRUE), "\n")
    }
    return(invisible(checks))
  }
  run_interval_self_tests(args$seed, smoke = isTRUE(args$smoke))
  run_covariate_study(args)
}

if (sys.nframe() == 0L) {
  interval_main()
}
