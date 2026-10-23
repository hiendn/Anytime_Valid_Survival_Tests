# Shared native-R simulation engine for the exact/interval survival studies.
#
# This file is intentionally source-safe: sourcing it defines functions and
# scenario objects but neither compiles kernels nor runs simulations.  The C++
# kernels are loaded lazily on the first computational call.

.simulation_current_file <- local({
  frames <- sys.frames()
  candidates <- vapply(
    frames,
    function(frame) {
      value <- frame$ofile
      if (is.null(value) || length(value) != 1L) "" else as.character(value)
    },
    character(1)
  )
  candidates <- candidates[nzchar(candidates)]
  if (length(candidates)) {
    normalizePath(candidates[[length(candidates)]], winslash = "/", mustWork = FALSE)
  } else {
    candidate <- file.path(getwd(), "code", "simulation_core.R")
    if (!file.exists(candidate)) candidate <- file.path(getwd(), "simulation_core.R")
    normalizePath(candidate, winslash = "/", mustWork = FALSE)
  }
})

SIMULATION_CODE_DIR <- dirname(.simulation_current_file)
source(file.path(SIMULATION_CODE_DIR, "simulation_config.R"), local = TRUE)
ALPHA <- 0.05
LOG_BOUNDARY <- log(1 / ALPHA)
METHODS <- c("mixed", "all_panel", "exact_only", "all_exact_oracle")
METHOD_LABELS <- c(
  mixed = "Unified mixed-resolution",
  all_panel = "All-panel coarsening",
  exact_only = "Exact subcohort only",
  all_exact_oracle = "All-exact oracle"
)
WEIGHTS <- c("logrank", "gray_tsiatis")


new_mixed_scenario <- function(
    name,
    description,
    target_weight,
    n,
    horizon,
    base_rate,
    cure_control,
    cure_treatment,
    hazard_ratio,
    dropout_probability) {
  scenario <- list(
    name = as.character(name),
    description = as.character(description),
    target_weight = as.character(target_weight),
    n = as.integer(n),
    horizon = as.integer(horizon),
    base_rate = as.numeric(base_rate),
    cure_control = as.numeric(cure_control),
    cure_treatment = as.numeric(cure_treatment),
    hazard_ratio = as.numeric(hazard_ratio),
    dropout_probability = as.numeric(dropout_probability)
  )
  if (length(scenario$name) != 1L || !nzchar(scenario$name)) {
    stop("scenario name must be a nonempty scalar", call. = FALSE)
  }
  if (!scenario$target_weight %in% WEIGHTS) {
    stop("target_weight must be one of: ", paste(WEIGHTS, collapse = ", "), call. = FALSE)
  }
  if (is.na(scenario$n) || scenario$n <= 0L || is.na(scenario$horizon) || scenario$horizon <= 0L) {
    stop("scenario n and horizon must be positive", call. = FALSE)
  }
  rates <- c(scenario$base_rate, scenario$hazard_ratio, scenario$dropout_probability)
  if (any(!is.finite(rates)) || any(rates <= 0)) {
    stop("base_rate, hazard_ratio, and dropout_probability must be positive", call. = FALSE)
  }
  probabilities <- c(scenario$cure_control, scenario$cure_treatment, scenario$dropout_probability)
  if (any(probabilities < 0 | probabilities > 1)) {
    stop("cure and dropout probabilities must lie in [0, 1]", call. = FALSE)
  }
  structure(scenario, class = c("mixed_scenario", "list"))
}


SCENARIOS <- list(
  latency_ph = new_mixed_scenario(
    name = "latency_ph",
    description = "Equal cure fractions and susceptible proportional-hazards ratio 0.65.",
    target_weight = "logrank",
    n = REFERENCE_N,
    horizon = 8L,
    base_rate = 0.21,
    cure_control = 0.15,
    cure_treatment = 0.15,
    hazard_ratio = 0.65,
    dropout_probability = 0.14
  ),
  pure_cure = new_mixed_scenario(
    name = "pure_cure",
    description = "Equal susceptible latency distributions and cure fractions 0.15 versus 0.34.",
    target_weight = "gray_tsiatis",
    n = REFERENCE_N,
    horizon = 10L,
    base_rate = 0.25,
    cure_control = 0.15,
    cure_treatment = 0.34,
    hazard_ratio = 1.0,
    dropout_probability = 0.12
  )
)


null_version <- function(sc) {
  do.call(
    new_mixed_scenario,
    within(unclass(sc), {
      cure_treatment <- cure_control
      hazard_ratio <- 1.0
      description <- "No treatment effect under the same observation design."
    })
  )
}


configure_simulation_rng <- function() {
  RNGkind(
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  invisible(RNGkind())
}


.simulation_kernels_loaded <- FALSE

ensure_simulation_kernels <- function(rebuild = FALSE) {
  required <- c(
    "table_loge_cpp", "calibration_curve_cpp", "row_cumsums_cpp",
    "build_updates_cpp"
  )
  if (!rebuild && .simulation_kernels_loaded &&
      all(vapply(required, exists, logical(1), mode = "function", inherits = TRUE))) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("The mixed-resolution engine requires the Rcpp package", call. = FALSE)
  }
  kernel_path <- file.path(SIMULATION_CODE_DIR, "simulation_kernels.cpp")
  if (!file.exists(kernel_path)) {
    stop("Cannot find simulation kernel source: ", kernel_path, call. = FALSE)
  }
  cache_dir <- Sys.getenv(
    "SIMULATION_RCPP_CACHE_DIR",
    unset = file.path(tempdir(), "mixed_resolution_rcpp_cache")
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  old_cppflags <- Sys.getenv("PKG_CPPFLAGS", unset = "")
  on.exit(Sys.setenv(PKG_CPPFLAGS = old_cppflags), add = TRUE)
  # Some macOS Command Line Tools installations keep libc++ headers inside
  # the active SDK but omit the usual forwarding headers.  Add the SDK path
  # only when it exists; this is harmless on complete installations and keeps
  # the project runnable on the current arm64 host.
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    sdk_path <- tryCatch(
      system2("xcrun", "--show-sdk-path", stdout = TRUE, stderr = FALSE),
      error = function(error) character()
    )
    if (length(sdk_path) && nzchar(sdk_path[[1L]])) {
      libcxx_path <- file.path(sdk_path[[1L]], "usr", "include", "c++", "v1")
      if (dir.exists(libcxx_path)) {
        extra <- paste("-isystem", shQuote(libcxx_path))
        Sys.setenv(PKG_CPPFLAGS = trimws(paste(old_cppflags, extra)))
      }
    }
  }
  Rcpp::sourceCpp(
      file = kernel_path,
      rebuild = isTRUE(rebuild),
      cacheDir = cache_dir,
      showOutput = FALSE,
      verbose = FALSE
    )
  missing <- required[!vapply(required, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) {
    stop("Rcpp kernels failed to export: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  .simulation_kernels_loaded <<- TRUE
  invisible(TRUE)
}


simulate_latent <- function(sc, reps, seed) {
  reps <- as.integer(reps)
  seed <- as.integer(seed)
  if (is.na(reps) || reps <= 0L) stop("reps must be positive", call. = FALSE)
  configure_simulation_rng()
  set.seed(seed)

  total <- reps * sc$n
  group <- matrix(as.integer(runif(total) < 0.5), nrow = reps, ncol = sc$n)
  cure_probability <- ifelse(group == 1L, sc$cure_treatment, sc$cure_control)
  cured <- matrix(runif(total), nrow = reps, ncol = sc$n) < cure_probability
  rates <- sc$base_rate * ifelse(group == 1L, sc$hazard_ratio, 1.0)
  event_time <- matrix(rexp(total, rate = as.vector(rates)), nrow = reps, ncol = sc$n)
  event_time[cured] <- Inf

  # The study definition uses 1 + Geometric(p), with the geometric law
  # starting at one. R's rgeom starts at zero, hence 2 + rgeom is equivalent.
  final_visit <- matrix(
    2L + rgeom(total, prob = sc$dropout_probability),
    nrow = reps,
    ncol = sc$n
  )
  last_visit <- matrix(
    as.integer(pmin(final_visit, sc$horizon)),
    nrow = reps,
    ncol = sc$n
  )
  regime_uniform <- matrix(runif(total), nrow = reps, ncol = sc$n)
  exit_time <- pmin(event_time, last_visit)
  exit_order <- t(vapply(
    seq_len(reps),
    function(rep) order(exit_time[rep, ], seq_len(sc$n), method = "radix"),
    integer(sc$n)
  ))
  storage.mode(exit_order) <- "integer"

  structure(
    list(
      group = group,
      event_time = event_time,
      last_visit = last_visit,
      regime_uniform = regime_uniform,
      exit_order = exit_order
    ),
    class = c("latent_survival_data", "list")
  )
}


.restore_shape <- function(values, template) {
  template_dim <- dim(template)
  if (!is.null(template_dim)) {
    dim(values) <- template_dim
    dimnames(values) <- dimnames(template)
  }
  values
}


table_loge <- function(n1, n0, d, x, u, tol = 1e-15) {
  ensure_simulation_kernels()
  lengths <- lengths(list(n1, n0, d, x, u))
  if (length(unique(lengths)) != 1L) {
    stop("n1, n0, d, x, and u must have identical lengths", call. = FALSE)
  }
  values <- table_loge_cpp(
    as.integer(n1), as.integer(n0), as.integer(d), as.integer(x),
    as.numeric(u), as.numeric(tol)
  )
  .restore_shape(values, n1)
}


row_cumsums <- function(x) {
  if (is.null(dim(x))) return(cumsum(as.numeric(x)))
  ensure_simulation_kernels()
  row_cumsums_cpp(matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x)))
}


build_updates <- function(data, sc, rho, method) {
  ensure_simulation_kernels()
  method <- match.arg(method, METHODS)
  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho < 0 || rho > 1) {
    stop("rho must be a scalar in [0, 1]", call. = FALSE)
  }
  required <- c("group", "event_time", "last_visit", "regime_uniform", "exit_order")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("latent data omit: ", paste(missing, collapse = ", "), call. = FALSE)
  method_code <- match(method, METHODS) - 1L
  build_updates_cpp(
    data$group,
    data$event_time,
    data$last_visit,
    data$regime_uniform,
    data$exit_order,
    rho,
    as.integer(sc$horizon),
    as.integer(method_code),
    as.numeric(sc$cure_control),
    as.numeric(sc$base_rate)
  )
}


weights_for <- function(updates, weight_name) {
  weight_name <- match.arg(weight_name, WEIGHTS)
  if (weight_name == "logrank") {
    return(matrix(1.0, nrow = nrow(updates$time), ncol = ncol(updates$time)))
  }
  updates$gt_weight
}


hypergeom_mean_var <- function(n1, n0, d) {
  n <- n1 + n0
  mean <- matrix(0.0, nrow = nrow(n), ncol = ncol(n))
  var <- matrix(0.0, nrow = nrow(n), ncol = ncol(n))
  ok <- n > 0
  mean[ok] <- d[ok] * n1[ok] / n[ok]
  ok_var <- n > 1 & d > 0 & d < n & n1 > 0 & n0 > 0
  nn <- as.numeric(n[ok_var])
  dd <- as.numeric(d[ok_var])
  probability <- as.numeric(n1[ok_var]) / nn
  var[ok_var] <- dd * probability * (1 - probability) * (nn - dd) / (nn - 1)
  list(mean = mean, var = var)
}


terminal_loge <- function(updates, eta, weight_name, tol = 1e-15) {
  weight <- weights_for(updates, weight_name)
  increments <- table_loge(
    updates$n1, updates$n0, updates$d, updates$x,
    as.numeric(eta) * weight,
    tol = tol
  )
  rowSums(increments)
}


calibrate_eta <- function(updates, weight_name, grid, tol = 1e-15) {
  ensure_simulation_kernels()
  grid <- as.numeric(grid)
  if (!length(grid) || any(!is.finite(grid))) stop("grid must contain finite values", call. = FALSE)
  weight <- weights_for(updates, weight_name)
  values <- calibration_curve_cpp(
    as.integer(updates$n1),
    as.integer(updates$n0),
    as.integer(updates$d),
    as.integer(updates$x),
    as.numeric(weight),
    grid,
    as.integer(nrow(updates$n1)),
    as.numeric(tol)
  )
  index <- which.max(values)
  list(
    eta = grid[[index]],
    curve = data.frame(
      eta = grid,
      mean_terminal_loge = as.numeric(values),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}


evaluate_updates <- function(updates, eta, weight_name, horizon, tol = 1e-15) {
  weight <- weights_for(updates, weight_name)
  increments <- table_loge(
    updates$n1, updates$n0, updates$d, updates$x,
    as.numeric(eta) * weight,
    tol = tol
  )
  loge <- row_cumsums(increments)
  cross <- loge >= LOG_BOUNDARY
  reject <- rowSums(cross) > 0L
  first_index <- max.col(cross, ties.method = "first")
  stop_time <- rep(as.numeric(horizon), length(reject))
  rejected_rows <- which(reject)
  if (length(rejected_rows)) {
    stop_time[rejected_rows] <- updates$time[cbind(
      rejected_rows,
      first_index[rejected_rows]
    )]
  }

  moments <- hypergeom_mean_var(updates$n1, updates$n0, updates$d)
  score_increments <- weight * (updates$x - moments$mean)
  variance_increments <- weight * weight * moments$var
  score <- row_cumsums(score_increments)
  variance <- row_cumsums(variance_increments)
  z_path <- matrix(0.0, nrow = nrow(score), ncol = ncol(score))
  positive <- variance > 1e-12
  z_path[positive] <- score[positive] / sqrt(variance[positive])

  terminal_z <- numeric(length(reject))
  valid_rows <- which(updates$count > 0L)
  if (length(valid_rows)) {
    terminal_z[valid_rows] <- z_path[cbind(valid_rows, updates$count[valid_rows])]
  }
  fixed_reject <- terminal_z <= stats::qnorm(ALPHA)
  naive_reject <- rowSums(z_path <= stats::qnorm(ALPHA)) > 0L

  list(
    reject = reject,
    stop_time = stop_time,
    terminal_loge = rowSums(increments),
    terminal_z = terminal_z,
    fixed_reject = fixed_reject,
    naive_reject = naive_reject,
    loge_path = loge
  )
}


# Public alias used by downstream R entry points.
evaluate <- evaluate_updates


observation_mix <- function(data, rho, horizon = NULL) {
  exact <- data$regime_uniform < rho
  panel <- !exact
  observed_exact_event <- exact & data$event_time <= data$last_visit
  panel_left <- panel & data$event_time <= 1.0
  panel_interval <- panel & data$event_time > 1.0 & data$event_time <= data$last_visit
  right <- data$event_time > data$last_visit
  total <- length(data$group)
  list(
    exact_event_fraction = sum(observed_exact_event) / total,
    left_censored_fraction = sum(panel_left) / total,
    interval_censored_fraction = sum(panel_interval) / total,
    right_censored_fraction = sum(right) / total,
    exact_monitoring_fraction = mean(exact)
  )
}


summarize_mixed_resolution <- function(
    scenario,
    truth,
    rho,
    method,
    weight,
    eta,
    eval_result,
    updates) {
  reps <- length(eval_result$reject)
  make_row <- function(procedure, rejection, stop, metric) {
    probability <- mean(rejection)
    data.frame(
      scenario = scenario,
      truth = truth,
      exact_fraction = rho,
      method = method,
      method_label = unname(METHOD_LABELS[[method]]),
      weight = weight,
      eta = eta,
      procedure = procedure,
      rejection_probability = probability,
      mc_se = sqrt(max(probability * (1 - probability) / reps, 0)),
      n_reps = reps,
      median_stop_given_reject = if (!is.null(stop) && any(rejection)) {
        stats::median(stop[rejection])
      } else {
        NA_real_
      },
      mean_stop_capped = if (is.null(stop)) NA_real_ else mean(stop),
      mean_terminal_metric = mean(metric),
      mean_updates = mean(updates$count),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  rbind(
    make_row("anytime_e", eval_result$reject, eval_result$stop_time, eval_result$terminal_loge),
    make_row("fixed_classical", eval_result$fixed_reject, NULL, eval_result$terminal_z),
    make_row("naive_repeated", eval_result$naive_reject, NULL, eval_result$terminal_z)
  )
}


paired_difference <- function(
    scenario,
    rho,
    a_name,
    a,
    b_name,
    b,
    confidence = 0.99) {
  difference <- as.numeric(a) - as.numeric(b)
  estimate <- mean(difference)
  standard_error <- if (length(difference) > 1L) {
    stats::sd(difference) / sqrt(length(difference))
  } else {
    NA_real_
  }
  z <- stats::qnorm(0.5 + confidence / 2)
  data.frame(
    scenario = scenario,
    exact_fraction = rho,
    method_a = a_name,
    method_b = b_name,
    difference_a_minus_b = estimate,
    mc_se = standard_error,
    confidence_level = confidence,
    ci_lower = estimate - z * standard_error,
    ci_upper = estimate + z * standard_error,
    n_reps = length(difference),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


validate_algebra <- function(seed = 9817L, iterations = 100L) {
  configure_simulation_rng()
  set.seed(as.integer(seed))
  max_expectation_error <- 0
  max_derivative_error <- 0
  for (iteration in seq_len(as.integer(iterations))) {
    n1 <- sample.int(18L, 1L) + 1L
    n0 <- sample.int(18L, 1L) + 1L
    n <- n1 + n0
    events <- sample.int(n - 1L, 1L)
    eta <- runif(1L, -1, 1)
    weight <- runif(1L, 0.2, 3)
    lo <- max(0L, events - n0)
    hi <- min(n1, events)
    support <- seq.int(lo, hi)
    null_probability <- stats::dhyper(support, n1, n0, events)
    e_factor <- exp(table_loge(
      rep.int(n1, length(support)),
      rep.int(n0, length(support)),
      rep.int(events, length(support)),
      support,
      rep.int(eta * weight, length(support))
    ))
    max_expectation_error <- max(
      max_expectation_error,
      abs(sum(null_probability * e_factor) - 1)
    )

    observed <- support[[ceiling(length(support) / 2)]]
    epsilon <- 1e-6
    plus <- table_loge(n1, n0, events, observed, epsilon * weight)
    minus <- table_loge(n1, n0, events, observed, -epsilon * weight)
    numeric_derivative <- (plus - minus) / (2 * epsilon)
    analytic_derivative <- weight * (observed - events * n1 / n)
    max_derivative_error <- max(
      max_derivative_error,
      abs(numeric_derivative - analytic_derivative)
    )
  }
  list(
    max_conditional_e_expectation_error = max_expectation_error,
    max_score_derivative_error = max_derivative_error,
    algebra_within_tolerance = isTRUE(
      max_expectation_error < 1e-10 && max_derivative_error < 1e-7
    )
  )
}


simulation_core_self_test <- function(stop_on_failure = TRUE) {
  algebra <- validate_algebra(iterations = 25L)
  scenario <- new_mixed_scenario(
    name = "self_test",
    description = "Small deterministic self-test scenario.",
    target_weight = "logrank",
    n = 24L,
    horizon = 4L,
    base_rate = 0.25,
    cure_control = 0.1,
    cure_treatment = 0.1,
    hazard_ratio = 0.75,
    dropout_probability = 0.2
  )
  latent <- simulate_latent(scenario, reps = 8L, seed = 19421L)
  rho_zero_mixed <- build_updates(latent, scenario, 0, "mixed")
  rho_zero_panel <- build_updates(latent, scenario, 0, "all_panel")
  rho_one_mixed <- build_updates(latent, scenario, 1, "mixed")
  rho_one_oracle <- build_updates(latent, scenario, 1, "all_exact_oracle")
  fields <- c("n1", "n0", "d", "x", "time", "count")
  zero_identity <- all(vapply(
    fields,
    function(field) identical(rho_zero_mixed[[field]], rho_zero_panel[[field]]),
    logical(1)
  ))
  one_identity <- all(vapply(
    fields,
    function(field) identical(rho_one_mixed[[field]], rho_one_oracle[[field]]),
    logical(1)
  ))
  eta_fit <- calibrate_eta(
    rho_zero_panel,
    "logrank",
    seq(-0.5, -0.05, length.out = 4L)
  )
  evaluation <- evaluate_updates(
    rho_zero_panel,
    eta_fit$eta,
    "logrank",
    scenario$horizon
  )
  dimensions_hold <- identical(dim(evaluation$loge_path), dim(rho_zero_panel$n1))
  finite_terminal <- all(is.finite(evaluation$terminal_loge))
  passed <- isTRUE(algebra$algebra_within_tolerance) && zero_identity && one_identity &&
    dimensions_hold && finite_terminal
  result <- c(
    algebra,
    list(
      rho0_mixed_equals_panel = zero_identity,
      rho1_mixed_equals_oracle = one_identity,
      evaluation_dimensions_hold = dimensions_hold,
      finite_terminal_loge = finite_terminal,
      passed = passed
    )
  )
  if (isTRUE(stop_on_failure) && !passed) {
    stop("simulation-core self-test failed", call. = FALSE)
  }
  result
}


.method_colours <- c(
  mixed = "#1B6CA8",
  all_panel = "#E07A1F",
  exact_only = "#3A923A",
  all_exact_oracle = "#8E5AA7"
)


.open_two_panel_pdf <- function(path) {
  grDevices::pdf(path, width = 11, height = 4.2, onefile = TRUE)
  # Reserve a small device-wide lower margin for a shared legend. This keeps
  # legends out of both panels and prevents them from obscuring trajectories.
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.5, 1.2), oma = c(1.8, 0, 0, 0))
  invisible(path)
}


.draw_method_lines <- function(
    data,
    value,
    title,
    ylab,
    ylim,
    horizontal = NULL,
    show_y = TRUE) {
  graphics::par(mar = c(4.2, if (show_y) 4.4 else 1.2, 2.5, 1.2))
  graphics::plot(
    NA_real_, NA_real_,
    xlim = c(0, 1), ylim = ylim,
    xlab = "Fraction continuously monitored",
    ylab = if (show_y) ylab else "",
    yaxt = if (show_y) "s" else "n",
    main = title
  )
  graphics::grid(col = "grey88")
  if (!is.null(horizontal)) graphics::abline(h = horizontal, lty = 2, col = "grey35")
  for (method in METHODS) {
    selected <- data[data$method == method, , drop = FALSE]
    selected <- selected[order(selected$exact_fraction), , drop = FALSE]
    graphics::lines(
      selected$exact_fraction,
      selected[[value]],
      type = "o", pch = 16, lwd = 1.4,
      col = .method_colours[[method]]
    )
  }
}


.add_device_legend <- function(legend, col, lty = NULL, pch = NULL, fill = NULL) {
  # Draw in a full-device overlay rather than relative to the active right
  # panel. The reserved outer margin keeps this legend clear of x-axis labels.
  graphics::par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  arguments <- list(
    x = "bottom",
    legend = legend,
    horiz = TRUE,
    cex = 0.70,
    bty = "n",
    xpd = NA,
    inset = c(0, 0.008)
  )
  if (!is.null(col)) arguments$col <- col
  if (!is.null(lty)) arguments$lty <- lty
  if (!is.null(pch)) arguments$pch <- pch
  if (!is.null(fill)) arguments$fill <- fill
  do.call(graphics::legend, arguments)
  invisible(NULL)
}


make_mixed_resolution_figures <- function(summary, observation, example_paths, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  alternative <- summary[summary$procedure == "anytime_e" & summary$truth == "H1", , drop = FALSE]
  null <- summary[summary$procedure == "anytime_e" & summary$truth == "H0", , drop = FALSE]
  scenario_titles <- c(
    latency_ph = "Latency proportional hazards",
    pure_cure = "Pure cure-rate difference"
  )

  power_path <- file.path(figures_dir, "mixed_resolution_power.pdf")
  .open_two_panel_pdf(power_path)
  for (index in seq_along(scenario_titles)) {
    scenario <- names(scenario_titles)[[index]]
    .draw_method_lines(
      alternative[alternative$scenario == scenario, , drop = FALSE],
      "rejection_probability",
      scenario_titles[[scenario]],
      "Finite-horizon anytime power",
      c(0, 1),
      show_y = index == 1L
    )
  }
  .add_device_legend(
    unname(METHOD_LABELS[METHODS]),
    col = .method_colours[METHODS], lty = 1, pch = 16
  )
  grDevices::dev.off()

  loge_path <- file.path(figures_dir, "mixed_resolution_loge.pdf")
  .open_two_panel_pdf(loge_path)
  terminal_range <- range(c(alternative$mean_terminal_metric, LOG_BOUNDARY), finite = TRUE)
  padding <- max(diff(terminal_range) * 0.08, 0.25)
  terminal_limits <- terminal_range + c(-padding, padding)
  for (index in seq_along(scenario_titles)) {
    scenario <- names(scenario_titles)[[index]]
    .draw_method_lines(
      alternative[alternative$scenario == scenario, , drop = FALSE],
      "mean_terminal_metric",
      scenario_titles[[scenario]],
      "Mean terminal log e-value",
      terminal_limits,
      horizontal = LOG_BOUNDARY,
      show_y = index == 1L
    )
  }
  .add_device_legend(
    unname(METHOD_LABELS[METHODS]),
    col = .method_colours[METHODS], lty = 1, pch = 16
  )
  grDevices::dev.off()

  null_path <- file.path(figures_dir, "mixed_resolution_null.pdf")
  .open_two_panel_pdf(null_path)
  null_titles <- c(latency_ph = "Latency design under H0", pure_cure = "Cure design under H0")
  for (index in seq_along(null_titles)) {
    scenario <- names(null_titles)[[index]]
    .draw_method_lines(
      null[null$scenario == scenario, , drop = FALSE],
      "rejection_probability",
      null_titles[[scenario]],
      "Null anytime crossing probability",
      c(0, 0.08),
      horizontal = ALPHA,
      show_y = index == 1L
    )
  }
  .add_device_legend(
    unname(METHOD_LABELS[METHODS]),
    col = .method_colours[METHODS], lty = 1, pch = 16
  )
  grDevices::dev.off()

  observation_path <- file.path(figures_dir, "mixed_resolution_observation_mix.pdf")
  .open_two_panel_pdf(observation_path)
  components <- c(
    "exact_event_fraction", "left_censored_fraction",
    "interval_censored_fraction", "right_censored_fraction"
  )
  component_labels <- c("Exact event", "Left-censored", "Interval-censored", "Right-censored")
  component_colours <- c("#2878B5", "#F4A259", "#6FB98F", "#9B8AC4")
  short_titles <- c(latency_ph = "Latency design", pure_cure = "Cure design")
  for (index in seq_along(short_titles)) {
    scenario <- names(short_titles)[[index]]
    selected <- observation[
      observation$scenario == scenario & observation$truth == "H1",
      , drop = FALSE
    ]
    selected <- selected[order(selected$exact_fraction), , drop = FALSE]
    heights <- t(as.matrix(selected[, components, drop = FALSE]))
    show_y <- index == 1L
    graphics::par(mar = c(4.2, if (show_y) 4.4 else 1.2, 2.5, 1.2))
    graphics::barplot(
      heights,
      names.arg = format(selected$exact_fraction, trim = TRUE),
      ylim = c(0, 1), col = component_colours, border = NA,
      xlab = "Fraction continuously monitored",
      ylab = if (show_y) "Proportion of subjects" else "",
      yaxt = if (show_y) "s" else "n",
      main = short_titles[[scenario]]
    )
  }
  .add_device_legend(
    component_labels,
    col = NULL,
    fill = component_colours
  )
  grDevices::dev.off()

  examples_path <- file.path(figures_dir, "mixed_resolution_example_paths.pdf")
  .open_two_panel_pdf(examples_path)
  example_titles <- c(latency_ph = "Latency path, 50% exact", pure_cure = "Cure path, 50% exact")
  for (index in seq_along(example_titles)) {
    scenario <- names(example_titles)[[index]]
    paths <- example_paths[[scenario]]
    all_times <- unlist(lapply(paths, function(path) path$time), use.names = FALSE)
    all_loge <- unlist(lapply(paths, function(path) path$loge), use.names = FALSE)
    x_limits <- if (length(all_times)) range(c(0, all_times), finite = TRUE) else c(0, 1)
    y_limits <- range(c(0, LOG_BOUNDARY, all_loge), finite = TRUE)
    padding <- max(diff(y_limits) * 0.08, 0.25)
    graphics::par(mar = c(4.2, 4.4, 2.5, 1.2))
    graphics::plot(
      NA_real_, NA_real_, xlim = x_limits, ylim = y_limits + c(-padding, padding),
      xlab = "Calendar time", ylab = "log e-value", main = example_titles[[scenario]]
    )
    graphics::grid(col = "grey88")
    graphics::abline(h = LOG_BOUNDARY, lty = 2, col = "grey35")
    for (method in METHODS) {
      path <- paths[[method]]
      if (length(path$time)) {
        graphics::lines(path$time, path$loge, type = "s", lwd = 1.4, col = .method_colours[[method]])
      }
    }
  }
  .add_device_legend(
    unname(METHOD_LABELS[METHODS]),
    col = .method_colours[METHODS], lty = 1
  )
  grDevices::dev.off()

  invisible(c(power_path, loge_path, null_path, observation_path, examples_path))
}


write_result_csv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(data, path, na = "", quote = "auto")
  } else {
    utils::write.table(
      data, path, sep = ",", row.names = FALSE, col.names = TRUE,
      quote = TRUE, na = "", qmethod = "double"
    )
  }
  invisible(path)
}


write_result_json <- function(data, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Writing validation metadata requires jsonlite", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    data, path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = 16
  )
  invisible(path)
}


collect_r_metadata <- function(packages = c("Rcpp", "data.table", "jsonlite")) {
  versions <- vapply(
    packages,
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    },
    character(1)
  )
  list(
    r = R.version.string,
    platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    rng_kind = unname(RNGkind()),
    packages = as.list(versions),
    session_info = capture.output(utils::sessionInfo())
  )
}
