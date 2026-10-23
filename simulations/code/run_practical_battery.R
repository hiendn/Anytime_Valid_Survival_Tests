#!/usr/bin/env Rscript

# Expanded, table-oriented simulation battery for the unified survival
# e-process. Every study uses the total-cohort grid in simulation_config.R.
#
# The computational core and the two standalone study engines are sourced into
# separate module environments.  Their command-line entry points are guarded,
# so sourcing them here does not start either standalone simulation.

script_path <- function() {
  # When sourced from another Rscript, --file names the outer driver. The
  # source frame identifies this module and must take priority.
  frames <- sys.frames()
  if (length(frames)) {
    for (i in rev(seq_along(frames))) {
      candidate <- frames[[i]]$ofile
      if (!is.null(candidate) && nzchar(candidate) &&
          identical(basename(candidate), "run_practical_battery.R")) {
        return(normalizePath(candidate, mustWork = FALSE))
      }
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE), mustWork = FALSE))
  }
  stop("Cannot determine the path of run_practical_battery.R", call. = FALSE)
}

SCRIPT_FILE <- script_path()
CODE_DIR <- dirname(SCRIPT_FILE)
ROOT <- dirname(CODE_DIR)
source(file.path(CODE_DIR, "simulation_config.R"), local = TRUE)

load_r_module <- function(path, parent) {
  if (!file.exists(path)) {
    stop(sprintf("Required R module does not exist: %s", path), call. = FALSE)
  }
  module <- new.env(parent = parent)
  # source(local=) supplies a reliable `ofile` frame to nested source-safe
  # modules; sys.source() does not do so consistently when the outer program
  # itself was launched with Rscript.
  source(path, local = module, echo = FALSE, chdir = FALSE,
         keep.source = TRUE)
  module
}

# Rcpp::sourceCpp exports kernels to the process global environment by default.
# Keep that environment in the parent chain while isolating the two engines'
# scenario registries from one another.
CORE <- load_r_module(file.path(CODE_DIR, "simulation_core.R"), globalenv())
ENGINE <- load_r_module(file.path(CODE_DIR, "run_exact_interval_primary.R"), CORE)

# The exact standalone script deliberately delegates its engine API to the
# shared core rather than redeclaring wrappers.  Materialize those inherited
# bindings in the module environment so ENGINE$name has ordinary module
# semantics.
.exact_core_api <- c(
  "new_mixed_scenario", "SCENARIOS", "null_version", "simulate_latent",
  "build_updates", "weights_for", "hypergeom_mean_var", "table_loge",
  "calibrate_eta", "row_cumsums", "collect_r_metadata",
  "write_result_csv", "write_result_json"
)
for (.api_name in .exact_core_api) {
  if (!exists(.api_name, envir = ENGINE, inherits = FALSE)) {
    assign(.api_name, get(.api_name, envir = CORE, inherits = TRUE),
           envir = ENGINE)
  }
}
rm(.api_name)

ALPHA <- 0.05
LOG_BOUNDARY <- log(1 / ALPHA)
Z_FIXED <- qnorm(ALPHA)
Z_OBF_CONTINUOUS <- -qnorm(1 - ALPHA / 2)
ETA_GRID <- seq(-1.25, -0.015, length.out = ETA_GRID_SIZE)

set_cell_seed <- function(seed) {
  CORE$configure_simulation_rng()
  set.seed(as.integer(seed))
  invisible(seed)
}

row_order_matrix <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  out <- matrix(0L, nrow = nr, ncol = nc)
  for (i in seq_len(nr)) {
    out[i, ] <- order(x[i, ], method = "radix")
  }
  out
}

row_any <- function(x) {
  rowSums(x) > 0L
}

first_true <- function(x) {
  max.col(x * 1L, ties.method = "first")
}

bind_rows <- function(rows) {
  if (!length(rows)) {
    return(data.frame())
  }
  as.data.frame(data.table::rbindlist(rows, fill = TRUE, use.names = TRUE),
                stringsAsFactors = FALSE, check.names = FALSE)
}

bind_frames <- function(frames) {
  if (!length(frames)) {
    return(data.frame())
  }
  as.data.frame(data.table::rbindlist(frames, fill = TRUE, use.names = TRUE),
                stringsAsFactors = FALSE, check.names = FALSE)
}

write_result_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, na = "", quote = "auto")
}

custom_latent <- function(sc, reps, seed, treatment_probability = 0.5) {
  set_cell_seed(seed)
  reps <- as.integer(reps)
  n <- as.integer(sc$n)
  dims <- c(reps, n)

  group <- matrix(as.integer(runif(reps * n) < treatment_probability),
                  nrow = reps, ncol = n)
  cure_probability <- ifelse(group == 1L, sc$cure_treatment, sc$cure_control)
  cured <- matrix(runif(reps * n), nrow = reps, ncol = n) < cure_probability
  rates <- sc$base_rate * ifelse(group == 1L, sc$hazard_ratio, 1)
  event_time <- matrix(rexp(reps * n, rate = as.vector(rates)),
                       nrow = reps, ncol = n)
  event_time[cured] <- Inf

  # NumPy's geometric distribution has support 1,2,..., whereas R's rgeom
  # has support 0,1,....  Thus Python's 1 + Geometric(p) is 2 + rgeom(p).
  last_visit <- matrix(2L + rgeom(reps * n, prob = sc$dropout_probability),
                       nrow = reps, ncol = n)
  last_visit <- matrix(as.integer(pmin(last_visit, sc$horizon)),
                       nrow = reps, ncol = n)
  regime_uniform <- matrix(runif(reps * n), nrow = reps, ncol = n)
  exit_time <- pmin(event_time, last_visit)
  dim(exit_time) <- dims
  exit_order <- row_order_matrix(exit_time)

  list(
    group = group,
    event_time = event_time,
    last_visit = last_visit,
    regime_uniform = regime_uniform,
    exit_order = exit_order
  )
}

custom_latent_continuous_censor <- function(sc, reps, seed,
                                            treatment_probability,
                                            visit_width, total_horizon,
                                            censor_rate) {
  set_cell_seed(seed)
  reps <- as.integer(reps)
  n <- as.integer(sc$n)
  dims <- c(reps, n)

  group <- matrix(as.integer(runif(reps * n) < treatment_probability),
                  nrow = reps, ncol = n)
  cure_probability <- ifelse(group == 1L, sc$cure_treatment, sc$cure_control)
  cured <- matrix(runif(reps * n), nrow = reps, ncol = n) < cure_probability
  calendar_rate <- sc$base_rate / visit_width
  rates <- calendar_rate * ifelse(group == 1L, sc$hazard_ratio, 1)
  event_calendar <- matrix(rexp(reps * n, rate = as.vector(rates)),
                           nrow = reps, ncol = n)
  event_calendar[cured] <- Inf
  censor_calendar <- matrix(rexp(reps * n, rate = censor_rate),
                            nrow = reps, ncol = n)
  censor_calendar <- pmin(censor_calendar, total_horizon)
  dim(censor_calendar) <- dims
  last_visit <- floor(censor_calendar / visit_width + 1e-12)
  last_visit <- pmax(0, pmin(last_visit, sc$horizon))
  dim(last_visit) <- dims
  storage.mode(last_visit) <- "integer"
  event_time <- event_calendar / visit_width
  regime_uniform <- matrix(runif(reps * n), nrow = reps, ncol = n)
  exit_time <- pmin(event_time, last_visit)
  dim(exit_time) <- dims

  list(
    group = group,
    event_time = event_time,
    last_visit = last_visit,
    regime_uniform = regime_uniform,
    exit_order = row_order_matrix(exit_time)
  )
}

custom_latent_piecewise <- function(sc, reps, seed, treatment_probability,
                                    hr_early, hr_late,
                                    change_time = Inf) {
  set_cell_seed(seed)
  reps <- as.integer(reps)
  n <- as.integer(sc$n)

  group <- matrix(as.integer(runif(reps * n) < treatment_probability),
                  nrow = reps, ncol = n)
  cure_probability <- ifelse(group == 1L, sc$cure_treatment, sc$cure_control)
  cured <- matrix(runif(reps * n), nrow = reps, ncol = n) < cure_probability
  e <- matrix(rexp(reps * n, rate = 1), nrow = reps, ncol = n)
  r1 <- sc$base_rate * ifelse(group == 1L, hr_early, 1)
  r2 <- sc$base_rate * ifelse(group == 1L, hr_late, 1)

  if (is.infinite(change_time) || abs(hr_early - hr_late) < 1e-14) {
    event_time <- e / r1
  } else {
    h_change <- r1 * change_time
    before <- e <= h_change
    event_time <- matrix(0, nrow = reps, ncol = n)
    event_time[before] <- e[before] / r1[before]
    event_time[!before] <- change_time +
      (e[!before] - h_change[!before]) / r2[!before]
  }
  event_time[cured] <- Inf

  last_visit <- matrix(2L + rgeom(reps * n, prob = sc$dropout_probability),
                       nrow = reps, ncol = n)
  last_visit <- matrix(as.integer(pmin(last_visit, sc$horizon)),
                       nrow = reps, ncol = n)
  regime_uniform <- matrix(runif(reps * n), nrow = reps, ncol = n)
  exit_time <- pmin(event_time, last_visit)
  dim(exit_time) <- c(reps, n)

  list(
    group = group,
    event_time = event_time,
    last_visit = last_visit,
    regime_uniform = regime_uniform,
    exit_order = row_order_matrix(exit_time)
  )
}

null_scenario <- function(sc) {
  ENGINE$null_version(sc)
}

build <- function(data, sc, scheme) {
  rho <- if (identical(scheme, "exact")) 1 else 0
  ENGINE$build_updates(data, sc, rho, "mixed")
}

score_arrays <- function(updates, weight_name) {
  w <- ENGINE$weights_for(updates, weight_name)
  moments <- ENGINE$hypergeom_mean_var(updates$n1, updates$n0, updates$d)
  u_inc <- w * (updates$x - moments$mean)
  v_inc <- w * w * moments$var
  U <- CORE$row_cumsums(u_inc)
  V <- CORE$row_cumsums(v_inc)
  Z <- matrix(0, nrow = nrow(U), ncol = ncol(U))
  ok <- V > 1e-12
  Z[ok] <- U[ok] / sqrt(V[ok])
  list(U = U, V = V, Z = Z)
}

# Cumulative information must be aligned by prespecified calendar look, not by
# the packed position of an observed event. Zero-event visits are retained as
# looks by carrying the cumulative score forward. This applies to both schemes.
calendar_scores <- function(updates, score, look_times) {
  nr <- nrow(score$U)
  U <- V <- matrix(0, nrow = nr, ncol = length(look_times))
  for (i in seq_len(nr)) {
    k <- as.integer(updates$count[[i]])
    if (k > 0L) {
      index <- findInterval(look_times, updates$time[i, seq_len(k)])
      positive <- index > 0L
      U[i, positive] <- score$U[i, index[positive]]
      V[i, positive] <- score$V[i, index[positive]]
    }
  }
  Z <- matrix(0, nrow = nr, ncol = length(look_times))
  positive <- V > 1e-12
  Z[positive] <- U[positive] / sqrt(V[positive])
  list(U = U, V = V, Z = Z)
}

planned_information_fraction <- function(null_updates, weight_name,
                                         look_times) {
  score <- calendar_scores(null_updates, score_arrays(null_updates, weight_name),
                           look_times)
  mean_v <- colMeans(score$V)
  final <- max(tail(mean_v, 1L), 1e-12)
  list(fraction = cummax(pmin(1, pmax(0, mean_v / final))),
       look_times = look_times)
}

evaluate_all <- function(updates, eta, weight_name, horizon_calendar,
                         time_scale, planned_t) {
  w <- ENGINE$weights_for(updates, weight_name)
  inc <- ENGINE$table_loge(updates$n1, updates$n0, updates$d, updates$x,
                           eta * w)
  loge <- CORE$row_cumsums(inc)
  e_cross_mat <- loge >= LOG_BOUNDARY
  e_reject <- row_any(e_cross_mat)
  e_first <- first_true(e_cross_mat)

  score <- score_arrays(updates, weight_name)
  V <- score$V
  Z <- score$Z
  counts <- as.integer(updates$count)
  idx <- seq_along(counts)
  last_index <- pmax(counts, 1L)
  terminal_z <- Z[cbind(idx, last_index)]
  terminal_z[counts <= 0L] <- 0
  fixed_reject <- terminal_z <= Z_FIXED

  positions <- matrix(rep(seq_len(ncol(Z)), each = nrow(Z)),
                      nrow = nrow(Z), ncol = ncol(Z))
  valid <- positions <= counts
  naive_mat <- (Z <= Z_FIXED) & valid & (V > 1e-12)
  naive_reject <- row_any(naive_mat)
  naive_first <- first_true(naive_mat)

  planned_score <- calendar_scores(updates, score, planned_t$look_times)
  t <- planned_t$fraction
  obf_boundary <- Z_OBF_CONTINUOUS / sqrt(pmax(t, 1e-10))
  obf_mat <- sweep(planned_score$Z, 2L, obf_boundary, "<=") &
    (planned_score$V > 1e-12)
  obf_reject <- row_any(obf_mat)
  obf_first <- first_true(obf_mat)
  obf_stop <- rep(as.numeric(horizon_calendar), length(obf_reject))
  obf_stop[obf_reject] <- planned_t$look_times[obf_first[obf_reject]] * time_scale

  update_times <- updates$time * time_scale
  stopping_time <- function(reject, first) {
    out <- rep(as.numeric(horizon_calendar), length(reject))
    rr <- which(reject)
    if (length(rr)) {
      out[rr] <- update_times[cbind(rr, first[rr])]
    }
    out
  }

  list(
    e_reject = e_reject,
    e_stop = stopping_time(e_reject, e_first),
    fixed_reject = fixed_reject,
    obf_reject = obf_reject,
    obf_stop = obf_stop,
    naive_reject = naive_reject,
    naive_stop = stopping_time(naive_reject, naive_first),
    terminal_loge = rowSums(inc),
    maximum_loge = pmax(0, matrixStats::rowMaxs(loge)),
    terminal_z = terminal_z,
    total_events = as.numeric(rowSums(updates$d)),
    n_updates = as.numeric(counts)
  )
}

p_se <- function(x) {
  p <- mean(x)
  list(p = as.numeric(p),
       se = sqrt(max(p * (1 - p) / length(x), 0)))
}

summarize_eval <- function(study, scenario, scheme, truth, settings, ev,
                           data, horizon_calendar) {
  rows <- list()
  right_censored <- data$event_time > data$last_visit
  rc <- mean(right_censored)
  methods <- list(
    c("Anytime e-process", "e_reject", "e_stop"),
    c("Single terminal score test", "fixed_reject", NA_character_),
    c("O'Brien--Fleming-type monitoring", "obf_reject", "obf_stop"),
    c("Terminal critical value repeated", "naive_reject", "naive_stop")
  )

  for (spec in methods) {
    label <- spec[[1L]]
    key <- spec[[2L]]
    stop_key <- spec[[3L]]
    reject <- as.logical(ev[[key]])
    ps <- p_se(reject)
    stop_time <- if (!is.na(stop_key)) as.numeric(ev[[stop_key]]) else NULL
    row <- list(
      study = study,
      scenario = scenario,
      scheme = scheme,
      truth = truth,
      procedure = label,
      rejection_probability = ps$p,
      mc_se = ps$se,
      median_stop_given_reject = if (!is.null(stop_time) && any(reject))
        median(stop_time[reject]) else NA_real_,
      mean_stop_capped = if (!is.null(stop_time)) mean(stop_time) else NA_real_,
      mean_terminal_loge = if (identical(label, "Anytime e-process"))
        mean(ev$terminal_loge) else NA_real_,
      mean_maximum_loge = if (identical(label, "Anytime e-process"))
        mean(ev$maximum_loge) else NA_real_,
      mean_total_events = mean(ev$total_events),
      mean_updates = mean(ev$n_updates),
      right_censored_fraction = rc,
      horizon = as.numeric(horizon_calendar),
      n_reps = length(reject)
    )
    rows[[length(rows) + 1L]] <- c(row, settings)
  }
  rows
}

calibrate_design <- function(sc, scheme, weight, treatment_probability, seed,
                             pilot_reps) {
  alt <- custom_latent(sc, pilot_reps, seed, treatment_probability)
  up_alt <- build(alt, sc, scheme)
  calibration <- ENGINE$calibrate_eta(up_alt, weight, ETA_GRID)
  null_sc <- null_scenario(sc)
  null <- custom_latent(null_sc, pilot_reps, seed + 1L,
                        treatment_probability)
  up_null <- build(null, null_sc, scheme)
  list(
    eta = calibration$eta,
    planned_t = planned_information_fraction(up_null, weight, seq_len(sc$horizon)),
    curve = calibration$curve
  )
}

calibrate_design_from_latent <- function(sc, scheme, weight, alt, null) {
  calibration <- ENGINE$calibrate_eta(build(alt, sc, scheme), weight, ETA_GRID)
  list(
    eta = calibration$eta,
    planned_t = planned_information_fraction(
      build(null, null_scenario(sc), scheme), weight, seq_len(sc$horizon)
    ),
    curve = calibration$curve
  )
}

run_primary <- function(args, seed_base) {
  rows <- list()
  calibration_rows <- list()
  scenarios <- list(
    "Latency proportional hazards" = ENGINE$SCENARIOS[["latency_ph"]],
    "Pure cure-rate difference" = ENGINE$SCENARIOS[["pure_cure"]]
  )

  original_seed <- seed_base
  for (n in SAMPLE_SIZES) {
    seed_base <- original_seed + match(n, SAMPLE_SIZES) * 1000000L
    message(sprintf("[primary] total N=%d", n))
    for (si0 in seq_along(scenarios) - 1L) {
      name <- names(scenarios)[[si0 + 1L]]
      sc <- scenarios[[si0 + 1L]]
      sc$n <- n
      null_sc <- null_scenario(sc)

      # Independent alternative and null pilots are generated for every
      # N/scenario. Both observation maps use the same latent subjects.
      alt_pilot <- custom_latent(sc, args$pilot, seed_base + 100L * si0, 0.5)
      null_pilot <- custom_latent(null_sc, args$pilot,
                                  seed_base + 100L * si0 + 1L, 0.5)
      designs <- list()
      for (scheme in c("interval", "exact")) {
        design <- calibrate_design_from_latent(
          sc, scheme, sc$target_weight, alt_pilot, null_pilot
        )
        designs[[scheme]] <- design
        curve <- design$curve
        curve$n <- n
        curve$study <- "primary"
        curve$scenario <- name
        curve$scheme <- scheme
        calibration_rows[[length(calibration_rows) + 1L]] <- curve
      }

      truths <- list(Null = null_sc, Alternative = sc)
      for (ti0 in seq_along(truths) - 1L) {
        truth <- names(truths)[[ti0 + 1L]]
        sc_data <- truths[[ti0 + 1L]]
        # Paired evaluation applies both observation maps to one latent dataset.
        data <- custom_latent(
          sc_data, args$eval_primary,
          seed_base + 1000L + 100L * si0 + 10L * ti0, 0.5
        )
        for (scheme in c("interval", "exact")) {
          design <- designs[[scheme]]
          updates <- build(data, sc_data, scheme)
          ev <- evaluate_all(
            updates, design$eta, sc$target_weight, as.numeric(sc$horizon),
            1, design$planned_t
          )
          rows <- c(rows, summarize_eval(
            "primary", name, scheme, truth,
            list(
              n = sc$n,
              allocation_treatment = 0.5,
              dropout_probability = sc$dropout_probability,
              visit_width = 1,
              target_weight = sc$target_weight,
              eta = design$eta
            ),
            ev, data, as.numeric(sc$horizon)
          ))
        }
      }
    }
  }
  list(
    data = bind_rows(rows),
    calibration = bind_frames(calibration_rows),
    exact_interval_common_latent = TRUE
  )
}

run_design_grid <- function(args, seed_base, base_eta) {
  rows <- list()
  original_seed <- seed_base
  for (n in SAMPLE_SIZES) {
    seed_base <- original_seed + match(n, SAMPLE_SIZES) * 1000000L
    message(sprintf("[design sensitivity] total N=%d", n))
    for (dropout in c(0.08, 0.20)) {
      for (alloc in c(0.50, 1 / 3)) {
        sc <- ENGINE$new_mixed_scenario(
          name = "design_ph",
          description = "PH design sensitivity",
          target_weight = "logrank",
          n = n,
          horizon = 8L,
          base_rate = 0.21,
          cure_control = 0.15,
          cure_treatment = 0.15,
          hazard_ratio = 0.70,
          dropout_probability = dropout
        )
        null_sc <- null_scenario(sc)
        for (scheme in c("interval", "exact")) {
          scheme_offset <- 0L # same latent draws under both observation schemes
          null_pilot <- custom_latent(
            null_sc, args$pilot_small,
            seed_base + n + as.integer(dropout * 1000) +
              as.integer(alloc * 100) + scheme_offset,
            alloc
          )
          planned_t <- planned_information_fraction(
            build(null_pilot, null_sc, scheme), "logrank", seq_len(sc$horizon)
          )
          eta <- base_eta[[paste("ph", scheme, sep = "_")]]
          truths <- list(Null = null_sc, Alternative = sc)
          for (ti0 in seq_along(truths) - 1L) {
            truth <- names(truths)[[ti0 + 1L]]
            sc_data <- truths[[ti0 + 1L]]
            seed <- seed_base + 100000L + n * 17L +
              as.integer(dropout * 10000) * 3L +
              as.integer(alloc * 1000) * 5L + ti0 * 11L + scheme_offset
            data <- custom_latent(sc_data, args$eval_grid, seed, alloc)
            ev <- evaluate_all(
              build(data, sc_data, scheme), eta, "logrank", 8, 1,
              planned_t
            )
            rows <- c(rows, summarize_eval(
              "design_grid", "Latency PH, HR=0.70", scheme, truth,
              list(
                n = n,
                allocation_treatment = alloc,
                dropout_probability = dropout,
                visit_width = 1,
                target_weight = "logrank",
                eta = eta
              ),
              ev, data, 8
            ))
          }
        }
      }
    }
  }
  bind_rows(rows)
}

run_visit_frequency <- function(args, seed_base, eta_interval) {
  rows <- list()
  total_horizon <- 8
  annual_dropout_hazard <- -log(1 - 0.14)
  original_seed <- seed_base
  for (n in SAMPLE_SIZES) {
    seed_base <- original_seed + match(n, SAMPLE_SIZES) * 1000000L
    message(sprintf("[visit frequency] total N=%d", n))
    for (delta in c(0.5, 1, 2)) {
      horizon_units <- as.integer(round(total_horizon / delta))
      dropout_p <- 1 - exp(-annual_dropout_hazard * delta)
      sc <- ENGINE$new_mixed_scenario(
        name = paste0("visit_", delta),
        description = "Visit-frequency sensitivity",
        target_weight = "logrank",
        n = n,
        horizon = horizon_units,
        base_rate = 0.21 * delta,
        cure_control = 0.15,
        cure_treatment = 0.15,
        hazard_ratio = 0.70,
        dropout_probability = dropout_p
      )
      null_sc <- null_scenario(sc)
      censor_rate <- 0.085
      null_pilot <- custom_latent_continuous_censor(
        null_sc, args$pilot_small, seed_base + as.integer(delta * 100),
        0.5, delta, total_horizon, censor_rate
      )
      planned_t <- planned_information_fraction(
        build(null_pilot, null_sc, "interval"), "logrank", seq_len(sc$horizon)
      )
      truths <- list(Null = null_sc, Alternative = sc)
      for (ti0 in seq_along(truths) - 1L) {
        truth <- names(truths)[[ti0 + 1L]]
        sc_data <- truths[[ti0 + 1L]]
        data <- custom_latent_continuous_censor(
          sc_data, args$eval_visit,
          seed_base + 1000L + as.integer(delta * 100) + ti0,
          0.5, delta, total_horizon, censor_rate
        )
        ev <- evaluate_all(
          build(data, sc_data, "interval"), eta_interval, "logrank",
          total_horizon, delta, planned_t
        )
        rows <- c(rows, summarize_eval(
          "visit_frequency", "Latency PH, HR=0.70", "interval", truth,
          list(
            n = n,
            allocation_treatment = 0.5,
            dropout_probability = dropout_p,
            visit_width = delta,
            n_assessments = horizon_units,
            target_weight = "logrank",
            eta = eta_interval
          ),
          ev, data, total_horizon
        ))
      }
    }
  }
  bind_rows(rows)
}

run_effect_and_tilt <- function(args, seed_base, baseline_eta) {
  rows <- list()
  calibration_rows <- list()
  specifications <- list()
  for (hr in c(0.80, 0.70, 0.60)) {
    specifications[[length(specifications) + 1L]] <- list(
      family = "Latency PH",
      effect = hr,
      sc = ENGINE$new_mixed_scenario(
        name = paste0("ph_", hr),
        description = "Effect-size sensitivity",
        target_weight = "logrank",
        n = REFERENCE_N,
        horizon = 8L,
        base_rate = 0.21,
        cure_control = 0.15,
        cure_treatment = 0.15,
        hazard_ratio = hr,
        dropout_probability = 0.14
      ),
      weight = "logrank",
      eta_family = "ph"
    )
  }
  for (cure_t in c(0.25, 0.34, 0.43)) {
    specifications[[length(specifications) + 1L]] <- list(
      family = "Pure cure-rate",
      effect = cure_t,
      sc = ENGINE$new_mixed_scenario(
        name = paste0("cure_", cure_t),
        description = "Effect-size sensitivity",
        target_weight = "gray_tsiatis",
        n = REFERENCE_N,
        horizon = 10L,
        base_rate = 0.25,
        cure_control = 0.15,
        cure_treatment = cure_t,
        hazard_ratio = 1,
        dropout_probability = 0.12
      ),
      weight = "gray_tsiatis",
      eta_family = "cure"
    )
  }

  original_seed <- seed_base
  for (n in SAMPLE_SIZES) {
    seed_base <- original_seed + match(n, SAMPLE_SIZES) * 1000000L
    message(sprintf("[effect and tilt] total N=%d", n))
    for (si0 in seq_along(specifications) - 1L) {
      spec <- specifications[[si0 + 1L]]
      sc <- spec$sc
      sc$n <- n
      for (scheme in c("interval", "exact")) {
        scheme_offset <- 0L # same latent draws under both observation schemes
        pilot <- custom_latent(
          sc, args$pilot_small, seed_base + 100L * si0 + scheme_offset, 0.5
        )
        calibration <- ENGINE$calibrate_eta(
          build(pilot, sc, scheme), spec$weight, ETA_GRID
        )
        curve <- calibration$curve
        curve$n <- n
        curve$study <- "effect_tilt"
        curve$scenario_family <- spec$family
        curve$effect <- spec$effect
        curve$scheme <- scheme
        calibration_rows[[length(calibration_rows) + 1L]] <- curve

        data <- custom_latent(
          sc, args$eval_effect,
          seed_base + 5000L + 100L * si0 + scheme_offset, 0.5
        )
        updates <- build(data, sc, scheme)
        null_sc <- null_scenario(sc)
        null_pilot <- custom_latent(
          null_sc, args$pilot_small,
          seed_base + 8000L + 100L * si0 + scheme_offset, 0.5
        )
        planned_t <- planned_information_fraction(
          build(null_pilot, null_sc, scheme), spec$weight, seq_len(sc$horizon)
        )
        eta_choices <- list(
          "Middle-effect design tilt" = baseline_eta[[paste(spec$eta_family, scheme, sep = "_")]],
          "Effect-specific oracle tilt" = calibration$eta
        )
        for (eta_label in names(eta_choices)) {
          eta <- eta_choices[[eta_label]]
          ev <- evaluate_all(
            updates, eta, spec$weight, as.numeric(sc$horizon), 1, planned_t
          )
          ps <- p_se(ev$e_reject)
          reject <- ev$e_reject
          rows[[length(rows) + 1L]] <- list(
            n = n,
            scenario_family = spec$family,
            effect = spec$effect,
            scheme = scheme,
            tilt_choice = eta_label,
            eta = eta,
            power = ps$p,
            mc_se = ps$se,
            median_stop_given_reject = if (any(reject))
              median(ev$e_stop[reject]) else NA_real_,
            mean_stop_capped = mean(ev$e_stop),
            mean_terminal_loge = mean(ev$terminal_loge),
            mean_maximum_loge = mean(ev$maximum_loge),
            right_censored_fraction = mean(data$event_time > data$last_visit),
            n_reps = length(reject)
          )
        }
      }
    }
  }
  list(data = bind_rows(rows), calibration = bind_frames(calibration_rows))
}

all_weight_matrices <- function(updates) {
  n <- updates$n1 + updates$n0
  storage.mode(n) <- "double"
  d <- updates$d
  storage.mode(d) <- "double"
  reps <- nrow(n)
  m <- ncol(n)
  logrank <- matrix(1, nrow = reps, ncol = m)
  initial <- pmax(n[, 1L], 1)
  gehan <- sweep(n, 1L, initial, "/")
  peto <- matrix(0, nrow = reps, ncol = m)
  fh_late <- matrix(0, nrow = reps, ncol = m)
  gray <- matrix(0, nrow = reps, ncol = m)
  surv <- rep(1, reps)
  for (j in seq_len(m)) {
    peto[, j] <- surv
    fh_late[, j] <- pmax(1 - surv, 0)
    gray[, j] <- 1 / pmax(surv, 0.12)
    frac <- rep(0, reps)
    ok <- n[, j] > 0
    frac[ok] <- d[ok, j] / n[ok, j]
    surv <- surv * pmax(1 - frac, 1e-10)
  }
  list(
    "Logrank" = logrank,
    "Gehan--Breslow" = gehan,
    "Peto--Peto" = peto,
    "Fleming--Harrington (0,1)" = fh_late,
    "Gray--Tsiatis" = gray
  )
}

scaled_weights <- function(null_updates) {
  mats <- all_weight_matrices(null_updates)
  moments <- ENGINE$hypergeom_mean_var(
    null_updates$n1, null_updates$n0, null_updates$d
  )
  base <- max(mean(rowSums(moments$var)), 1e-12)
  scales <- vapply(names(mats), function(name) {
    info <- mean(rowSums(mats[[name]] * mats[[name]] * moments$var))
    sqrt(max(info / base, 1e-12))
  }, numeric(1L))
  list(matrices = mats, scales = scales)
}

calibrate_weight_eta <- function(alt_updates, scales) {
  CORE$ensure_simulation_kernels()
  mats <- all_weight_matrices(alt_updates)
  curves <- list()
  chosen <- vapply(names(mats), function(name) {
    ww <- mats[[name]] / scales[[name]]
    values <- calibration_curve_cpp(
      as.integer(alt_updates$n1), as.integer(alt_updates$n0),
      as.integer(alt_updates$d), as.integer(alt_updates$x),
      as.numeric(ww), ETA_GRID, nrow(alt_updates$n1), 1e-15
    )
    curves[[name]] <<- data.frame(weight = name, scale = scales[[name]],
                                 eta = ETA_GRID, mean_terminal_loge = values)
    ETA_GRID[[which.max(values)]]
  }, numeric(1L))
  list(eta = chosen, curve = bind_frames(curves))
}

evaluate_weight_e <- function(updates, name, scale, eta, horizon_calendar,
                              time_scale = 1, weight_matrix = NULL) {
  if (is.null(weight_matrix)) weight_matrix <- all_weight_matrices(updates)[[name]]
  w <- weight_matrix / scale
  inc <- ENGINE$table_loge(
    updates$n1, updates$n0, updates$d, updates$x, eta * w
  )
  loge <- CORE$row_cumsums(inc)
  cross <- loge >= LOG_BOUNDARY
  reject <- row_any(cross)
  first <- first_true(cross)
  stop_time <- rep(as.numeric(horizon_calendar), length(reject))
  rr <- which(reject)
  if (length(rr)) {
    stop_time[rr] <- updates$time[cbind(rr, first[rr])] * time_scale
  }
  list(
    reject = reject,
    stop = stop_time,
    terminal = rowSums(inc),
    path = loge
  )
}

logaddexp_pair <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  both_neg_inf <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
  out[both_neg_inf] <- -Inf
  out
}

run_two_scheme_weight_battery <- function(args, seed_base) {
  rows <- list()
  calibration_rows <- list()
  scenario_specs <- list(
    list("Proportional hazards", REFERENCE_N, 8L, 1, 0.21, 0.15, 0.15,
         0.65, 0.65, Inf, 0.14),
    list("Pure cure-rate", REFERENCE_N, 10L, 1, 0.25, 0.15, 0.34,
         1, 1, Inf, 0.12),
    list("Early effect", REFERENCE_N, 12L, 0.5, 0.20 * 0.5, 0, 0,
         0.65, 1, 2.5 / 0.5, 1 - exp(-0.07 * 0.5)),
    list("Delayed effect", REFERENCE_N, 16L, 0.5, 0.18 * 0.5, 0, 0,
         1, 0.55, 3 / 0.5, 1 - exp(-0.06 * 0.5)),
    list("Joint cure--latency", REFERENCE_N, 10L, 1, 0.23, 0.20, 0.30,
         0.82, 0.82, Inf, 0.035)
  )

  original_seed <- seed_base
  for (n in SAMPLE_SIZES) {
    seed_base <- original_seed + match(n, SAMPLE_SIZES) * 1000000L
    message(sprintf("[two-scheme weights] total N=%d", n))
    for (si0 in seq_along(scenario_specs) - 1L) {
      spec <- scenario_specs[[si0 + 1L]]
      label <- spec[[1L]]
      horizon <- spec[[3L]]
      time_scale <- spec[[4L]]
      base_rate <- spec[[5L]]
      c0 <- spec[[6L]]
      c1 <- spec[[7L]]
      hr1 <- spec[[8L]]
      hr2 <- spec[[9L]]
      change <- spec[[10L]]
      dropout <- spec[[11L]]
      sc <- ENGINE$new_mixed_scenario(
        name = paste0("weights_", si0),
        description = label,
        target_weight = "logrank",
        n = n,
        horizon = horizon,
        base_rate = base_rate,
        cure_control = c0,
        cure_treatment = c1,
        hazard_ratio = 1,
        dropout_probability = dropout
      )
      sc_null <- null_scenario(sc)

      for (scheme in c("interval", "exact")) {
        scheme_offset <- 0L # same latent draws under both observation schemes
        pilot_null_data <- custom_latent_piecewise(
          sc_null, args$pilot_small,
          seed_base + si0 * 100L + scheme_offset,
          0.5, 1, 1, change
        )
        pilot_alt_data <- custom_latent_piecewise(
          sc, args$pilot_small,
          seed_base + 5000L + si0 * 100L + scheme_offset,
          0.5, hr1, hr2, change
        )
        pilot_null <- build(pilot_null_data, sc_null, scheme)
        pilot_alt <- build(pilot_alt_data, sc, scheme)
        scaled <- scaled_weights(pilot_null)
        scales <- scaled$scales
        calibrated <- calibrate_weight_eta(pilot_alt, scales)
        etas <- calibrated$eta
        curve <- calibrated$curve
        curve$n <- n
        curve$scenario <- label
        curve$scheme <- scheme
        calibration_rows[[length(calibration_rows) + 1L]] <- curve

        truths <- list(
          Null = list(sc = sc_null, e1 = 1, e2 = 1, offset = 0L),
          Alternative = list(sc = sc, e1 = hr1, e2 = hr2, offset = 1L)
        )
        for (truth in names(truths)) {
          truth_spec <- truths[[truth]]
          data <- custom_latent_piecewise(
            truth_spec$sc, args$eval_weight,
            seed_base + 10000L + si0 * 100L + truth_spec$offset * 10L +
              scheme_offset,
            0.5, truth_spec$e1, truth_spec$e2, change
          )
          updates <- build(data, truth_spec$sc, scheme)
          paths <- list()
          weight_matrices <- all_weight_matrices(updates)
          method_names <- c(
            "Logrank", "Gehan--Breslow", "Peto--Peto",
            "Fleming--Harrington (0,1)", "Gray--Tsiatis"
          )
          for (method_name in method_names) {
            evaluated <- evaluate_weight_e(
              updates, method_name, scales[[method_name]], etas[[method_name]],
              horizon * time_scale, time_scale, weight_matrices[[method_name]]
            )
            if (method_name %in% c("Logrank", "Gray--Tsiatis"))
              paths[[method_name]] <- evaluated$path
            ps <- p_se(evaluated$reject)
            rows[[length(rows) + 1L]] <- list(
              scenario = label,
              scheme = scheme,
              truth = truth,
              weight = method_name,
              eta = etas[[method_name]],
              rejection_probability = ps$p,
              mc_se = ps$se,
              median_stop_given_reject = if (any(evaluated$reject))
                median(evaluated$stop[evaluated$reject]) else NA_real_,
              mean_stop_capped = mean(evaluated$stop),
              mean_terminal_loge = mean(evaluated$terminal),
              mean_maximum_loge = mean(pmax(0, matrixStats::rowMaxs(evaluated$path))),
              n_reps = length(evaluated$reject),
              n = n,
              horizon = horizon * time_scale
            )
          }

          combo <- logaddexp_pair(
            paths[["Logrank"]] + log(0.5),
            paths[["Gray--Tsiatis"]] + log(0.5)
          )
          combo_cross <- combo >= LOG_BOUNDARY
          reject <- row_any(combo_cross)
          first <- first_true(combo_cross)
          stop_time <- rep(horizon * time_scale, length(reject))
          rr <- which(reject)
          if (length(rr)) {
            stop_time[rr] <- updates$time[cbind(rr, first[rr])] * time_scale
          }
          ps <- p_se(reject)
          rows[[length(rows) + 1L]] <- list(
            scenario = label,
            scheme = scheme,
            truth = truth,
            weight = "Equal logrank--Gray--Tsiatis average",
            eta = NA_real_,
            rejection_probability = ps$p,
            mc_se = ps$se,
            median_stop_given_reject = if (any(reject))
              median(stop_time[reject]) else NA_real_,
            mean_stop_capped = mean(stop_time),
            mean_terminal_loge = mean(combo[, ncol(combo)]),
            mean_maximum_loge = mean(pmax(0, matrixStats::rowMaxs(combo))),
            n_reps = length(reject),
            n = n,
            horizon = horizon * time_scale
          )
        }
      }
    }
  }
  list(data = bind_rows(rows), calibration = bind_frames(calibration_rows))
}

null_check <- function(x, alpha = ALPHA, z = qnorm(0.975)) {
  if (!nrow(x)) {
    return(list(
      all_estimates_at_or_below_alpha = NA,
      all_estimates_within_mc_tolerance = NA,
      max_estimate = NA_real_,
      max_standardized_excess = NA_real_
    ))
  }
  p <- x$rejection_probability
  se <- x$mc_se
  standardized <- ifelse(
    se > 0,
    pmax(p - alpha, 0) / se,
    ifelse(p <= alpha, 0, Inf)
  )
  list(
    all_estimates_at_or_below_alpha = all(p <= alpha),
    all_estimates_within_mc_tolerance = all(p <= alpha + z * se),
    max_estimate = max(p),
    max_standardized_excess = max(standardized)
  )
}

validate_outputs <- function(primary, design, visit, effect) {
  primary_null <- primary[
    primary$truth == "Null" & primary$procedure == "Anytime e-process",
    , drop = FALSE
  ]
  design_null <- design[
    design$truth == "Null" & design$procedure == "Anytime e-process",
    , drop = FALSE
  ]
  visit_null <- visit[
    visit$truth == "Null" & visit$procedure == "Anytime e-process",
    , drop = FALSE
  ]
  pcheck <- null_check(primary_null)
  dcheck <- null_check(design_null)
  vcheck <- null_check(visit_null)
  all_probabilities <- c(
    primary$rejection_probability,
    design$rejection_probability,
    visit$rejection_probability,
    effect$power
  )
  list(
    # Retained names provide direct continuity with the Python validation file.
    primary_null_e_below_0_05 = pcheck$all_estimates_at_or_below_alpha,
    design_null_e_below_0_05 = dcheck$all_estimates_at_or_below_alpha,
    visit_null_e_below_0_05 = vcheck$all_estimates_at_or_below_alpha,
    # These are the preferred descriptive, Monte-Carlo-aware checks.
    primary_null_e_within_1_96_mc_se = pcheck$all_estimates_within_mc_tolerance,
    design_null_e_within_1_96_mc_se = dcheck$all_estimates_within_mc_tolerance,
    visit_null_e_within_1_96_mc_se = vcheck$all_estimates_within_mc_tolerance,
    primary_null_diagnostics = pcheck,
    design_null_diagnostics = dcheck,
    visit_null_diagnostics = vcheck,
    effect_power_in_unit_interval = all(effect$power >= 0 & effect$power <= 1),
    all_probabilities_in_unit_interval = all(
      all_probabilities >= 0 & all_probabilities <= 1
    ),
    null_checks_are_descriptive = TRUE,
    null_check_interpretation = paste(
      "Raw <= 0.05 flags are retained for comparability; preferred flags allow",
      "1.96 Monte Carlo standard errors above 0.05 and are not algebraic",
      "validity assertions."
    )
  )
}

validate_smoke_outputs <- function(primary, primary_cal, design, visit,
                                   effect, effect_cal, weights) {
  expected_rows <- c(
    expanded_primary = 32L * length(SAMPLE_SIZES),
    expanded_primary_eta_calibration = 4L * length(ETA_GRID) * length(SAMPLE_SIZES),
    design_sensitivity = 64L * length(SAMPLE_SIZES),
    visit_frequency = 24L * length(SAMPLE_SIZES),
    effect_tilt_sensitivity = 24L * length(SAMPLE_SIZES),
    effect_tilt_eta_calibration = 12L * length(ETA_GRID) * length(SAMPLE_SIZES),
    two_scheme_weight_battery = 120L * length(SAMPLE_SIZES)
  )
  observed_rows <- c(
    expanded_primary = nrow(primary),
    expanded_primary_eta_calibration = nrow(primary_cal),
    design_sensitivity = nrow(design),
    visit_frequency = nrow(visit),
    effect_tilt_sensitivity = nrow(effect),
    effect_tilt_eta_calibration = nrow(effect_cal),
    two_scheme_weight_battery = nrow(weights)
  )
  row_counts_match <- identical(as.integer(observed_rows),
                                as.integer(expected_rows))
  names(row_counts_match) <- NULL

  schemas <- list(
    expanded_primary = c(
      "study", "scenario", "scheme", "truth", "procedure",
      "rejection_probability", "mc_se", "median_stop_given_reject",
      "mean_stop_capped", "mean_terminal_loge", "mean_maximum_loge", "mean_total_events",
      "mean_updates", "right_censored_fraction", "horizon", "n_reps",
      "n", "allocation_treatment", "dropout_probability", "visit_width",
      "target_weight", "eta"
    ),
    design_sensitivity = c(
      "study", "scenario", "scheme", "truth", "procedure",
      "rejection_probability", "mc_se", "median_stop_given_reject",
      "mean_stop_capped", "mean_terminal_loge", "mean_maximum_loge", "mean_total_events",
      "mean_updates", "right_censored_fraction", "horizon", "n_reps",
      "n", "allocation_treatment", "dropout_probability", "visit_width",
      "target_weight", "eta"
    ),
    visit_frequency = c(
      "study", "scenario", "scheme", "truth", "procedure",
      "rejection_probability", "mc_se", "median_stop_given_reject",
      "mean_stop_capped", "mean_terminal_loge", "mean_maximum_loge", "mean_total_events",
      "mean_updates", "right_censored_fraction", "horizon", "n_reps",
      "n", "allocation_treatment", "dropout_probability", "visit_width",
      "n_assessments", "target_weight", "eta"
    ),
    effect_tilt_sensitivity = c(
      "n", "scenario_family", "effect", "scheme", "tilt_choice", "eta",
      "power", "mc_se", "median_stop_given_reject", "mean_stop_capped",
      "mean_terminal_loge", "mean_maximum_loge", "right_censored_fraction", "n_reps"
    ),
    two_scheme_weight_battery = c(
      "scenario", "scheme", "truth", "weight", "eta",
      "rejection_probability", "mc_se", "median_stop_given_reject",
      "mean_stop_capped", "mean_terminal_loge", "mean_maximum_loge", "n_reps", "n", "horizon"
    )
  )
  observed <- list(
    expanded_primary = primary,
    design_sensitivity = design,
    visit_frequency = visit,
    effect_tilt_sensitivity = effect,
    two_scheme_weight_battery = weights
  )
  schema_match <- vapply(names(schemas), function(name) {
    identical(names(observed[[name]]), schemas[[name]])
  }, logical(1L))
  probabilities_finite <- all(is.finite(c(
    primary$rejection_probability,
    design$rejection_probability,
    visit$rejection_probability,
    effect$power,
    weights$rejection_probability
  )))
  sample_sizes_match <- all(vapply(c(observed, list(primary_cal, effect_cal)),
    function(x) identical(sort(unique(as.integer(x$n))), sort(SAMPLE_SIZES)), logical(1L)))
  passed <- isTRUE(row_counts_match) && all(schema_match) &&
    probabilities_finite && sample_sizes_match
  if (!passed) {
    stop(
      paste0(
        "Reduced practical-battery smoke validation failed. Row counts: ",
        paste(names(observed_rows), observed_rows, sep = "=", collapse = ", "),
        "; schema matches: ",
        paste(names(schema_match), schema_match, sep = "=", collapse = ", "),
        "; finite probabilities: ", probabilities_finite
      ),
      call. = FALSE
    )
  }
  list(
    passed = TRUE,
    expected_rows = as.list(expected_rows),
    observed_rows = as.list(observed_rows),
    schema_match = as.list(schema_match),
    probabilities_finite = probabilities_finite,
    sample_sizes_match = sample_sizes_match
  )
}

package_versions <- function(packages) {
  out <- lapply(packages, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }
  })
  names(out) <- packages
  out
}

run_battery <- function(args) {
  start <- proc.time()[["elapsed"]]
  output_dir <- normalizePath(args$output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  primary_result <- run_primary(args, args$seed)
  primary <- primary_result$data
  primary_cal <- primary_result$calibration
  write_result_csv(primary, file.path(output_dir, "expanded_primary.csv"))
  write_result_csv(primary_cal,
                   file.path(output_dir, "expanded_primary_eta_calibration.csv"))

  base_eta <- list()
  baseline_calibration <- list()
  base_specs <- list(
    list(
      family = "ph",
      sc = ENGINE$new_mixed_scenario(
        "ph_base", "", "logrank", REFERENCE_N, 8L, 0.21, 0.15, 0.15,
        0.70, 0.14
      ),
      weight = "logrank"
    ),
    list(
      family = "cure",
      sc = ENGINE$new_mixed_scenario(
        "cure_base", "", "gray_tsiatis", REFERENCE_N, 10L, 0.25, 0.15,
        0.34, 1, 0.12
      ),
      weight = "gray_tsiatis"
    )
  )
  for (base_spec in base_specs) {
    for (scheme in c("interval", "exact")) {
      design <- calibrate_design(
        base_spec$sc, scheme, base_spec$weight, 0.5,
        args$seed + 10000000L + match(base_spec$family, c("ph", "cure")) * 1000L, args$pilot
      )
      base_eta[[paste(base_spec$family, scheme, sep = "_")]] <- design$eta
      curve <- design$curve
      curve$n <- REFERENCE_N
      curve$family <- base_spec$family
      curve$scheme <- scheme
      baseline_calibration[[length(baseline_calibration) + 1L]] <- curve
    }
  }

  write_result_csv(bind_frames(baseline_calibration),
                   file.path(output_dir, "baseline_eta_calibration.csv"))

  design <- run_design_grid(args, args$seed + 20000000L, base_eta)
  write_result_csv(design, file.path(output_dir, "design_sensitivity.csv"))

  visit <- run_visit_frequency(
    args, args$seed + 30000000L, base_eta[["ph_interval"]]
  )
  write_result_csv(visit, file.path(output_dir, "visit_frequency.csv"))

  effect_result <- run_effect_and_tilt(args, args$seed + 40000000L, base_eta)
  effect <- effect_result$data
  effect_cal <- effect_result$calibration
  write_result_csv(effect,
                   file.path(output_dir, "effect_tilt_sensitivity.csv"))
  write_result_csv(effect_cal,
                   file.path(output_dir, "effect_tilt_eta_calibration.csv"))

  weight_result <- run_two_scheme_weight_battery(
    args, args$seed + 50000000L
  )
  two_scheme_weights <- weight_result$data
  write_result_csv(weight_result$calibration,
                   file.path(output_dir, "two_scheme_weight_eta_calibration.csv"))
  write_result_csv(two_scheme_weights,
                   file.path(output_dir, "two_scheme_weight_battery.csv"))
  validation <- validate_outputs(primary, design, visit, effect)
  weight_null <- two_scheme_weights[
    two_scheme_weights$truth == "Null", , drop = FALSE
  ]
  weight_check <- null_check(weight_null)
  validation$two_scheme_weight_null_below_0_05 <-
    weight_check$all_estimates_at_or_below_alpha
  validation$two_scheme_weight_null_within_1_96_mc_se <-
    weight_check$all_estimates_within_mc_tolerance
  validation$two_scheme_weight_null_diagnostics <- weight_check
  validation$primary_exact_interval_common_latent <-
    primary_result$exact_interval_common_latent
  validation$intentional_design_correction <- paste(
    "Exact and interval analyses apply both observation maps to the same",
    "latent pilot and evaluation datasets within each N/scenario/truth cell.",
    "OBF information planning and monitoring align on scheduled calendar looks."
  )
  {
    validation$structure <- validate_smoke_outputs(
      primary, primary_cal, design, visit, effect, effect_cal,
      two_scheme_weights
    )
  }
  validation$r_version <- R.version.string
  validation$platform <- paste(R.version$platform, R.version$os, sep = "; ")
  validation$rng_kind <- as.list(RNGkind())
  names(validation$rng_kind) <- c("kind", "normal_kind", "sample_kind")
  validation$packages <- package_versions(
    c("Rcpp", "data.table", "jsonlite", "matrixStats")
  )
  validation$seed <- args$seed
  validation$settings <- args
  validation$baseline_eta <- base_eta
  validation$sample_sizes <- SAMPLE_SIZES
  validation$reference_n <- REFERENCE_N
  validation$weight_calibration_rows <- nrow(weight_result$calibration)
  stopifnot(nrow(weight_result$calibration) == 50L * length(ETA_GRID) * length(SAMPLE_SIZES),
            nrow(bind_frames(baseline_calibration)) == 4L * length(ETA_GRID))
  validation$paired_schemes_all_studies <- TRUE
  validation$obf_looks <- "Prespecified scheduled assessment calendar times for both observation schemes; cumulative scores carried forward through event-free visits."
  validation$pilot_design <- "Independent null pilot for weight scaling/information; independent alternative pilot for tilt; independent evaluation. Sensitivity baseline tilts fixed at reference N across N and effect cells."
  validation$elapsed_seconds <- proc.time()[["elapsed"]] - start
  validation$session_info <- capture.output(utils::sessionInfo())

  validation_path <- file.path(output_dir, "expanded_validation.json")
  jsonlite::write_json(
    validation, validation_path, pretty = TRUE, auto_unbox = TRUE,
    digits = NA, na = "null"
  )
  cat(paste(readLines(validation_path, warn = FALSE), collapse = "\n"), "\n")
  invisible(list(
    primary = primary,
    design = design,
    visit = visit,
    effect = effect,
    two_scheme_weights = two_scheme_weights,
    validation = validation
  ))
}

usage <- function() {
  cat(paste0(
    "Usage: Rscript code/run_practical_battery.R [options]\n\n",
    "Options:\n",
    "  --seed INTEGER\n",
    "  --pilot INTEGER\n",
    "  --pilot-small INTEGER\n",
    "  --eval-primary INTEGER\n",
    "  --eval-grid INTEGER\n",
    "  --eval-visit INTEGER\n",
    "  --eval-effect INTEGER\n",
    "  --eval-weight INTEGER\n",
    "  --output-dir PATH\n",
    "  --smoke-test       Cap all replication counts at small test values.\n",
    "  --help\n"
  ))
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    seed = 20260901L,
    pilot = PILOT_REPS,
    pilot_small = PILOT_REPS,
    eval_primary = EVALUATION_REPS,
    eval_grid = EVALUATION_REPS,
    eval_visit = EVALUATION_REPS,
    eval_effect = EVALUATION_REPS,
    eval_weight = EVALUATION_REPS,
    output_dir = file.path(ROOT, "results"),
    smoke_test = FALSE
  )
  option_names <- list(
    seed = "seed",
    pilot = "pilot",
    `pilot-small` = "pilot_small",
    `eval-primary` = "eval_primary",
    `eval-grid` = "eval_grid",
    `eval-visit` = "eval_visit",
    `eval-effect` = "eval_effect",
    `eval-weight` = "eval_weight",
    `output-dir` = "output_dir"
  )
  args <- defaults
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (token == "--help") {
      usage()
      quit(save = "no", status = 0)
    }
    if (token == "--smoke-test") {
      args$smoke_test <- TRUE
      i <- i + 1L
      next
    }
    if (!startsWith(token, "--")) {
      stop(sprintf("Unexpected positional argument: %s", token), call. = FALSE)
    }
    stripped <- substring(token, 3L)
    if (grepl("=", stripped, fixed = TRUE)) {
      pieces <- strsplit(stripped, "=", fixed = TRUE)[[1L]]
      key <- pieces[[1L]]
      value <- paste(pieces[-1L], collapse = "=")
    } else {
      key <- stripped
      i <- i + 1L
      if (i > length(argv)) {
        stop(sprintf("Option --%s requires a value", key), call. = FALSE)
      }
      value <- argv[[i]]
    }
    target <- option_names[[key]]
    if (is.null(target)) {
      stop(sprintf("Unknown option: --%s", key), call. = FALSE)
    }
    if (is.integer(defaults[[target]])) {
      parsed <- suppressWarnings(as.integer(value))
      if (is.na(parsed) || parsed <= 0L) {
        stop(sprintf("--%s must be a positive integer", key), call. = FALSE)
      }
      args[[target]] <- parsed
    } else {
      args[[target]] <- value
    }
    i <- i + 1L
  }

  if (args$smoke_test) {
    args$pilot <- min(args$pilot, 16L)
    args$pilot_small <- min(args$pilot_small, 12L)
    args$eval_primary <- min(args$eval_primary, 24L)
    args$eval_grid <- min(args$eval_grid, 16L)
    args$eval_visit <- min(args$eval_visit, 20L)
    args$eval_effect <- min(args$eval_effect, 20L)
    args$eval_weight <- min(args$eval_weight, 20L)
  }
  args
}

main <- function() {
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  args <- parse_args()
  run_battery(args)
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  main()
}
