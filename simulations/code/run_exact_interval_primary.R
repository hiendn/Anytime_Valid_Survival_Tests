#!/usr/bin/env Rscript

# Native-R mixed-resolution simulation study for weighted survival e-processes.
#
# The entry point is source-safe: functions are defined when sourced, and the
# production run starts only when this file is invoked directly with Rscript.

.exact_interval_script <- local({
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
  matching <- candidates[basename(candidates) == "run_exact_interval_primary.R"]
  if (length(matching)) {
    normalizePath(matching[[length(matching)]], winslash = "/", mustWork = FALSE)
  } else {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg)) {
      normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE), winslash = "/", mustWork = FALSE)
    } else {
      normalizePath(file.path(getwd(), "code", "run_exact_interval_primary.R"), winslash = "/", mustWork = FALSE)
    }
  }
})

EXACT_INTERVAL_ROOT <- dirname(dirname(.exact_interval_script))
source(file.path(dirname(.exact_interval_script), "simulation_config.R"), local = environment())
if (!exists("new_mixed_scenario", mode = "function", inherits = TRUE)) {
  source(
    file.path(dirname(.exact_interval_script), "simulation_core.R"),
    local = environment()
  )
}


exact_interval_usage <- function() {
  paste(
    "Usage: Rscript code/run_exact_interval_primary.R [options]",
    "",
    "Core options:",
    "  --seed INTEGER            Principal seed (default: 20260830)",
    "  --pilot-reps INTEGER      Alternative pilot replications (default: 1000)",
    "  --eval-reps INTEGER       Evaluation replications per truth (default: 5000)",
    "  --sample-sizes CSV        Subset of the common N grid (default: 100,200,500,1000,2000)",
    "  --eval-chunk-size INTEGER Evaluation batch size (default: 250)",
    "  --fractions CSV           Exact fractions (default: 0,0.25,0.5,0.75,1)",
    "  --eta-min NUMBER          Lower eta-grid endpoint (default: -1.25)",
    "  --eta-max NUMBER          Upper eta-grid endpoint (default: -0.015)",
    "  --eta-grid INTEGER        Number of eta-grid points (default: 50)",
    "  --scenarios CSV           Scenario names (default: latency_ph,pure_cure)",
    "",
    "Reproducible presets and validation:",
    "  --production              Manuscript settings: pilot 1000, eval 5000, grid 50",
    "  --smoke                   Small full-pipeline settings for development",
    "  --self-test               Run focused core tests and exit",
    "  --identity-reps INTEGER   Replications for endpoint identity checks (default: 200)",
    "  --validation-iterations N Algebra-validation iterations (default: 100)",
    "",
    "Output options:",
    "  --output-dir PATH         CSV/JSON directory (default: results/)",
    "  --figures-dir PATH        PDF directory (default: manuscript/figures/)",
    "  --skip-figures            Do not generate the five PDF figures",
    "  --help                    Show this help",
    sep = "\n"
  )
}


.parse_scalar <- function(arguments, index, option) {
  if (index >= length(arguments)) stop(option, " requires a value", call. = FALSE)
  arguments[[index + 1L]]
}


parse_exact_interval_args <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  production <- "--production" %in% arguments
  smoke <- "--smoke" %in% arguments
  if (production && smoke) stop("--production and --smoke are mutually exclusive", call. = FALSE)

  settings <- list(
    seed = 20260830L,
    pilot_reps = if (smoke) 12L else PILOT_REPS,
    eval_reps = if (smoke) 24L else EVALUATION_REPS,
    sample_sizes = paste(SAMPLE_SIZES, collapse = ","),
    eval_chunk_size = 250L,
    fractions = "0,0.25,0.5,0.75,1",
    eta_min = -1.25,
    eta_max = -0.015,
    eta_grid = if (smoke) 5L else ETA_GRID_SIZE,
    scenarios = "latency_ph,pure_cure",
    identity_reps = if (smoke) 12L else 200L,
    validation_iterations = if (smoke) 20L else 100L,
    output_dir = file.path(EXACT_INTERVAL_ROOT, "results"),
    figures_dir = file.path(EXACT_INTERVAL_ROOT, "manuscript", "figures"),
    skip_figures = FALSE,
    production = production,
    smoke = smoke,
    self_test = "--self-test" %in% arguments,
    help = "--help" %in% arguments
  )

  value_options <- c(
    "--seed", "--pilot-reps", "--eval-reps", "--fractions", "--sample-sizes", "--eval-chunk-size",
    "--eta-min", "--eta-max", "--eta-grid", "--scenarios",
    "--identity-reps", "--validation-iterations", "--output-dir", "--figures-dir"
  )
  flag_options <- c("--production", "--smoke", "--self-test", "--skip-figures", "--help")
  index <- 1L
  while (index <= length(arguments)) {
    option <- arguments[[index]]
    if (option %in% flag_options) {
      if (option == "--skip-figures") settings$skip_figures <- TRUE
      index <- index + 1L
      next
    }
    if (!option %in% value_options) stop("Unknown option: ", option, call. = FALSE)
    value <- .parse_scalar(arguments, index, option)
    key <- switch(
      option,
      "--seed" = "seed",
      "--sample-sizes" = "sample_sizes",
      "--eval-chunk-size" = "eval_chunk_size",
      "--pilot-reps" = "pilot_reps",
      "--eval-reps" = "eval_reps",
      "--fractions" = "fractions",
      "--eta-min" = "eta_min",
      "--eta-max" = "eta_max",
      "--eta-grid" = "eta_grid",
      "--scenarios" = "scenarios",
      "--identity-reps" = "identity_reps",
      "--validation-iterations" = "validation_iterations",
      "--output-dir" = "output_dir",
      "--figures-dir" = "figures_dir"
    )
    if (key %in% c("seed", "pilot_reps", "eval_reps", "eta_grid", "identity_reps", "validation_iterations", "eval_chunk_size")) {
      converted <- suppressWarnings(as.integer(value))
      if (is.na(converted)) stop(option, " requires an integer", call. = FALSE)
      settings[[key]] <- converted
    } else if (key %in% c("eta_min", "eta_max")) {
      converted <- suppressWarnings(as.numeric(value))
      if (!is.finite(converted)) stop(option, " requires a finite number", call. = FALSE)
      settings[[key]] <- converted
    } else {
      settings[[key]] <- value
    }
    index <- index + 2L
  }

  positive <- c("pilot_reps", "eval_reps", "eta_grid", "identity_reps", "validation_iterations", "eval_chunk_size")
  if (any(vapply(settings[positive], function(value) value <= 0L, logical(1)))) {
    stop("replication counts, eta-grid size, and validation counts must be positive", call. = FALSE)
  }
  ns <- suppressWarnings(as.numeric(strsplit(settings$sample_sizes, ",", fixed = TRUE)[[1L]]))
  if (!length(ns) || any(!is.finite(ns)) || any(!ns %in% SAMPLE_SIZES) || anyDuplicated(ns)) {
    stop("--sample-sizes must be a distinct subset of the common sample-size grid: ", paste(SAMPLE_SIZES, collapse = ","), call. = FALSE)
  }
  settings$sample_size_values <- sort(as.integer(ns))
  fractions <- suppressWarnings(as.numeric(strsplit(settings$fractions, ",", fixed = TRUE)[[1L]]))
  if (!length(fractions) || any(!is.finite(fractions)) || any(fractions < 0 | fractions > 1)) {
    stop("--fractions must be a comma-separated list in [0, 1]", call. = FALSE)
  }
  if (anyDuplicated(fractions)) stop("--fractions must not contain duplicates", call. = FALSE)
  if (!settings$skip_figures && !any(abs(fractions - 0.5) < 1e-12)) {
    stop("figure generation requires an exact fraction of 0.5; add it or use --skip-figures", call. = FALSE)
  }
  settings$fraction_values <- fractions
  scenario_names <- trimws(strsplit(settings$scenarios, ",", fixed = TRUE)[[1L]])
  if (!length(scenario_names) || any(!scenario_names %in% names(SCENARIOS))) {
    stop("unknown scenario in --scenarios; available: ", paste(names(SCENARIOS), collapse = ", "), call. = FALSE)
  }
  if (!settings$skip_figures && !all(names(SCENARIOS) %in% scenario_names)) {
    stop("all five figures require both standard scenarios; use both or --skip-figures", call. = FALSE)
  }
  settings$scenario_names <- scenario_names
  settings$output_dir <- normalizePath(settings$output_dir, winslash = "/", mustWork = FALSE)
  settings$figures_dir <- normalizePath(settings$figures_dir, winslash = "/", mustWork = FALSE)
  settings
}


.chosen_key <- function(rho, method) paste0(sprintf("%.17g", rho), "|", method)


.append_columns <- function(data, values) {
  for (name in names(values)) data[[name]] <- values[[name]]
  data
}


.seed_schedule <- function(seed, scenario_names, sample_sizes, n_chunks) {
  schedule <- list(algebra_validation = 9817L, scenario = list())
  for (name in scenario_names) {
    scenario_id <- match(name, c("latency_ph", "pure_cure")) - 1L
    schedule$scenario[[name]] <- list()
    for (n in sample_sizes) {
      base <- as.double(seed) + scenario_id * 10000000 + n * 1000
      seeds <- list(
        pilot_alternative = as.integer(base + 1),
        evaluation_null = as.integer(base + 1000 + seq_len(n_chunks)),
        evaluation_alternative = as.integer(base + 2000 + seq_len(n_chunks)),
        endpoint_identity = as.integer(base + 9001)
      )
      if (anyNA(unlist(seeds))) stop("seed exceeds the R integer range", call. = FALSE)
      schedule$scenario[[name]][[as.character(n)]] <- seeds
    }
  }
  all_seeds <- unlist(schedule$scenario, use.names = FALSE)
  if (anyDuplicated(all_seeds)) stop("seed schedule contains collisions", call. = FALSE)
  schedule
}


# Only small vectors are retained across batches; every median is computed from
# all evaluation replications, rather than by averaging batch medians.
.mixed_evaluation <- function(updates, eta, scenario, keep_path = FALSE) {
  evaluation <- evaluate(updates, eta, scenario$target_weight, scenario$horizon)
  path <- NULL
  if (keep_path) {
    columns <- seq_len(updates$count[[1L]])
    path <- list(time = as.numeric(updates$time[1L, columns]),
                 loge = as.numeric(evaluation$loge_path[1L, columns]))
  }
  evaluation$loge_path <- NULL
  list(evaluation = evaluation, count = updates$count, path = path)
}


.join_mixed_batches <- function(batches) {
  fields <- names(batches[[1L]]$evaluation)
  list(
    evaluation = setNames(lapply(fields, function(field) {
      unlist(lapply(batches, function(batch) batch$evaluation[[field]]), use.names = FALSE)
    }), fields),
    count = unlist(lapply(batches, function(batch) batch$count), use.names = FALSE)
  )
}


run_exact_interval_primary <- function(settings) {
  start <- proc.time()[["elapsed"]]
  configure_simulation_rng()
  ensure_simulation_kernels()
  fractions <- settings$fraction_values
  sample_sizes <- settings$sample_size_values
  eta_grid <- seq(settings$eta_min, settings$eta_max, length.out = settings$eta_grid)
  scenarios <- SCENARIOS[settings$scenario_names]
  n_chunks <- ceiling(settings$eval_reps / settings$eval_chunk_size)
  chunk_sizes <- pmin(settings$eval_chunk_size,
                      settings$eval_reps - (seq_len(n_chunks) - 1L) * settings$eval_chunk_size)
  seed_schedule <- .seed_schedule(settings$seed, settings$scenario_names, sample_sizes, n_chunks)
  summary_rows <- calibration_rows <- paired_rows <- observation_rows <- list()
  examples <- list()
  validation <- validate_algebra(seed = seed_schedule$algebra_validation,
                                 iterations = settings$validation_iterations)

  for (scenario_name in names(scenarios)) for (n in sample_sizes) {
    scenario <- scenarios[[scenario_name]]
    scenario$n <- n
    message("Scenario ", scenario_name, ", total N = ", n,
            ": calibration (", settings$pilot_reps, " replications)")
    scenario_seeds <- seed_schedule$scenario[[scenario_name]][[as.character(n)]]
    pilot <- simulate_latent(scenario, settings$pilot_reps, scenario_seeds$pilot_alternative)
    chosen <- list()
    calibrate <- function(rho, method, reported_rho = rho) {
      updates <- build_updates(pilot, scenario, rho, method)
      fit <- calibrate_eta(updates, scenario$target_weight, eta_grid)
      calibration_rows[[length(calibration_rows) + 1L]] <<- .append_columns(
        fit$curve, list(scenario = scenario_name, n = n, exact_fraction = reported_rho,
                        method = method, weight = scenario$target_weight,
                        n_reps = settings$pilot_reps,
                        seed = scenario_seeds$pilot_alternative))
      fit$eta
    }
    eta_panel <- calibrate(0, "all_panel", -1)
    eta_oracle <- calibrate(1, "all_exact_oracle", -1)
    for (rho in fractions) {
      chosen[[.chosen_key(rho, "all_panel")]] <- eta_panel
      chosen[[.chosen_key(rho, "all_exact_oracle")]] <- eta_oracle
      if (rho <= 0) {
        chosen[[.chosen_key(rho, "mixed")]] <- eta_panel
        chosen[[.chosen_key(rho, "exact_only")]] <- eta_grid[[floor(length(eta_grid) / 2) + 1L]]
      } else if (rho >= 1) {
        chosen[[.chosen_key(rho, "mixed")]] <- eta_oracle
        chosen[[.chosen_key(rho, "exact_only")]] <- eta_oracle
      } else for (method in c("mixed", "exact_only")) {
        chosen[[.chosen_key(rho, method)]] <- calibrate(rho, method)
      }
    }
    rm(pilot)
    invisible(gc(FALSE))
    for (truth in c("H0", "H1")) {
      data_scenario <- if (truth == "H0") null_version(scenario) else scenario
      evaluation_seeds <- if (truth == "H0") scenario_seeds$evaluation_null else scenario_seeds$evaluation_alternative
      batches <- list()
      observation_batches <- setNames(lapply(fractions, function(rho) list()), as.character(fractions))
      for (chunk in seq_len(n_chunks)) {
        message("  ", truth, ", batch ", chunk, "/", n_chunks,
                " (", chunk_sizes[[chunk]], " replications)")
        latent <- simulate_latent(data_scenario, chunk_sizes[[chunk]], evaluation_seeds[[chunk]])
        keep_path <- truth == "H1" && chunk == 1L
        panel_updates <- build_updates(latent, data_scenario, 0, "all_panel")
        panel <- .mixed_evaluation(panel_updates, eta_panel, scenario, keep_path)
        rm(panel_updates)
        oracle_updates <- build_updates(latent, data_scenario, 1, "all_exact_oracle")
        oracle <- .mixed_evaluation(oracle_updates, eta_oracle, scenario, keep_path)
        rm(oracle_updates)
        for (rho in fractions) {
          mix <- observation_mix(latent, rho, scenario$horizon)
          observation_batches[[as.character(rho)]][[chunk]] <- unlist(mix)
          if (rho <= 0) {
            mixed <- panel
            updates <- build_updates(latent, data_scenario, 0, "exact_only")
            exact <- .mixed_evaluation(updates, chosen[[.chosen_key(rho, "exact_only")]], scenario, keep_path)
            rm(updates)
          } else if (rho >= 1) {
            mixed <- oracle
            exact <- oracle
          } else {
            updates <- build_updates(latent, data_scenario, rho, "mixed")
            mixed <- .mixed_evaluation(updates, chosen[[.chosen_key(rho, "mixed")]], scenario, keep_path)
            rm(updates)
            updates <- build_updates(latent, data_scenario, rho, "exact_only")
            exact <- .mixed_evaluation(updates, chosen[[.chosen_key(rho, "exact_only")]], scenario, keep_path)
            rm(updates)
          }
          evaluations <- list(mixed = mixed, all_panel = panel, exact_only = exact, all_exact_oracle = oracle)
          if (keep_path && abs(rho - 0.5) < 1e-8) {
            examples[[paste(scenario_name, n, sep = "|")]] <- lapply(evaluations, function(result) result$path)
          }
          for (method in METHODS) {
            key <- .chosen_key(rho, method)
            result <- evaluations[[method]]
            result$path <- NULL
            batches[[key]][[chunk]] <- result
          }
          rm(mixed, exact, evaluations)
        }
        rm(latent, panel, oracle)
        invisible(gc(FALSE))
      }
      for (rho in fractions) {
        obs_matrix <- do.call(rbind, observation_batches[[as.character(rho)]])
        mix <- colSums(obs_matrix * as.numeric(chunk_sizes)) / settings$eval_reps
        observation_rows[[length(observation_rows) + 1L]] <- data.frame(
          scenario = scenario_name, n = n, truth = truth, exact_fraction = rho,
          as.list(mix), n_reps = settings$eval_reps, check.names = FALSE)
        evaluations <- list()
        for (method in METHODS) {
          joined <- .join_mixed_batches(batches[[.chosen_key(rho, method)]])
          evaluations[[method]] <- joined$evaluation
          summary_rows[[length(summary_rows) + 1L]] <- .append_columns(
            summarize_mixed_resolution(scenario_name, truth, rho, method, scenario$target_weight,
              chosen[[.chosen_key(rho, method)]], joined$evaluation, list(count = joined$count)),
            list(n = n))
        }
        if (truth == "H1") for (comparison in list(c("mixed", "all_panel"), c("mixed", "exact_only"), c("all_exact_oracle", "mixed"))) {
          a <- comparison[[1L]]
          b <- comparison[[2L]]
          paired_rows[[length(paired_rows) + 1L]] <- .append_columns(
            paired_difference(scenario_name, rho, a, evaluations[[a]]$reject, b, evaluations[[b]]$reject), list(n = n))
        }
      }
      rm(batches, observation_batches)
      invisible(gc(FALSE))
    }
    identity_data <- simulate_latent(scenario, settings$identity_reps, scenario_seeds$endpoint_identity)
    for (check in list(list(0, "mixed", "all_panel", "rho0_mixed_equals_panel"),
                       list(1, "mixed", "all_exact_oracle", "rho1_mixed_equals_oracle"))) {
      first <- build_updates(identity_data, scenario, check[[1L]], check[[2L]])
      second <- build_updates(identity_data, scenario, check[[1L]], check[[3L]])
      equal <- all(vapply(c("n1", "n0", "d", "x", "time", "count"),
                           function(field) identical(first[[field]], second[[field]]), logical(1)))
      validation[[paste0(scenario_name, "_n", n, "_", check[[4L]])]] <- equal
      rm(first, second)
    }
    rm(identity_data)
    invisible(gc(FALSE))
  }
  summary <- do.call(rbind, summary_rows)
  calibration <- do.call(rbind, calibration_rows)
  paired <- do.call(rbind, paired_rows)
  observation <- do.call(rbind, observation_rows)
  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  write_result_csv(summary, file.path(settings$output_dir, "mixed_resolution_summary.csv"))
  write_result_csv(calibration, file.path(settings$output_dir, "mixed_resolution_eta_calibration.csv"))
  write_result_csv(paired, file.path(settings$output_dir, "mixed_resolution_paired_differences.csv"))
  write_result_csv(observation, file.path(settings$output_dir, "mixed_resolution_observation_mix.csv"))
  saveRDS(examples, file.path(settings$output_dir, "mixed_resolution_example_paths.rds"))
  figure_paths <- character()
  if (!settings$skip_figures) {
    source(file.path(dirname(.exact_interval_script), "mixed_resolution_figures.R"), local = environment())
    figure_paths <- make_mixed_resolution_grid_figures(summary, observation, examples, settings$figures_dir)
  }
  equality_names <- grep("equals", names(validation), value = TRUE)
  validation$all_endpoint_identities_hold <- all(unlist(validation[equality_names], use.names = FALSE))
  if (!isTRUE(validation$algebra_within_tolerance) || !isTRUE(validation$all_endpoint_identities_hold)) {
    stop("mixed-resolution algebra or endpoint validation failed", call. = FALSE)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  validation <- c(validation, collect_r_metadata(c("Rcpp", "data.table", "jsonlite")), list(
    seed = settings$seed, seed_schedule = seed_schedule,
    settings = list(seed = settings$seed, sample_sizes = sample_sizes,
      sample_size_definition = "total randomized subjects; independent Bernoulli(0.5) allocation",
      pilot_reps = settings$pilot_reps, eval_reps = settings$eval_reps,
      eval_chunk_size = settings$eval_chunk_size, evaluation_chunk_sizes = chunk_sizes,
      fractions = settings$fractions, eta_min = settings$eta_min, eta_max = settings$eta_max,
      eta_grid = settings$eta_grid, scenarios = settings$scenarios,
      identity_reps = settings$identity_reps, validation_iterations = settings$validation_iterations,
      production = settings$production, smoke = settings$smoke),
    elapsed_seconds = elapsed, generated_figures = unname(basename(figure_paths))))
  write_result_json(validation, file.path(settings$output_dir, "mixed_resolution_validation.json"))
  message(sprintf("Completed in %.1f seconds", elapsed))
  invisible(list(summary = summary, calibration = calibration, paired = paired,
                 observation = observation, validation = validation, figure_paths = figure_paths))
}


.run_exact_interval_as_main <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) return(FALSE)
  invoked <- normalizePath(
    gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE),
    winslash = "/",
    mustWork = FALSE
  )
  if (!identical(basename(invoked), "run_exact_interval_primary.R")) return(FALSE)
  identical(invoked, normalizePath(.exact_interval_script, winslash = "/", mustWork = FALSE))
}


exact_interval_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  settings <- parse_exact_interval_args(arguments)
  if (settings$help) {
    cat(exact_interval_usage(), "\n")
    return(invisible(NULL))
  }
  if (settings$self_test) {
    result <- simulation_core_self_test(stop_on_failure = TRUE)
    print(result)
    return(invisible(result))
  }
  run_exact_interval_primary(settings)
}


if (.run_exact_interval_as_main()) {
  tryCatch(
    exact_interval_main(),
    error = function(error) {
      message("ERROR: ", conditionMessage(error))
      quit(save = "no", status = 1L, runLast = FALSE)
    }
  )
}
