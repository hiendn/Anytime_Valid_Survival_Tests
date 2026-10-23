#!/usr/bin/env Rscript
# Isolated SIM-1/SIM-2 experiment. Does not edit manuscripts or old results.

session_root <- function() {
  ff <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(ff)) stop("Run this script with Rscript")
  dirname(dirname(normalizePath(sub("^--file=", "", ff[[1L]]))))
}

scenario_specs <- function() list(
  S1 = list(label = "Sustained susceptible-hazard effect", K = 8L, delta = 1,
            rate = .21, c0 = .15, c1 = .15, early = .65, late = .65,
            change = Inf, dropout = .14),
  S2 = list(label = "Pure cure-probability effect", K = 10L, delta = 1,
            rate = .25, c0 = .15, c1 = .34, early = 1, late = 1,
            change = Inf, dropout = .12),
  S3 = list(label = "Early effect", K = 12L, delta = .5,
            rate = .20 * .5, c0 = 0, c1 = 0, early = .65, late = 1,
            change = 2.5 / .5, dropout = 1 - exp(-.07 * .5)),
  S4 = list(label = "Delayed effect", K = 16L, delta = .5,
            rate = .18 * .5, c0 = 0, c1 = 0, early = 1, late = .55,
            change = 3 / .5, dropout = 1 - exp(-.06 * .5)),
  S5 = list(label = "Joint latency and cure effect", K = 10L, delta = 1,
            rate = .23, c0 = .20, c1 = .30, early = .82, late = .82,
            change = Inf, dropout = .035)
)

parse_session_args <- function() {
  argv <- commandArgs(TRUE)
  out <- list(mode = "smoke", output = NULL, reps = NULL, pilot = NULL,
              batch = 250L, n = NULL, scenarios = paste0("S", 1:5), resume = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[[i]])
    if (key == "resume") { out$resume <- TRUE; i <- i + 1L; next }
    if (!key %in% c("mode", "output", "reps", "pilot", "batch", "n", "scenarios")) stop("Unknown argument: ", argv[[i]])
    i <- i + 1L
    if (i > length(argv)) stop("Missing value for ", key)
    out[[key]] <- argv[[i]]
    i <- i + 1L
  }
  stopifnot(out$mode %in% c("smoke", "production"))
  if (is.null(out$reps)) out$reps <- if (out$mode == "production") 5000L else 32L
  if (is.null(out$pilot)) out$pilot <- if (out$mode == "production") 1000L else 80L
  if (is.null(out$n)) out$n <- if (out$mode == "production") "100,200,500,1000,2000" else "100,500"
  out$n <- as.integer(strsplit(out$n, ",", fixed = TRUE)[[1L]])
  out$scenarios <- strsplit(paste(out$scenarios, collapse = ","), ",", fixed = TRUE)[[1L]]
  for (k in c("reps", "pilot", "batch")) out[[k]] <- as.integer(out[[k]])
  stopifnot(all(out$n %in% c(100L, 200L, 500L, 1000L, 2000L)),
            all(out$scenarios %in% paste0("S", 1:5)), out$reps > 0L,
            out$pilot > 0L, out$batch > 0L, !is.null(out$output))
  out
}

trim_updates <- function(up) {
  keep <- seq_len(max(1L, up$count))
  for (nm in names(up)) if (is.matrix(up[[nm]])) up[[nm]] <- up[[nm]][, keep, drop = FALSE]
  up
}

wilson <- function(k, n) {
  z <- qnorm(.975); p <- k / n; den <- 1 + z * z / n
  center <- (p + z * z / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
  c(lower = max(0, center - half), upper = min(1, center + half))
}

main <- function() {
  root <- session_root()
  args <- parse_session_args()
  if (!startsWith(args$output, "/")) args$output <- file.path(root, args$output)
  out <- normalizePath(args$output, mustWork = FALSE)
  allowed <- normalizePath(file.path(root, "results", "session"), mustWork = TRUE)
  if (!startsWith(paste0(out, "/"), paste0(allowed, "/"))) stop("Outputs must be inside results/session/")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "batches"), showWarnings = FALSE)
  if (file.exists(file.path(out, "design.json")) && !args$resume) stop("Output already exists; use a fresh directory or --resume")
  lib <- file.path(root, ".local", "R-library")
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(matrixStats); library(digest) })
  Sys.setenv(SIMULATION_RCPP_CACHE_DIR = file.path(root, "build", "rcpp_cache"))
  B <- new.env(parent = globalenv())
  source(file.path(root, "simulations", "code", "run_practical_battery.R"), local = B)
  source(file.path(root, "scripts", "sim1_benchmark.R"), local = globalenv())
  source(file.path(root, "scripts", "sim2_comparators.R"), local = globalenv())
  B$CORE$ensure_simulation_kernels()
  specs <- scenario_specs()
  weight_names <- c(LR = "Logrank", GB = "Gehan--Breslow", PP = "Peto--Peto",
                    FH = "Fleming--Harrington (0,1)", GT = "Gray--Tsiatis")
  eta_fixed <- log(.7)
  eta_mix <- log(c(.5, .7, .9))
  make_sc <- function(sp, n, truth) B$ENGINE$new_mixed_scenario(
    name = "session", description = sp$label, target_weight = "logrank",
    n = n, horizon = sp$K, base_rate = sp$rate, cure_control = sp$c0,
    cure_treatment = if (truth == "null") sp$c0 else sp$c1,
    hazard_ratio = 1, dropout_probability = sp$dropout)
  generate <- function(sp, n, truth, reps, seed) {
    sc <- make_sc(sp, n, truth)
    data <- B$custom_latent_piecewise(sc, reps, seed, .5,
              if (truth == "null") 1 else sp$early,
              if (truth == "null") 1 else sp$late, sp$change)
    list(sc = sc, data = data)
  }
  seed_base <- 71000000L
  seed_for <- function(n, scenario, truth, batch) as.integer(seed_base +
    match(n, c(100, 200, 500, 1000, 2000)) * 10000000L +
    match(scenario, paste0("S", 1:5)) * 1000000L +
    match(truth, c("null", "alternative")) * 100000L + batch)
  source_files <- c(file.path(root, "scripts", c("run_remaining_simulations.R", "sim1_benchmark.R", "sim2_comparators.R")),
                    list.files(file.path(root, "scripts"), "\\.cpp$", full.names = TRUE),
                    file.path(root, "simulations", "code", c("simulation_core.R", "simulation_kernels.cpp", "run_practical_battery.R", "run_exact_interval_primary.R", "simulation_config.R")))
  hashes <- setNames(lapply(source_files, function(p) digest(p, algo = "sha256", file = TRUE)), substring(source_files, nchar(root) + 2L))
  design <- list(mode = args$mode, n = args$n, scenarios = args$scenarios,
    reps_per_truth_cell = args$reps, batch_size = args$batch,
    pilot_reps = args$pilot, pilot_n = 500L, pilot_scenario = "S1 null",
    pilot_seed_base = 900000L, evaluation_seed_base = seed_base,
    alpha = .05, fixed_eta = eta_fixed, mixture_eta = eta_mix,
    mixture_weights = c("LR", "FH", "GT"), mixture_masses = "1/9 for each weight/coefficient pair",
    gray_survival_floor = .12, scenarios_parameters = specs, source_sha256 = hashes,
    observation_pairing = "Same latent draw and dropout records for both schemes and all procedures",
    prequential = "Unrestricted original smoothed MLE plus separately labeled protective restriction; grouped FNCH extension",
    manuscripts_modified = FALSE)
  if (args$resume) {
    saved <- fromJSON(file.path(out, "design.json"), simplifyVector = FALSE)
    if (!identical(toJSON(saved, auto_unbox = TRUE, null = "null"),
                   toJSON(design, auto_unbox = TRUE, null = "null"))) stop("Resume design/source hashes differ")
  } else write_json(design, file.path(out, "design.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))

  # Null-only design stage, evaluated before and independently of every trial.
  scales_file <- file.path(out, "frozen_scales.rds")
  if (args$resume && file.exists(scales_file)) scales <- readRDS(scales_file) else {
    information <- lapply(c("interval", "exact"), function(x) setNames(rep(0, 5), names(weight_names)))
    names(information) <- c("interval", "exact")
    for (bb in seq_len(ceiling(args$pilot / args$batch))) {
      nr <- min(args$batch, args$pilot - (bb - 1L) * args$batch)
      gen <- generate(specs$S1, 500L, "null", nr, 900000L + bb)
      for (scheme in names(information)) {
        up <- trim_updates(B$build(gen$data, gen$sc, scheme))
        mats <- B$all_weight_matrices(up)
        vv <- B$ENGINE$hypergeom_mean_var(up$n1, up$n0, up$d)$var
        for (w in names(weight_names)) information[[scheme]][[w]] <- information[[scheme]][[w]] + sum(mats[[weight_names[[w]]]]^2 * vv)
      }
      rm(gen, up, mats, vv); gc(FALSE)
    }
    scales <- lapply(information, function(info) sqrt(info / info[["LR"]]))
    saveRDS(scales, scales_file)
    fwrite(rbindlist(lapply(names(scales), function(s) data.frame(scheme = s, weight = names(scales[[s]]), scale = scales[[s]], normalization_multiplier = 1 / scales[[s]], pilot_n = 500, pilot_reps = args$pilot))), file.path(out, "frozen_scales.csv"))
  }
  print(scales)

  record <- function(method, reject, stop, terminal = NA_real_, maximum = NA_real_, pvalue = NA_real_, score = NA_real_, variance = NA_real_) {
    data.table(method = method, replicate = seq_along(reject), reject = as.logical(reject),
      stop_time = stop, terminal_loge = terminal, maximum_loge = maximum,
      pvalue = pvalue, score = score, variance = variance)
  }
  # Mixtures average complete e-processes, including the initial value one.
  from_path <- function(path, up, horizon, delta) {
    cross <- path >= log(20)
    rejected <- rowSums(cross) > 0
    first <- max.col(cross * 1L, ties.method = "first")
    stop <- rep(horizon, length(rejected))
    rr <- which(rejected)
    if (length(rr)) stop[rr] <- up$time[cbind(rr, first[rr])] * delta
    list(reject = rejected, stop = stop, terminal = path[, ncol(path)], maximum = pmax(0, rowMaxs(path)))
  }
  evaluate_scheme <- function(gen, sp, scheme) {
    up <- trim_updates(B$build(gen$data, gen$sc, scheme))
    init1 <- rowSums(gen$data$group); init0 <- gen$sc$n - init1
    mats <- B$all_weight_matrices(up)
    horizon <- sp$K * sp$delta
    rows <- list(); component_paths <- list()
    for (w in names(weight_names)) {
      ev <- B$evaluate_weight_e(up, weight_names[[w]], scales[[scheme]][[w]], eta_fixed,
                              horizon, sp$delta, mats[[weight_names[[w]]]])
      method <- if (w == "LR") "LR_fixed_original" else paste0(w, "_fixed")
      rows[[method]] <- record(method, ev$reject, ev$stop, ev$terminal,
                              pmax(0, rowMaxs(ev$path)))
      if (w %in% c("LR", "FH", "GT")) {
        for (ei in seq_along(eta_mix)) {
          path <- if (abs(eta_mix[[ei]] - eta_fixed) < 1e-14) ev$path else B$CORE$row_cumsums(
            B$ENGINE$table_loge(up$n1, up$n0, up$d, up$x,
              eta_mix[[ei]] * mats[[weight_names[[w]]]] / scales[[scheme]][[w]]))
          component_paths[[paste(w, ei, sep = "_")]] <- path
        }
      }
    }
    for (mixture in c("LR_coefficient_mixture", "LR_FH_GT_mixture")) {
      selected <- if (mixture == "LR_coefficient_mixture") component_paths[startsWith(names(component_paths), "LR_")] else component_paths
      logsum <- selected[[1L]]
      for (j in seq.int(2L, length(selected))) logsum <- B$logaddexp_pair(logsum, selected[[j]])
      ev <- from_path(logsum - log(length(selected)), up, horizon, sp$delta)
      rows[[mixture]] <- record(mixture, ev$reject, ev$stop, ev$terminal, ev$maximum)
    }
    comparators <- sim2_comparators(up, initial_n1 = init1, initial_n0 = init0,
                                   horizon = sp$K, eta = eta_fixed)
    if (max(abs(comparators$fixed$terminal_loge - rows$LR_fixed_original$terminal_loge)) > 1e-7 ||
        any(comparators$fixed$reject != rows$LR_fixed_original$reject) ||
        max(abs(comparators$fixed$stop_time * sp$delta - rows$LR_fixed_original$stop_time)) > 1e-8)
      stop("Original fixed-alternative/LR identity failed")
    for (direction in c("unrestricted", "protective")) {
      ev <- comparators[[if (direction == "protective") "prequential_protective" else "prequential"]]
      method <- paste0("prequential_", direction)
      rows[[method]] <- record(method, ev$reject, ev$stop_time * sp$delta,
                              ev$terminal_loge, ev$maximum_loge)
    }
    score <- B$score_arrays(up, "logrank")
    ix <- cbind(seq_along(init1), pmax(up$count, 1L))
    U <- score$U[ix]; V <- score$V[ix]
    p <- ifelse(V > 0, pnorm(U / sqrt(pmax(V, 1e-300))), 1)
    rows$fixed_score_LR <- record("fixed_score_LR", p <= .05 & V > 0,
                                rep(horizon, length(p)), pvalue = p, score = U, variance = V)
    if (scheme == "interval") {
      ev <- sim1_fast_updates(up, n1 = init1, n0 = init0, alternative = "less")
      rows$Fay_Shaw_Sun_pclt <- record("Fay_Shaw_Sun_pclt", ev$reject,
                        rep(horizon, length(p)), pvalue = ev$pvalue, score = ev$score, variance = ev$variance)
      if (max(abs(ev$score - U)) > 1e-8) stop("Fay-Shaw/common-grid Fixed score identity failed")
    }
    result <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    result[, scheme := scheme]
    result[, right_censored_fraction := mean(gen$data$event_time > gen$data$last_visit)]
    result
  }

  start <- proc.time()[["elapsed"]]
  for (n in args$n) for (scenario in args$scenarios) for (truth in c("null", "alternative")) {
    sp <- specs[[scenario]]
    for (bb in seq_len(ceiling(args$reps / args$batch))) {
      filename <- file.path(out, "batches", sprintf("n%d_%s_%s_b%03d.rds", n, scenario, truth, bb))
      if (args$resume && file.exists(filename)) next
      nr <- min(args$batch, args$reps - (bb - 1L) * args$batch)
      seed <- seed_for(n, scenario, truth, bb)
      gen <- generate(sp, n, truth, nr, seed)
      ans <- rbindlist(lapply(c("interval", "exact"), function(s) evaluate_scheme(gen, sp, s)))
      ans[, `:=`(n = n, scenario = scenario, truth = truth, batch = bb, seed = seed,
                  replicate = replicate + (bb - 1L) * args$batch, horizon = sp$K * sp$delta)]
      saveRDS(ans, paste0(filename, ".tmp"), compress = TRUE)
      if (!file.rename(paste0(filename, ".tmp"), filename)) stop("Could not commit batch")
      message(sprintf("[%s] n=%d %s %s batch %d/%d; elapsed %.1fs",
                  format(Sys.time(), "%H:%M:%S"), n, scenario, truth, bb,
                  ceiling(args$reps / args$batch), proc.time()[["elapsed"]] - start))
      rm(gen, ans); gc(FALSE)
    }
  }
  all <- rbindlist(lapply(list.files(file.path(out, "batches"), "\\.rds$", full.names = TRUE), readRDS))
  if (anyDuplicated(all[, .(n, scenario, truth, scheme, method, replicate)])) stop("Duplicate replicate outcomes")
  summary <- all[, {
    count <- sum(reject); ci <- wilson(count, .N)
    list(replicates = .N, rejections = count, rejection_probability = mean(reject),
      mc_se = sqrt(mean(reject) * (1 - mean(reject)) / .N), ci95_low = ci[[1L]], ci95_high = ci[[2L]],
      mean_stop_capped = mean(stop_time), stop_mc_se = sd(stop_time) / sqrt(.N),
      median_stop_given_reject = if (any(reject)) median(stop_time[reject]) else NA_real_,
      mean_terminal_loge = if (all(is.na(terminal_loge))) NA_real_ else mean(terminal_loge),
      mean_maximum_loge = if (all(is.na(maximum_loge))) NA_real_ else mean(maximum_loge),
      horizon = horizon[[1L]], mean_right_censored_fraction = mean(right_censored_fraction))
  }, by = .(n, scenario, truth, scheme, method)]
  if (any(summary$replicates != args$reps)) stop("Incomplete replicate counts")
  setorder(summary, n, scenario, truth, scheme, method)
  fwrite(summary, file.path(out, "summary.csv"))
  pairs <- list()
  for (reference in c("LR_fixed_original", "prequential_unrestricted", "prequential_protective", "fixed_score_LR")) {
    ref <- all[method == reference, .(n, scenario, truth, scheme, replicate,
                  ref_reject = as.numeric(reject), ref_stop = stop_time)]
    joined <- merge(all[method != reference], ref, by = c("n", "scenario", "truth", "scheme", "replicate"))
    z <- joined[, .(difference_rejection = mean(as.numeric(reject) - ref_reject),
      difference_rejection_mc_se = sd(as.numeric(reject) - ref_reject) / sqrt(.N),
      difference_capped_stop = mean(stop_time - ref_stop),
      difference_capped_stop_mc_se = sd(stop_time - ref_stop) / sqrt(.N),
      rejection_disagreements = sum(as.numeric(reject) != ref_reject), replicates = .N),
      by = .(n, scenario, truth, scheme, method)]
    z[, reference := reference]
    pairs[[reference]] <- z
  }
  fwrite(rbindlist(pairs), file.path(out, "paired_differences.csv"))
  fwrite(unique(all[, .(n, scenario, truth, batch, seed)]), file.path(out, "seed_ledger.csv"))
  sim1 <- merge(all[method == "Fay_Shaw_Sun_pclt"], all[method == "fixed_score_LR" & scheme == "interval"],
                by = c("n", "scenario", "truth", "scheme", "replicate"))
  identity <- sim1[, .(maximum_abs_score_difference = max(abs(score.x - score.y)),
                maximum_abs_variance_difference = max(abs(variance.x - variance.y)),
                rejection_disagreements = sum(reject.x != reject.y), replicates = .N), by = .(n, scenario, truth)]
  fwrite(identity, file.path(out, "sim1_score_equivalence.csv"))
  write_json(list(status = "complete", completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
     stochastic_datasets = args$reps * length(args$n) * length(args$scenarios) * 2L,
     evaluation_cells = length(args$n) * length(args$scenarios) * 2L * 2L,
     summary_rows = nrow(summary), outcome_rows = nrow(all),
     maximum_sim1_score_difference = max(identity$maximum_abs_score_difference),
     total_elapsed_seconds = proc.time()[["elapsed"]] - start,
     manuscript_modified = FALSE), file.path(out, "completion.json"), pretty = TRUE, auto_unbox = TRUE)
  message("Completed: ", out)
}

if (sys.nframe() == 0L) main()
