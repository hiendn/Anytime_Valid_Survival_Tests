#!/usr/bin/env Rscript
# Recomputable structural/statistical audit of the full mixed-resolution output.
validate_mixed_results <- function(output_dir, require_production = TRUE) {
  get_csv <- function(name) utils::read.csv(file.path(output_dir, name), check.names = FALSE)
  summary <- get_csv('mixed_resolution_summary.csv')
  calibration <- get_csv('mixed_resolution_eta_calibration.csv')
  paired <- get_csv('mixed_resolution_paired_differences.csv')
  observation <- get_csv('mixed_resolution_observation_mix.csv')
  manifest <- jsonlite::fromJSON(file.path(output_dir, 'mixed_resolution_validation.json'))
  settings <- manifest$settings
  ns <- as.integer(settings$sample_sizes)
  scenarios <- trimws(strsplit(settings$scenarios, ',', fixed = TRUE)[[1L]])
  fractions <- as.numeric(strsplit(settings$fractions, ',', fixed = TRUE)[[1L]])
  methods <- c('mixed', 'all_panel', 'exact_only', 'all_exact_oracle')
  procedures <- c('anytime_e', 'fixed_classical', 'naive_repeated')
  stopifnot(isTRUE(manifest$algebra_within_tolerance), isTRUE(manifest$all_endpoint_identities_hold))
  if (require_production) stopifnot(
    identical(ns, c(100L,200L,500L,1000L,2000L)), settings$pilot_reps == 1000L,
    settings$eval_reps == 5000L, settings$eta_grid == 50L,
    identical(fractions, c(0,0.25,0.5,0.75,1)), length(scenarios) == 2L,
    isTRUE(settings$production), !isTRUE(settings$smoke))
  key <- function(data, fields) do.call(paste, c(data[fields], sep = '|'))
  expected <- expand.grid(scenario = scenarios, n = ns, truth = c('H0','H1'),
                          exact_fraction = fractions, method = methods, procedure = procedures,
                          stringsAsFactors = FALSE)
  fields <- names(expected)
  stopifnot(!anyDuplicated(key(summary,fields)),
            setequal(key(summary,fields),key(expected,fields)),
            all(summary$n_reps == settings$eval_reps),
            all(summary$rejection_probability >= 0 & summary$rejection_probability <= 1),
            all(is.finite(summary$mean_terminal_metric)), all(is.finite(summary$mean_updates)),
            all(summary$mean_updates >= 0),
            max(abs(summary$mc_se - sqrt(summary$rejection_probability *
                   (1-summary$rejection_probability)/summary$n_reps))) < 1e-12)
  compare <- function(rho, a, b) {
    x <- summary[summary$exact_fraction == rho & summary$method == a,]
    y <- summary[summary$exact_fraction == rho & summary$method == b,]
    fields <- c('scenario','n','truth','procedure')
    i <- match(key(x,fields),key(y,fields))
    metrics <- c('rejection_probability','mc_se','median_stop_given_reject',
                 'mean_stop_capped','mean_terminal_metric','mean_updates')
    stopifnot(isTRUE(all.equal(unname(as.matrix(x[metrics])),
                               unname(as.matrix(y[i,metrics])), tolerance = 1e-12)))
  }
  if (0 %in% fractions) {
    compare(0,'mixed','all_panel')
    empty <- summary[summary$exact_fraction == 0 & summary$method == 'exact_only',]
    stopifnot(all(empty$rejection_probability == 0), all(empty$mean_updates == 0),
              all(empty$mean_terminal_metric == 0))
  }
  if (1 %in% fractions) {
    compare(1,'mixed','all_exact_oracle')
    compare(1,'exact_only','all_exact_oracle')
  }
  # Panel and oracle are paired reference analyses of the same data at each rho.
  for (method in c('all_panel','all_exact_oracle')) {
    a <- summary[summary$method == method,]
    groups <- split(a, key(a,c('scenario','n','truth','procedure')))
    stopifnot(all(vapply(groups, function(g) length(unique(g$rejection_probability)) == 1L &&
                                      length(unique(g$mean_terminal_metric)) == 1L, logical(1))))
  }
  stopifnot(all(calibration$n_reps == settings$pilot_reps),
            !anyDuplicated(key(calibration,c('scenario','n','exact_fraction','method','eta'))),
            setequal(calibration$n,ns))
  curve_groups <- split(calibration,key(calibration,c('scenario','n','exact_fraction','method')))
  stopifnot(all(vapply(curve_groups, nrow, integer(1)) == settings$eta_grid))
  for (row in seq_len(nrow(summary))) {
    x <- summary[row,]
    if (x$method == 'exact_only' && x$exact_fraction == 0) next
    method <- x$method
    rho <- x$exact_fraction
    if (method == 'mixed' && rho == 0) method <- 'all_panel'
    if (method %in% c('mixed','exact_only') && rho == 1) method <- 'all_exact_oracle'
    if (method %in% c('all_panel','all_exact_oracle')) rho <- -1
    curve <- calibration[calibration$scenario == x$scenario & calibration$n == x$n &
                          calibration$method == method & calibration$exact_fraction == rho,]
    stopifnot(nrow(curve) == settings$eta_grid,
              abs(x$eta - curve$eta[which.max(curve$mean_terminal_loge)]) < 1e-12)
  }
  obs_expected <- expand.grid(scenario=scenarios,n=ns,truth=c('H0','H1'),exact_fraction=fractions)
  fields <- names(obs_expected)
  components <- c('exact_event_fraction','left_censored_fraction',
                  'interval_censored_fraction','right_censored_fraction')
  stopifnot(!anyDuplicated(key(observation,fields)),
            setequal(key(observation,fields), key(obs_expected,fields)),
            all(observation$n_reps == settings$eval_reps),
            max(abs(rowSums(observation[components])-1)) < 1e-12)
  alternative <- summary[summary$truth == 'H1' & summary$procedure == 'anytime_e',]
  stopifnot(nrow(paired) == length(ns)*length(scenarios)*length(fractions)*3L,
            all(paired$n_reps == settings$eval_reps), all(paired$mc_se >= 0))
  for (row in seq_len(nrow(paired))) {
    p <- paired[row,]
    selected <- alternative[alternative$scenario == p$scenario & alternative$n == p$n &
                             alternative$exact_fraction == p$exact_fraction,]
    difference <- selected$rejection_probability[selected$method==p$method_a] -
                  selected$rejection_probability[selected$method==p$method_b]
    stopifnot(length(difference)==1L, abs(difference - p$difference_a_minus_b)<1e-12,
              abs(p$ci_lower - (difference-qnorm(.5+p$confidence_level/2)*p$mc_se))<1e-12,
              abs(p$ci_upper - (difference+qnorm(.5+p$confidence_level/2)*p$mc_se))<1e-12)
  }
  seeds <- unlist(manifest$seed_schedule$scenario,use.names=FALSE)
  stopifnot(!anyDuplicated(seeds))
  saved_paths <- readRDS(file.path(output_dir,'mixed_resolution_example_paths.rds'))
  if (.5 %in% fractions) stopifnot(length(saved_paths) == length(ns)*length(scenarios))
  null <- summary[summary$truth == 'H0' & summary$procedure == 'anytime_e',]
  list(passed=TRUE, sample_sizes=ns, summary_rows=nrow(summary), calibration_rows=nrow(calibration),
       paired_rows=nrow(paired), observation_rows=nrow(observation),
       null_crossing_range=range(null$rejection_probability),
       unique_stream_seeds=length(seeds))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  script <- gsub("~+~", " ", sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1L]), fixed = TRUE)
  root <- dirname(dirname(normalizePath(script)))
  directory <- if(length(args)) args[[1L]] else file.path(root,'results')
  print(validate_mixed_results(directory,require_production=!'--allow-smoke' %in% args))
}
