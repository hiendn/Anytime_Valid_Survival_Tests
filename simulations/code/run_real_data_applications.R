#!/usr/bin/env Rscript

# Deterministic native-R applications to the real msm::cav and msm::psor
# data. The output includes every fixed-margin table used by
# the e-process, complete monitored paths for six weight specifications,
# manuscript summaries, validation metadata, and a two-panel LR diagnostic.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) {
    normalizePath(gsub("~+~", " ", sub("^--file=", "", hit[[1L]]), fixed = TRUE), mustWork = TRUE)
  } else {
    normalizePath("code/run_real_data_applications.R", mustWork = TRUE)
  }
}

ROOT <- dirname(dirname(script_path()))
CODE <- file.path(ROOT, "code")
DEFAULT_RESULTS <- file.path(ROOT, "results")
DEFAULT_FIGURES <- file.path(ROOT, "manuscript", "figures")
source(file.path(CODE, "simulation_core.R"), local = TRUE)

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  value_after <- function(flag, default) {
    hit <- which(argv == flag)
    if (!length(hit)) return(default)
    if (
      length(hit) != 1L || hit == length(argv) ||
        startsWith(argv[[hit + 1L]], "--")
    ) {
      stop("Missing value after ", flag)
    }
    argv[[hit + 1L]]
  }
  valued <- c("--output-dir", "--figures-dir")
  positions <- which(argv %in% valued)
  consumed <- unique(c(positions, positions + 1L))
  bare <- if (length(consumed)) argv[-consumed] else argv
  unknown <- setdiff(bare, "--skip-figures")
  if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = ", "))
  list(
    output_dir = normalizePath(
      value_after("--output-dir", DEFAULT_RESULTS),
      winslash = "/", mustWork = FALSE
    ),
    figures_dir = normalizePath(
      value_after("--figures-dir", DEFAULT_FIGURES),
      winslash = "/", mustWork = FALSE
    ),
    skip_figures = "--skip-figures" %in% argv
  )
}

log_mean_exp <- function(z) {
  anchor <- max(z)
  anchor + log(mean(exp(z - anchor)))
}

eta_text <- function(eta) {
  paste(formatC(eta, digits = 12L, format = "fg", flag = "#"), collapse = ";")
}

mixture_path <- function(cells, eta, analysis_id, display_name,
                         weight = rep(1, nrow(cells))) {
  stopifnot(
    nrow(cells) > 0L, length(eta) > 0L,
    all(c("look_time", "n1", "n0", "d", "x") %in% names(cells)),
    length(weight) == nrow(cells), all(is.finite(weight)),
    all(weight >= 0)
  )
  component_increment <- vapply(
    eta,
    function(value) {
      as.numeric(table_loge(
        cells$n1, cells$n0, cells$d, cells$x,
        value * weight
      ))
    },
    numeric(nrow(cells))
  )
  if (is.null(dim(component_increment))) {
    component_increment <- matrix(component_increment, ncol = length(eta))
  }
  component_names <- paste0("component_", seq_along(eta))
  increment_dt <- as.data.table(component_increment)
  setnames(increment_dt, component_names)
  increment_dt[, look_time := cells$look_time]
  by_look <- increment_dt[
    , lapply(.SD, sum), by = look_time, .SDcols = component_names
  ]
  setorder(by_look, look_time)
  component_paths <- apply(
    as.matrix(by_look[, ..component_names]), 2L, cumsum
  )
  if (is.null(dim(component_paths))) {
    component_paths <- matrix(component_paths, ncol = length(eta))
  }
  log_e <- apply(component_paths, 1L, log_mean_exp)
  path <- data.table(
    analysis = analysis_id,
    display_name = display_name,
    look_index = seq_len(nrow(by_look)),
    look_time = by_look$look_time,
    e_value = exp(log_e),
    log_e_value = log_e
  )
  # The monitored process begins at E_0 = 1 before the first update.
  path[, max_e_to_date := cummax(pmax(1, e_value))]
  path[, crossed_20 := max_e_to_date >= 20]
  path[, eta_values := eta_text(eta)]
  path
}

# Each weight is fixed before the group-labelled event count of its cell is
# revealed. The pooled survival recursion uses prognostic bands, not the
# endpoint-specific fixed-margin cells. All cells at one release time use the
# same left-limit estimate within their band.
add_application_weights <- function(cells, baseline_band_counts,
                                    s_min = 0.12) {
  stopifnot(
    all(c("weight_band", "look_time", "n1", "n0", "d") %in% names(cells)),
    all(c("weight_band", "baseline_band_n") %in%
          names(baseline_band_counts)),
    !anyDuplicated(baseline_band_counts$weight_band), s_min > 0
  )
  cells <- copy(cells)
  cells[baseline_band_counts, on = "weight_band",
        baseline_band_n := i.baseline_band_n]
  if (anyNA(cells$baseline_band_n) ||
      any(cells$baseline_band_n <= 0L)) {
    stop("Missing or zero baseline prognostic-band count")
  }
  releases <- cells[, .(
    pooled_n = sum(n0 + n1), pooled_d = sum(d)
  ), by = .(weight_band, look_time)]
  if (any(releases$pooled_d > releases$pooled_n) ||
      any(releases$pooled_n <= 0L)) {
    stop("Invalid pooled release-time event fraction")
  }
  setorder(releases, weight_band, look_time)
  releases[, pooled_s_before := shift(
    cumprod(1 - pooled_d / pooled_n), fill = 1
  ), by = weight_band]
  cells[releases, on = .(weight_band, look_time),
        pooled_s_before := i.pooled_s_before]
  cells[, c("weight_LR", "weight_GB", "weight_PP", "weight_FH",
            "weight_GT") := list(
    1,
    (n0 + n1) / baseline_band_n,
    pooled_s_before,
    1 - pooled_s_before,
    1 / pmax(pooled_s_before, s_min)
  )]
  if (anyNA(cells$pooled_s_before) ||
      any(cells$weight_GB > 1 + 1e-12) ||
      any(!is.finite(as.matrix(cells[, .(
        weight_LR, weight_GB, weight_PP, weight_FH, weight_GT
      )])))) {
    stop("Invalid pre-allocation application weights")
  }
  cells
}

application_weight_paths <- function(result) {
  cells <- result$cells
  eta <- result$eta
  ordinary <- lapply(c("LR", "GB", "PP", "FH", "GT"), function(code) {
    path <- mixture_path(
      cells, eta, unique(cells$analysis), unique(cells$display_name),
      weight = cells[[paste0("weight_", code)]]
    )
    path[, weight_code := code]
    path
  })
  names(ordinary) <- c("LR", "GB", "PP", "FH", "GT")
  average <- copy(ordinary$LR)
  average[, e_value := (ordinary$LR$e_value + ordinary$GT$e_value) / 2]
  average[, c("log_e_value", "max_e_to_date", "crossed_20",
              "weight_code") := list(
    log(e_value), cummax(pmax(1, e_value)),
    cummax(pmax(1, e_value)) >= 20, "LR--GT"
  )]
  all_paths <- rbindlist(c(ordinary, list(average)))
  all_paths[, dataset := unique(cells$dataset)]
  summaries <- all_paths[, {
    peak <- which.max(e_value)
    crossing <- which(e_value >= 20)
    .(
      dataset = dataset[[1L]],
      is_primary = result$summary$is_primary,
      q = result$summary$q,
      time_unit = result$summary$time_unit,
      fixed_margin_cells = nrow(cells),
      looks = .N,
      max_e = e_value[[peak]],
      max_e_time = look_time[[peak]],
      terminal_e = e_value[[.N]],
      terminal_time = look_time[[.N]],
      first_crossing_20 = if (length(crossing)) {
        look_time[[crossing[[1L]]]]
      } else NA_real_,
      crossed_20 = length(crossing) > 0L,
      eta_values = eta_text(eta)
    )
  }, by = .(analysis, weight_code)]
  list(paths = all_paths, summary = summaries)
}

cav_intervals <- function(cav) {
  cav <- as.data.table(copy(cav))
  setorder(cav, PTNUM, years)
  baseline <- cav[, .SD[1L], by = PTNUM]
  eligible <- baseline[as.character(pdiag) %chin% c("IDC", "IHD")]
  if (any(eligible$state != 1L) || any(eligible$years != 0)) {
    stop("The eligible CAV baseline records are not all event-free at transplant")
  }

  interval_rows <- vector("list", nrow(eligible) * 4L)
  interval_count <- 0L
  subject_rows <- vector("list", nrow(eligible))
  for (ii in seq_len(nrow(eligible))) {
    id <- eligible$PTNUM[[ii]]
    records <- cav[PTNUM == id]
    event_type <- "censored"
    completed <- 0L
    if (nrow(records) >= 2L) {
      for (kk in 2:nrow(records)) {
        current_state <- records$state[[kk]]
        if (current_state == 4L) {
          event_type <- "death_before_CAV"
          break
        }
        if (!current_state %in% 1:3) {
          stop("Unexpected CAV state for patient ", id)
        }
        completed <- completed + 1L
        interval_count <- interval_count + 1L
        if (interval_count > length(interval_rows)) {
          length(interval_rows) <- length(interval_rows) * 2L
        }
        interval_rows[[interval_count]] <- data.table(
          id = id,
          group = as.integer(records$pdiag[[1L]] == "IHD"),
          group_label = as.character(records$pdiag[[1L]]),
          left = records$years[[kk - 1L]],
          right = records$years[[kk]],
          event = as.integer(current_state %in% c(2L, 3L)),
          sex = records$sex[[1L]],
          recipient_age = records$age[[1L]],
          donor_age = records$dage[[1L]]
        )
        if (current_state %in% c(2L, 3L)) {
          event_type <- "CAV"
          break
        }
      }
    }
    subject_rows[[ii]] <- data.table(
      id = id,
      group = as.integer(eligible$pdiag[[ii]] == "IHD"),
      group_label = as.character(eligible$pdiag[[ii]]),
      outcome = event_type,
      completed_intervals = completed
    )
  }
  intervals <- rbindlist(interval_rows[seq_len(interval_count)])
  subjects <- rbindlist(subject_rows)
  if (any(intervals$right <= intervals$left)) stop("Nonpositive raw CAV interval")
  if (anyDuplicated(intervals[event == 1L, id])) stop("A CAV event was counted twice")
  list(
    baseline = eligible,
    subjects = subjects,
    intervals = intervals,
    excluded = nrow(baseline) - nrow(eligible)
  )
}

cav_analysis <- function(intervals, subjects, baseline, q, primary = FALSE) {
  x <- copy(intervals)
  x[, c(
    "rounded_left", "rounded_right", "recipient_age_band",
    "donor_age_band", "sex_band"
  ) := list(
    q * round(left / q),
    q * round(right / q),
    cut(
      recipient_age, c(-Inf, 45, 55, Inf), right = FALSE,
      labels = c("<45", "45--54", "55+")
    ),
    cut(
      donor_age, c(-Inf, 30, Inf), right = FALSE,
      labels = c("<30", "30+")
    ),
    fifelse(sex == 1L, "female", "male")
  )]
  omitted <- x[rounded_right <= rounded_left]
  if (nrow(omitted) && any(omitted$event != 0L)) {
    stop("A first CAV detection would be omitted by the zero-width rule")
  }
  x <- x[rounded_right > rounded_left]
  if (anyDuplicated(x[, .(id, rounded_right)])) {
    stop("A subject appears twice in one rounded CAV release batch")
  }
  cells <- x[, .(
    n1 = sum(group),
    n0 = .N - sum(group),
    d = sum(event),
    x = sum(group * event),
    interval_records = .N,
    actual_right_max = max(right)
  ), by = .(
    rounded_left, rounded_right, sex_band,
    recipient_age_band, donor_age_band
  )]
  setorder(
    cells, rounded_right, rounded_left, sex_band,
    recipient_age_band, donor_age_band
  )
  analysis_id <- sprintf(
    "cav_q%.2f%s", q, if (primary) "_primary" else ""
  )
  display_name <- sprintf(
    "CAV, q = %.2f y%s", q, if (primary) " (primary)" else ""
  )
  cells[, c(
    "dataset", "analysis", "display_name", "observation_scheme",
    "look_time", "stratum", "mixture", "q", "q_unit"
  ) := list(
    "cav",
    analysis_id,
    display_name,
    "interval",
    rounded_right + q / 2,
    paste(
      sprintf("(%.2f,%.2f]", rounded_left, rounded_right),
      sex_band, recipient_age_band, donor_age_band, sep = " | "
    ),
    "equal two-sided logrank mixture: eta = +/-log(1.5), +/-log(2)",
    q,
    "years"
  )]
  # The rounded endpoint names the working comparison cell. Its fixed upper
  # bin edge is a predeclared inspection time after every assessment assigned
  # to that bin, unlike the retrospectively known last observed assessment.
  if (any(cells$look_time < cells$actual_right_max - 1e-10)) {
    stop("A CAV assessment lies beyond its planned release time")
  }
  cells[, update_id := seq_len(.N)]
  cells[, weight_band := paste(
    sex_band, recipient_age_band, donor_age_band, sep = " | "
  )]
  baseline <- as.data.table(copy(baseline))
  baseline[, weight_band := paste(
    fifelse(sex == 1L, "female", "male"),
    cut(
      age, c(-Inf, 45, 55, Inf), right = FALSE,
      labels = c("<45", "45--54", "55+")
    ),
    cut(dage, c(-Inf, 30, Inf), right = FALSE,
        labels = c("<30", "30+")),
    sep = " | "
  )]
  baseline_band_counts <- baseline[, .(
    baseline_band_n = .N
  ), by = weight_band]
  cells <- add_application_weights(cells, baseline_band_counts)
  eta <- c(-log(2), -log(1.5), log(1.5), log(2))
  path <- mixture_path(cells, eta, analysis_id, display_name)
  path[, dataset := "cav"]
  crossing <- which(path$e_value >= 20)
  max_row <- which.max(path$e_value)
  group_counts <- subjects[, .N, by = group]
  n0_total <- group_counts[group == 0L, N]
  n1_total <- group_counts[group == 1L, N]
  outcome_counts <- subjects[, .N, by = outcome]
  event_group_counts <- subjects[outcome == "CAV", .N, by = group]
  count_outcome <- function(name) {
    value <- outcome_counts[outcome == name, N]
    if (length(value)) value else 0L
  }
  summary <- data.table(
    dataset = "cav",
    analysis = analysis_id,
    display_name = display_name,
    is_primary = primary,
    observation_scheme = "Interval-censored first CAV detection",
    comparison = "IHD versus IDC (IHD is group 1)",
    q = q,
    q_unit = "years",
    time_unit = "years",
    n0 = n0_total,
    n1 = n1_total,
    n_total = n0_total + n1_total,
    events = sum(cells$d),
    events_n0 = event_group_counts[group == 0L, N],
    events_n1 = event_group_counts[group == 1L, N],
    censorings = count_outcome("censored") + count_outcome("death_before_CAV"),
    censorings_n0 = subjects[outcome != "CAV" & group == 0L, .N],
    censorings_n1 = subjects[outcome != "CAV" & group == 1L, .N],
    death_before_event = count_outcome("death_before_CAV"),
    raw_interval_records = nrow(intervals),
    analysed_interval_records = nrow(x),
    omitted_zero_width_negative = nrow(omitted),
    fixed_margin_cells = nrow(cells),
    informative_event_cells = cells[d > 0L & n1 > 0L & n0 > 0L, .N],
    informative_events = cells[n1 > 0L & n0 > 0L, sum(d)],
    looks = nrow(path),
    mixture = unique(cells$mixture),
    max_e = path$e_value[[max_row]],
    max_e_time = path$look_time[[max_row]],
    terminal_e = path$e_value[[nrow(path)]],
    terminal_time = path$look_time[[nrow(path)]],
    first_crossing_20 = if (length(crossing)) {
      path$look_time[[crossing[[1L]]]]
    } else {
      NA_real_
    },
    crossed_20 = length(crossing) > 0L,
    secondary_hr = NA_real_,
    secondary_p = NA_real_
  )
  list(cells = cells, path = path, summary = summary, eta = eta)
}

psor_cohort <- function(psor) {
  psor <- as.data.table(copy(psor))
  setorder(psor, ptnum, months)
  if (psor[, any(diff(months) <= 0), by = ptnum][V1 == TRUE, .N]) {
    stop("The psor visit times are not strictly increasing within patient")
  }
  if (psor[, any(diff(state) < 0), by = ptnum][V1 == TRUE, .N]) {
    stop("The psor damage states are not monotone within patient")
  }

  baseline <- psor[, .SD[1L], by = ptnum]
  eligible <- baseline[state == 1L & !is.na(ollwsdrt)]
  eligible[, c("group", "baseline_effusion", "baseline_months") := list(
    as.integer(ollwsdrt == 1L), as.integer(hieffusn == 1L), months
  )]
  subject_map <- eligible[, .(
    id = ptnum, group, baseline_effusion, baseline_months
  )]
  records <- merge(
    psor[ptnum %in% eligible$ptnum], subject_map,
    by.x = "ptnum", by.y = "id", all.x = TRUE, sort = FALSE
  )
  setorder(records, ptnum, months)
  records[, aligned_months := months - baseline_months]
  records <- records[aligned_months > 0]
  if (records[, uniqueN(ptnum)] != nrow(eligible)) {
    stop("Every eligible psor patient must have a post-baseline assessment")
  }

  raw_intervals <- records[, {
    first_positive <- which(state > 1L)
    last_index <- if (length(first_positive)) first_positive[[1L]] else .N
    out <- copy(.SD[seq_len(last_index)])
    out[, event := 0L]
    if (length(first_positive)) out[.N, event := 1L]
    out[, left_actual := shift(aligned_months, fill = 0)]
    out
  }, by = ptnum]
  setnames(raw_intervals, "ptnum", "id")
  if (any(raw_intervals$aligned_months <= raw_intervals$left_actual)) {
    stop("A nonpositive psor visit interval was constructed")
  }
  if (anyDuplicated(raw_intervals[event == 1L, id])) {
    stop("A psoriatic-arthritis progression was counted twice")
  }

  subjects <- raw_intervals[, .(
    group = group[[1L]],
    baseline_effusion = baseline_effusion[[1L]],
    outcome = if (any(event == 1L)) "progression" else "censored",
    completed_intervals = .N
  ), by = id]
  list(
    source = psor,
    baseline = baseline,
    eligible = eligible,
    records = records,
    raw_intervals = raw_intervals,
    subjects = subjects,
    excluded_initial_damage = baseline[state != 1L, .N],
    excluded_missing_baseline_esr = baseline[state == 1L & is.na(ollwsdrt), .N]
  )
}

psor_analysis <- function(records, raw_intervals, subjects, baseline,
                          q, primary = FALSE) {
  x <- copy(records)
  x[, release_bin := pmax(1L, as.integer(ceiling(aligned_months / q)))]
  setorder(x, ptnum, release_bin, aligned_months)
  # The last visit in a nominal window is selected by time alone, before the
  # recorded state is used to define first progression.
  x <- x[, .SD[.N], by = .(ptnum, release_bin)]
  x <- x[, {
    first_positive <- which(state > 1L)
    last_index <- if (length(first_positive)) first_positive[[1L]] else .N
    out <- copy(.SD[seq_len(last_index)])
    out[, event := 0L]
    if (length(first_positive)) out[.N, event := 1L]
    out
  }, by = ptnum]
  setnames(x, "ptnum", "id")
  setorder(x, id, release_bin)
  x[, c("previous_release_bin", "previous_actual") := list(
    shift(release_bin, fill = 0L),
    shift(aligned_months, fill = 0)
  ), by = id]
  x[, c("rounded_left", "rounded_right", "effusion_band") := list(
    q * previous_release_bin,
    q * release_bin,
    fifelse(baseline_effusion == 1L, "high effusion", "not high effusion")
  )]
  if (anyDuplicated(x[, .(id, rounded_right)])) {
    stop("A psor patient appears twice in one nominal release batch")
  }
  if (sum(x$event) != sum(raw_intervals$event)) {
    stop("Nominal psor windows did not preserve every first progression")
  }

  cells <- x[, .(
    n1 = sum(group),
    n0 = .N - sum(group),
    d = sum(event),
    x = sum(group * event),
    interval_records = .N,
    actual_left_min = min(previous_actual),
    actual_left_max = max(previous_actual),
    actual_right_min = min(aligned_months),
    actual_right_max = max(aligned_months)
  ), by = .(rounded_left, rounded_right, baseline_effusion, effusion_band)]
  setorder(cells, rounded_right, rounded_left, baseline_effusion)
  analysis_id <- sprintf(
    "psor_q%d%s", as.integer(q), if (primary) "_primary" else ""
  )
  display_name <- sprintf(
    "PsA progression, q = %d months%s",
    as.integer(q), if (primary) " (primary)" else ""
  )
  cells[, c(
    "dataset", "analysis", "display_name", "observation_scheme",
    "look_time", "stratum", "mixture", "q", "q_unit",
    "sex_band", "recipient_age_band", "donor_age_band"
  ) := list(
    "psor",
    analysis_id,
    display_name,
    "nominal interval",
    rounded_right,
    paste(
      sprintf("(%d,%d] months", rounded_left, rounded_right),
      effusion_band, sep = " | "
    ),
    "equal protective logrank mixture: conditional odds ratios = 0.50, 0.67, 0.80",
    q,
    "months",
    NA_character_,
    NA_character_,
    NA_character_
  )]
  cells[, update_id := seq_len(.N)]
  cells[, weight_band := effusion_band]
  baseline <- as.data.table(copy(baseline))
  baseline[, weight_band := fifelse(
    baseline_effusion == 1L, "high effusion", "not high effusion"
  )]
  baseline_band_counts <- baseline[, .(
    baseline_band_n = .N
  ), by = weight_band]
  cells <- add_application_weights(cells, baseline_band_counts)
  eta <- log(c(0.50, 0.67, 0.80))
  path <- mixture_path(cells, eta, analysis_id, display_name)
  path[, dataset := "psor"]
  crossing <- which(path$e_value >= 20)
  max_row <- which.max(path$e_value)
  group_counts <- subjects[, .N, by = group]
  event_group_counts <- subjects[outcome == "progression", .N, by = group]
  censor_group_counts <- subjects[outcome == "censored", .N, by = group]
  summary <- data.table(
    dataset = "psor",
    analysis = analysis_id,
    display_name = display_name,
    is_primary = primary,
    observation_scheme = "Interval-censored first joint-damage progression",
    comparison = "Low versus not-low baseline ESR (low is group 1)",
    q = q,
    q_unit = "months",
    time_unit = "months",
    n0 = group_counts[group == 0L, N],
    n1 = group_counts[group == 1L, N],
    n_total = nrow(subjects),
    events = sum(x$event),
    events_n0 = event_group_counts[group == 0L, N],
    events_n1 = event_group_counts[group == 1L, N],
    censorings = sum(subjects$outcome == "censored"),
    censorings_n0 = censor_group_counts[group == 0L, N],
    censorings_n1 = censor_group_counts[group == 1L, N],
    death_before_event = NA_integer_,
    raw_interval_records = nrow(raw_intervals),
    analysed_interval_records = nrow(x),
    collapsed_within_window = nrow(raw_intervals) - nrow(x),
    omitted_zero_width_negative = 0L,
    fixed_margin_cells = nrow(cells),
    informative_event_cells = cells[d > 0L & n1 > 0L & n0 > 0L, .N],
    informative_events = cells[n1 > 0L & n0 > 0L, sum(d)],
    looks = nrow(path),
    mixture = unique(cells$mixture),
    max_e = path$e_value[[max_row]],
    max_e_time = path$look_time[[max_row]],
    terminal_e = path$e_value[[nrow(path)]],
    terminal_time = path$look_time[[nrow(path)]],
    first_crossing_20 = if (length(crossing)) {
      path$look_time[[crossing[[1L]]]]
    } else {
      NA_real_
    },
    crossed_20 = length(crossing) > 0L,
    secondary_hr = NA_real_,
    secondary_p = NA_real_
  )
  list(data = x, cells = cells, path = path, summary = summary, eta = eta)
}

application_figure <- function(paths, output_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required unless --skip-figures is used")
  }
  cav <- copy(paths[analysis == "cav_q0.25_primary"])
  psor <- copy(paths[analysis == "psor_q6_primary"])
  cav[, panel := "(a) CAV: IHD versus IDC"]
  psor[, panel := "(b) PsA: low versus not-low baseline ESR"]
  cav_initial <- copy(cav[1L])
  cav_initial[, c(
    "look_index", "look_time", "e_value", "log_e_value",
    "max_e_to_date", "crossed_20"
  ) := list(0L, 0, 1, 0, 1, FALSE)]
  psor_initial <- copy(psor[1L])
  psor_initial[, c(
    "look_index", "look_time", "e_value", "log_e_value",
    "max_e_to_date", "crossed_20"
  ) := list(0L, 0, 1, 0, 1, FALSE)]
  plot_data <- rbindlist(list(cav_initial, cav, psor_initial, psor))
  threshold_data <- unique(plot_data[, .(panel)])
  crossing <- cav[e_value >= 20][1L]

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = look_time, y = e_value)
  ) +
    ggplot2::geom_hline(
      data = threshold_data,
      ggplot2::aes(yintercept = 20),
      colour = "#A33A2B", linewidth = 0.55, linetype = "longdash"
    ) +
    ggplot2::geom_hline(
      data = threshold_data,
      ggplot2::aes(yintercept = 1),
      colour = "grey55", linewidth = 0.35, linetype = "dotted"
    ) +
    ggplot2::geom_step(
      colour = "#173B57", linewidth = 0.65, direction = "hv"
    ) +
    ggplot2::geom_point(
      data = crossing,
      ggplot2::aes(x = look_time, y = e_value),
      inherit.aes = FALSE,
      shape = 21, size = 2.5, stroke = 0.55,
      fill = "white", colour = "#A33A2B"
    ) +
    ggplot2::facet_wrap(~panel, scales = "free_x", nrow = 1L) +
    ggplot2::scale_y_log10(
      limits = c(0.03, 600),
      breaks = c(0.05, 0.1, 0.5, 1, 5, 20, 100, 500),
      labels = c("0.05", "0.1", "0.5", "1", "5", "20", "100", "500")
    ) +
    ggplot2::labs(
      x = "Aligned follow-up (years for CAV; months for PsA)",
      y = "Mixture e-value (log scale)"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey94", colour = "grey65"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(5, 6, 4, 5),
      axis.title = ggplot2::element_text(size = 9)
    )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(output_dir, "real_data_eprocess.pdf"),
    p, width = 7.2, height = 3.35, units = "in",
    device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(output_dir, "real_data_eprocess.png"),
    p, width = 7.2, height = 3.35, units = "in",
    dpi = 300, bg = "white"
  )
}

main <- function() {
  args <- parse_args()
  required <- "msm"
  missing <- required[
    !vapply(required, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  data_env <- new.env(parent = emptyenv())
  utils::data("msmdata", package = "msm", envir = data_env)
  if (!all(c("cav", "psor") %in% ls(data_env))) {
    stop("The msmdata archive did not provide both cav and psor")
  }
  cav_raw <- data_env$cav
  cav <- cav_intervals(cav_raw)
  cav_results <- list(
    cav_analysis(cav$intervals, cav$subjects, cav$baseline,
                 0.10, primary = FALSE),
    cav_analysis(cav$intervals, cav$subjects, cav$baseline,
                 0.25, primary = TRUE),
    cav_analysis(cav$intervals, cav$subjects, cav$baseline,
                 0.50, primary = FALSE)
  )
  psor_raw <- data_env$psor
  psor <- psor_cohort(psor_raw)
  psor_results <- list(
    psor_analysis(
      psor$records, psor$raw_intervals, psor$subjects, psor$eligible,
      3, primary = FALSE
    ),
    psor_analysis(
      psor$records, psor$raw_intervals, psor$subjects, psor$eligible,
      6, primary = TRUE
    ),
    psor_analysis(
      psor$records, psor$raw_intervals, psor$subjects, psor$eligible,
      12, primary = FALSE
    )
  )

  summaries <- rbindlist(
    c(
      lapply(cav_results, function(z) z[["summary"]]),
      lapply(psor_results, function(z) z[["summary"]])
    ),
    fill = TRUE
  )
  cells <- rbindlist(
    c(
      lapply(cav_results, function(z) z[["cells"]]),
      lapply(psor_results, function(z) z[["cells"]])
    ),
    fill = TRUE
  )
  paths <- rbindlist(
    c(
      lapply(cav_results, function(z) z[["path"]]),
      lapply(psor_results, function(z) z[["path"]])
    ),
    fill = TRUE
  )
  weighted <- lapply(c(cav_results, psor_results), application_weight_paths)
  weight_summaries <- rbindlist(lapply(weighted, `[[`, "summary"))
  weight_paths <- rbindlist(lapply(weighted, `[[`, "paths"))
  # The manuscript reports every weight at every assessment width. Check the
  # Cartesian product, not only its size, so a duplicate cannot hide an omission.
  expected_weight_grid <- CJ(
    analysis=summaries$analysis,
    weight_code=c("LR","GB","PP","FH","GT","LR--GT")
  )
  observed_weight_grid <- weight_summaries[, .(analysis,weight_code)]
  weight_grid_complete <- !anyDuplicated(observed_weight_grid) &&
    fsetequal(observed_weight_grid,expected_weight_grid)
  path_counts <- weight_paths[, .(path_looks=.N),by=.(analysis,weight_code)]
  path_checks <- merge(weight_summaries,path_counts,
    by=c("analysis","weight_code"),all=TRUE)
  weight_paths_complete <- nrow(path_checks)==nrow(expected_weight_grid) &&
    !anyNA(path_checks[,.(looks,path_looks)]) &&
    all(path_checks$looks==path_checks$path_looks)
  if(!weight_grid_complete || !weight_paths_complete)
    stop("Incomplete application weight-by-window results")
  setcolorder(
    cells,
    c(
      "dataset", "analysis", "display_name", "observation_scheme",
      "update_id", "look_time", "stratum", "n1", "n0", "d", "x",
      "mixture", "q", "q_unit", "rounded_left", "rounded_right",
      "baseline_effusion", "effusion_band",
      "sex_band", "recipient_age_band", "donor_age_band",
      "interval_records", "actual_left_min", "actual_left_max",
      "actual_right_min", "actual_right_max", "weight_band",
      "baseline_band_n", "pooled_s_before", "weight_LR", "weight_GB",
      "weight_PP", "weight_FH", "weight_GT"
    )
  )
  fwrite(
    summaries, file.path(args$output_dir, "application_summary.csv"), na = ""
  )
  fwrite(
    cells,
    file.path(args$output_dir, "application_fixed_margin_updates.csv"),
    na = ""
  )
  fwrite(
    paths, file.path(args$output_dir, "application_eprocess_paths.csv"), na = ""
  )
  fwrite(
    weight_summaries,
    file.path(args$output_dir, "application_weight_summary.csv"), na = ""
  )
  fwrite(
    weight_paths,
    file.path(args$output_dir, "application_weight_eprocess_paths.csv"),
    na = ""
  )

  checks <- list(
    cav_baseline_total = nrow(cav$baseline) == 583L,
    cav_group_counts = identical(
      as.integer(cav$subjects[order(group), .N, by = group]$N),
      c(270L, 313L)
    ),
    cav_event_total = sum(cav$intervals$event) == 214L,
    cav_events_counted_once = !anyDuplicated(cav$intervals[event == 1L, id]),
    cav_primary_margin_bounds = with(
      cav_results[[2L]]$cells,
      all(x >= pmax(0L, d - n0) & x <= pmin(n1, d))
    ),
    cav_zero_width_rows_are_negative =
      cav_results[[2L]]$summary$omitted_zero_width_negative == 2L,
    cav_all_q_sensitivities_cross =
      all(summaries[dataset == "cav", crossed_20]),
    psor_source_dimensions = identical(dim(psor_raw), c(806L, 5L)),
    psor_source_subjects = uniqueN(psor_raw$ptnum) == 305L,
    psor_initial_state_one = psor$baseline[state == 1L, .N] == 212L,
    psor_missing_baseline_esr = psor$excluded_missing_baseline_esr == 27L,
    psor_eligible_total = nrow(psor$subjects) == 185L,
    psor_group_counts = identical(
      as.integer(psor$subjects[order(group), .N, by = group]$N),
      c(120L, 65L)
    ),
    psor_progressions = sum(psor$raw_intervals$event) == 72L,
    psor_events_counted_once =
      !anyDuplicated(psor$raw_intervals[event == 1L, id]),
    psor_primary_records =
      psor_results[[2L]]$summary$analysed_interval_records == 193L,
    psor_primary_margin_bounds = with(
      psor_results[[2L]]$cells,
      all(x >= pmax(0L, d - n0) & x <= pmin(n1, d))
    ),
    psor_all_windows_preserve_events =
      all(summaries[dataset == "psor", events] == 72L),
    psor_all_windows_do_not_cross =
      all(!summaries[dataset == "psor", crossed_20]),
    all_e_values_finite_positive =
      all(is.finite(paths$e_value) & paths$e_value > 0),
    all_weighted_e_values_finite_positive =
      all(is.finite(weight_paths$e_value) & weight_paths$e_value > 0),
    all_weight_window_combinations_present = weight_grid_complete,
    all_weight_window_paths_complete = weight_paths_complete,
    weighted_lr_reproduces_original = isTRUE(all.equal(
      weight_paths[weight_code == "LR", e_value],
      paths$e_value, tolerance = 1e-12
    )),
    arithmetic_lr_gt_mixture = isTRUE(all.equal(
      weight_paths[weight_code == "LR--GT", e_value],
      (weight_paths[weight_code == "LR", e_value] +
       weight_paths[weight_code == "GT", e_value]) / 2,
      tolerance = 1e-12
    ))
  )
  outcome_counts <- cav$subjects[, .N, by = outcome]
  validation <- list(
    analysis = "Real-data applications",
    deterministic = TRUE,
    r_version = R.version.string,
    rng_kind = as.list(RNGkind()),
    packages = list(
      msm = as.character(utils::packageVersion("msm")),
      data_table = as.character(utils::packageVersion("data.table")),
      Rcpp = as.character(utils::packageVersion("Rcpp"))
    ),
    sources = list(
      cav = "msm data archive msmdata; Sharples et al. (2003)",
      psor = "msm data archive msmdata; Gladman and Farewell (1999)"
    ),
    weight_specifications = list(
      comparison_status = "descriptive retrospective sensitivity analyses",
      number_of_window_weight_combinations = nrow(weight_summaries),
      assessment_widths = list(cav_years=as.list(c(0.10,0.25,0.50)),
                               psor_months=as.list(c(3,6,12))),
      coefficient_rule = "retain the original dataset-specific eta mixture for every weight",
      pooled_survival_rule = paste(
        "Within each prognostic band, all cells at a release use the",
        "left-limit pooled survival; after release multiply it by",
        "1 - sum(events)/sum(completed interval records)."
      ),
      gb_rule = "cell completed-interval count / eligible baseline count in its prognostic band",
      gt_survival_floor = 0.12,
      lr_gt_rule = "arithmetic average of the complete LR and GT e-processes"
    ),
    cav = list(
      source_rows = nrow(cav_raw),
      source_subjects = length(unique(cav_raw$PTNUM)),
      eligible_subjects = nrow(cav$subjects),
      excluded_other_or_missing_diagnosis = cav$excluded,
      completed_interval_records = nrow(cav$intervals),
      first_CAV_detections = sum(cav$intervals$event),
      outcome_counts = as.list(setNames(outcome_counts$N, outcome_counts$outcome)),
      primary_window_q_years = 0.25,
      application_weight_codes = as.list(c(
        "LR", "GB", "PP", "FH", "GT", "LR--GT"
      )),
      pooled_release_survival_floor = 0.12,
      primary_strata = c(
        "rounded interval endpoints", "sex",
        "baseline recipient age <45/45--54/55+",
        "donor age <30/30+"
      ),
      mixture_eta = as.list(c(-log(2), -log(1.5), log(1.5), log(2)))
    ),
    psor = list(
      source_rows = nrow(psor_raw),
      source_subjects = uniqueN(psor_raw$ptnum),
      subjects_initially_state_one = psor$baseline[state == 1L, .N],
      excluded_initial_damage = psor$excluded_initial_damage,
      excluded_missing_baseline_esr = psor$excluded_missing_baseline_esr,
      eligible_subjects = nrow(psor$subjects),
      completed_raw_intervals = nrow(psor$raw_intervals),
      first_progressions = sum(psor$raw_intervals$event),
      primary_window_q_months = 6,
      application_weight_codes = as.list(c(
        "LR", "GB", "PP", "FH", "GT", "LR--GT"
      )),
      pooled_release_survival_floor = 0.12,
      primary_strata = c(
        "previous nominal release", "current nominal release",
        "baseline high-effusion status"
      ),
      mixture_conditional_odds_ratios = as.list(c(0.50, 0.67, 0.80))
    ),
    checks = checks,
    all_checks_pass = all(unlist(checks, use.names = FALSE)),
    session_info = capture.output(utils::sessionInfo())
  )
  write_json(
    validation,
    file.path(args$output_dir, "application_validation.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"
  )
  if (!args$skip_figures) application_figure(paths, args$figures_dir)

  print(summaries[, .(
    analysis, n_total, events, fixed_margin_cells, looks,
    max_e, max_e_time, terminal_e, first_crossing_20
  )])
  if (!validation$all_checks_pass) stop("Application validation checks failed")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main()
