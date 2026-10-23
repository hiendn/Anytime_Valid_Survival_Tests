# Source-safe figures for the common total-sample-size grid. Each row is one N
# and each column one data-generating scenario; no cells pool sample sizes.
make_mixed_resolution_grid_figures <- function(summary, observation, example_paths, figures_dir) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  ns <- sort(unique(summary$n))
  scenario_titles <- c(latency_ph = "Latency proportional hazards", pure_cure = "Pure cure-rate difference")
  scenarios <- intersect(names(scenario_titles), unique(summary$scenario))
  open_grid <- function(filename) {
    path <- file.path(figures_dir, filename)
    grDevices::pdf(path, width = 11, height = 3.15 * length(ns) + 0.45, onefile = TRUE)
    graphics::par(mfrow = c(length(ns), length(scenarios)),
                  mar = c(3.9, 4.5, 2.1, 1.2), oma = c(1.8, 0, 0, 0),
                  cex = 0.92, cex.main = 1.1)
    path
  }
  close_methods <- function() {
    .add_device_legend(unname(METHOD_LABELS[METHODS]), col = .method_colours[METHODS], lty = 1, pch = 16)
    grDevices::dev.off()
  }
  title <- function(scenario, n) paste0(scenario_titles[[scenario]], "; N = ", n)
  alternative <- summary[summary$procedure == "anytime_e" & summary$truth == "H1", , drop = FALSE]
  null <- summary[summary$procedure == "anytime_e" & summary$truth == "H0", , drop = FALSE]
  paths <- character()
  for (definition in list(
    list(file = "mixed_resolution_power.pdf", data = alternative,
         value = "rejection_probability", label = "Finite-horizon anytime power", limits = c(0, 1)),
    list(file = "mixed_resolution_loge.pdf", data = alternative,
         value = "mean_terminal_metric", label = "Mean terminal log e", horizontal = LOG_BOUNDARY),
    list(file = "mixed_resolution_null.pdf", data = null,
         value = "rejection_probability", label = "Null crossing probability", limits = c(0, max(0.08, 1.05 * max(null$rejection_probability))), horizontal = ALPHA)
  )) {
    paths <- c(paths, open_grid(definition$file))
    for (n in ns) {
      row_values <- definition$data[definition$data$n == n, definition$value]
      limits <- definition$limits
      if (is.null(limits)) {
        limits <- range(c(row_values, definition$horizontal), finite = TRUE)
        pad <- max(diff(limits) * 0.08, 0.25)
        limits <- limits + c(-pad, pad)
      }
      for (scenario in scenarios) {
        selected <- definition$data[definition$data$scenario == scenario & definition$data$n == n, , drop = FALSE]
        .draw_method_lines(selected, definition$value, title(scenario, n), definition$label,
                           limits, horizontal = definition$horizontal, show_y = TRUE)
      }
    }
    close_methods()
  }
  paths <- c(paths, open_grid("mixed_resolution_observation_mix.pdf"))
  components <- c("exact_event_fraction", "left_censored_fraction", "interval_censored_fraction", "right_censored_fraction")
  labels <- c("Exact event", "Left-censored", "Interval-censored", "Right-censored")
  colours <- c("#2878B5", "#F4A259", "#6FB98F", "#9B8AC4")
  for (n in ns) for (scenario in scenarios) {
    selected <- observation[observation$n == n & observation$scenario == scenario & observation$truth == "H1", , drop = FALSE]
    selected <- selected[order(selected$exact_fraction), , drop = FALSE]
    graphics::barplot(t(as.matrix(selected[, components, drop = FALSE])),
                      names.arg = format(selected$exact_fraction, trim = TRUE),
                      ylim = c(0, 1), col = colours, border = NA,
                      xlab = "Fraction continuously monitored", ylab = "Proportion of subjects",
                      main = title(scenario, n))
  }
  .add_device_legend(labels, col = NULL, fill = colours)
  grDevices::dev.off()
  paths <- c(paths, open_grid("mixed_resolution_example_paths.pdf"))
  for (n in ns) for (scenario in scenarios) {
    one <- example_paths[[paste(scenario, n, sep = "|")]]
    all_times <- unlist(lapply(one, function(path) path$time), use.names = FALSE)
    all_loge <- unlist(lapply(one, function(path) path$loge), use.names = FALSE)
    xlimits <- if (length(all_times)) range(c(0, all_times), finite = TRUE) else c(0, 1)
    ylimits <- range(c(0, LOG_BOUNDARY, all_loge), finite = TRUE)
    pad <- max(diff(ylimits) * 0.08, 0.25)
    graphics::plot(NA_real_, NA_real_, xlim = xlimits, ylim = ylimits + c(-pad, pad),
                   xlab = "Calendar time (50% exact)", ylab = "log e-value", main = title(scenario, n))
    graphics::grid(col = "grey88")
    graphics::abline(h = LOG_BOUNDARY, lty = 2, col = "grey35")
    for (method in METHODS) if (length(one[[method]]$time)) {
      graphics::lines(one[[method]]$time, one[[method]]$loge, type = "s", lwd = 1.4,
                     col = .method_colours[[method]])
    }
  }
  close_methods()
  invisible(paths)
}
