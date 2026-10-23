#!/usr/bin/env Rscript
# Replot saved, previously audited results. This script does not simulate data.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(getwd())
primary_path <- file.path(root, "simulations/results/expanded_primary.csv")
visit_path <- file.path(root, "simulations/results/visit_frequency.csv")
primary <- read.csv(primary_path, check.names = FALSE)
visit <- read.csv(visit_path, check.names = FALSE)
stopifnot(nrow(primary) == 160L, nrow(visit) == 120L,
          all(primary$n_reps == 5000L), all(visit$n_reps == 5000L))
out <- file.path(root, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
ns <- c(100, 200, 500, 1000, 2000)
cols <- c("#246A92", "#BC5C24", "#5F754B")
proc <- c("Anytime e-process", "Single terminal score test", "O'Brien--Fleming-type monitoring")
scenarios <- c("Latency proportional hazards", "Pure cure-rate difference")
draw_axes <- function(ylim, ylab) {
  plot(range(log10(ns)) + c(-0.04, 0.04), ylim, type = "n", xlab = "Sample size", ylab = ylab,
       xaxt = "n", xaxs = "i", yaxs = "i", bty = "l")
  axis(1, at = log10(ns), labels = ns)
  abline(h = pretty(ylim), col = "#EEEEEE", lwd = 0.7)
}
pdf(file.path(out, "referee_sim_sample_size.pdf"), width = 7.2, height = 6.4,
    family = "Times", useDingbats = FALSE)
layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.16))
par(mar = c(3.1, 3.6, 2.1, 0.9), oma = c(0, 0, 1.6, 0),
    mgp = c(2.2, 0.65, 0), tcl = -0.25, cex = 0.9)
for (s in seq_along(scenarios)) {
  draw_axes(c(0, 1.03), "Power")
  title(paste0("S", s, if (s == 1L) ": timing difference" else ": cure difference"), cex.main = 1)
  for (p in seq_along(proc)) for (scheme in c("interval", "exact")) {
    z <- primary[primary$scenario == scenarios[[s]] & primary$truth == "Alternative" &
                   primary$procedure == proc[[p]] & primary$scheme == scheme, ]
    z <- z[order(z$n), ]
    stopifnot(identical(as.numeric(z$n), ns))
    lines(log10(z$n), z$rejection_probability, col = cols[[p]],
          lty = if (scheme == "interval") 1 else 2, lwd = 2.8)
    points(log10(z$n), z$rejection_probability, col = cols[[p]],
           pch = if (scheme == "interval") 16 else 1, cex = 0.65)
  }
}
for (s in seq_along(scenarios)) {
  draw_axes(c(0, 0.062), "E-process null crossing probability")
  title(paste0("S", s, ": null"), cex.main = 1)
  abline(h = 0.05, col = "#777777", lty = 3)
  for (scheme in c("interval", "exact")) {
    z <- primary[primary$scenario == scenarios[[s]] & primary$truth == "Null" &
                   primary$procedure == proc[[1L]] & primary$scheme == scheme, ]
    z <- z[order(z$n), ]
    lines(log10(z$n), z$rejection_probability, col = cols[[1L]],
          lty = if (scheme == "interval") 1 else 2, lwd = 2.8)
    arrows(log10(z$n), pmax(0, z$rejection_probability - z$mc_se),
           log10(z$n), z$rejection_probability + z$mc_se,
           code = 3, angle = 90, length = 0.025, col = cols[[1L]], lwd = 1.1)
    points(log10(z$n), z$rejection_probability, col = cols[[1L]],
           pch = if (scheme == "interval") 16 else 1, cex = 0.65)
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", c("E-process", "Fixed", "OBF (1.96)", "Interval", "Exact"),
       col = c(cols, "#333333", "#333333"), lty = c(1, 1, 1, 1, 2),
       pch = c(NA, NA, NA, 16, 1), lwd = 2.8, horiz = TRUE, bty = "n", cex = 0.83)
mtext("OBF(1.96): illustrative conservative reference", side = 3, outer = TRUE,
      line = 0.35, cex = 0.85, font = 2)
dev.off()

pdf(file.path(out, "referee_sim_assessment_frequency.pdf"), width = 7.2, height = 3.5,
    family = "Times", useDingbats = FALSE)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.20))
par(mar = c(3.1, 3.6, 2.1, 0.9), oma = c(0, 0, 0, 0),
    mgp = c(2.2, 0.65, 0), tcl = -0.25, cex = 0.9)
widths <- c(0.5, 1, 2)
for (metric in c("rejection_probability", "median_stop_given_reject")) {
  power <- metric == "rejection_probability"
  draw_axes(if (power) c(0, 1.03) else c(0, 6.3), if (power) "Power" else "Years")
  title(if (power) "E-process power" else "Median stop among rejections", cex.main = 1)
  for (k in seq_along(widths)) {
    z <- visit[visit$truth == "Alternative" & visit$procedure == proc[[1L]] &
                 visit$visit_width == widths[[k]], ]
    z <- z[order(z$n), ]
    stopifnot(identical(as.numeric(z$n), ns))
    lines(log10(z$n), z[[metric]], col = cols[[k]], lty = k, lwd = 2.8)
    points(log10(z$n), z[[metric]], col = cols[[k]], pch = c(16, 17, 15)[[k]], cex = 0.65)
  }
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", c("6-month assessments", "Annual assessments", "2-year assessments"),
       col = cols, lty = 1:3, pch = c(16, 17, 15), lwd = 2.8, horiz = TRUE, bty = "n", cex = 0.9)
dev.off()
cat("Two figures regenerated from existing CSVs; no simulations run.\n")
