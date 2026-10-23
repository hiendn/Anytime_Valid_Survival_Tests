#!/usr/bin/env Rscript
# Create a report and scientific figures from a completed isolated run.
args <- commandArgs(TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/report_remaining_simulations.R RUN_DIRECTORY")
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(ggplot2)})
run <- normalizePath(args[[1L]], mustWork = TRUE)
done <- fromJSON(file.path(run, "completion.json"))
if (!identical(done$status, "complete")) stop("Run is incomplete")
s <- fread(file.path(run, "summary.csv"))
p <- fread(file.path(run, "paired_differences.csv"))
eq <- fread(file.path(run, "sim1_score_equivalence.csv"))
scales <- fread(file.path(run, "frozen_scales.csv"))
design <- fromJSON(file.path(run, "design.json"))
if (file.exists(file.path(run, "design_full_precision.json")))
  design <- fromJSON(file.path(run, "design_full_precision.json"))
format_value <- function(x, digits = 3L) ifelse(is.na(x), "—", formatC(x, format = "f", digits = digits))
table_md <- function(x, digits = 3L) {
  x <- copy(x)
  for (nm in names(x)) if (is.numeric(x[[nm]])) x[[nm]] <- format_value(x[[nm]], digits)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |")))
}
labels <- c(LR_fixed_original = "LR fixed/original fixed", GB_fixed = "GB fixed",
  PP_fixed = "PP fixed", FH_fixed = "FH fixed", GT_fixed = "GT fixed",
  LR_coefficient_mixture = "LR coefficient mixture", LR_FH_GT_mixture = "LR/FH/GT mixture",
  prequential_unrestricted = "Prequential (grouped extension for intervals)",
  prequential_protective = "Protective prequential adaptation",
  fixed_score_LR = "Fixed LR terminal", Fay_Shaw_Sun_pclt = "Fay–Shaw/Sun terminal")
null <- s[truth == "null"]
primary <- c("LR_fixed_original", "LR_coefficient_mixture", "LR_FH_GT_mixture",
             "prequential_unrestricted", "prequential_protective")
e_null <- null[!method %in% c("fixed_score_LR", "Fay_Shaw_Sun_pclt")]
delayed <- s[n == 500 & scenario == "S4" & scheme == "interval" & truth == "alternative"]
value <- function(method0, field) delayed[method == method0][[field]][[1L]]
gain <- p[n == 500 & scenario == "S4" & scheme == "interval" & truth == "alternative" &
           method == "FH_fixed" & reference == "LR_fixed_original"]
text <- c("# Remaining simulation TODOs: completed experiment", "",
  sprintf("Mode: **%s**. %s independent null and %s alternative datasets per scenario/sample-size cell, shared by all procedures and both observation maps.", design$mode, design$reps_per_truth_cell, design$reps_per_truth_cell), "",
  sprintf("%s independent latent datasets; %s evaluation cells across S1–S5, the selected sample sizes and exact/interval observations. No clinical preprocessing was rerun and the manuscript remains unchanged.", done$stochastic_datasets, done$evaluation_cells), "",
  "## Findings", "",
  sprintf("- SIM-1: the generalized logrank benchmark and Fixed LR scores coincide, but their variances differ. Their rejection decisions differed on %s of %s datasets (%.2f%%).", sum(eq$rejection_disagreements), sum(eq$replicates), 100*sum(eq$rejection_disagreements)/sum(eq$replicates)),
  sprintf("- SIM-2, delayed effect S4 at n=500 with interval observation: fixed FH power was %.4f versus %.4f for fixed LR; the LR/FH/GT mixture had power %.4f. The FH–LR paired power difference was %.4f (MCSE %.4f).", value("FH_fixed","rejection_probability"), value("LR_fixed_original","rejection_probability"), value("LR_FH_GT_mixture","rejection_probability"), gain$difference_rejection, gain$difference_rejection_mc_se),
  sprintf("- In that delayed-effect case, FH's mean capped stopping time was %.4f years versus %.4f for LR. Unrestricted grouped prequential power was %.4f; the protective adaptation had power %.4f. These comparisons apply to the frozen designs tested.", value("FH_fixed","mean_stop_capped"), value("LR_fixed_original","mean_stop_capped"), value("prequential_unrestricted","rejection_probability"), value("prequential_protective","rejection_probability")),
  sprintf("- All e-process null rejection estimates were at or below %.4f across %s cells; the asymptotic terminal procedures are audited separately. This empirical check supports calibration in these designs, without establishing a new guarantee.", max(e_null$rejection_probability), nrow(e_null)),
  "- The fixed LR/FH/GT weight mixture improved delayed-effect power over fixed LR in this example but lost power in several other patterns. The LR coefficient mixture did not improve this delayed-effect case. There is no single best specification across all the tested patterns.", "",
  "## Design frozen before evaluation", "",
  "The five data-generating mechanisms, geometric censoring, horizons and assessment spacings match the manuscript. S3/S4 stopping times are converted from half-year engine units to calendar years. Coefficients and mixture probabilities do not change across scenarios or sample sizes.", "",
  sprintf("All single-weight methods use eta=log(0.7)=%.8f. Weight scales come from %s independent S1 null pilots at n=500 per observation scheme, then stay fixed. Pooled survival is evaluated before each allocation; GT uses a survival floor of 0.12.", design$fixed_eta, design$pilot_reps), "",
  "The LR coefficient mixture averages three complete processes at odds ratios 0.5, 0.7 and 0.9. The LR/FH/GT mixture averages nine complete processes with equal probability. Alternative-specific pilot tuning is not used in this experiment.", "",
  table_md(scales[, .(scheme, weight, scale, normalization_multiplier)], 5L), "",
  "## SIM-1: established interval-censored benchmark", "",
  "The benchmark is Fay and Shaw's interval implementation of Sun's generalized logrank score, with `scores=\"logrank1\"`, `method=\"pclt\"`, and a protective one-sided alternative. This permutation asymptotic calibration assumes exchangeable group labels under the null, including the observation mechanism, as in these simulations.", "",
  "On the common assessment grid its score equals the unweighted Fixed LR score. Its permutation variance generally differs from Fixed's accumulated hypergeometric variance, so the calibrated procedures need not reject the same datasets. A common-grid reduction was checked directly against the official R package before production; that reduction is used for the full run.", "",
  sprintf("Across the full evaluation, the largest absolute score discrepancy was %.3g; the largest absolute variance difference was %.5g. There were %s paired rejection disagreements across %s interval-record datasets.", max(eq$maximum_abs_score_difference), max(eq$maximum_abs_variance_difference), sum(eq$rejection_disagreements), sum(eq$replicates)), "")
for (truth0 in c("null", "alternative")) {
  x <- dcast(s[n == 500 & scheme == "interval" & truth == truth0 & method %in% c("fixed_score_LR", "Fay_Shaw_Sun_pclt")], scenario ~ method, value.var = "rejection_probability")
  text <- c(text, paste0("At n=500, ", if (truth0 == "null") "null rejection" else "terminal power", ":"), "", table_md(x), "")
}
text <- c(text, "## SIM-2: fixed practical designs and prequential comparisons", "",
  "The fixed unit-weight row is algebraically the original fixed-alternative procedure at the same coefficient and comparison sets. The implementation checked its terminal evidence, rejection decisions and stopping times against our LR row on every simulated dataset; it is not reported twice. Interval comparisons use the manuscript's completed-window risk sets, which differ from the original paper's events-before-censoring convention. The interval coefficient specifies conditional event odds and is not automatically a population hazard ratio.", "",
  "The unrestricted prequential comparator uses the published two-virtual-event smoothed MLE from strictly earlier release times and the initial enrolled group counts. For grouped event counts, fitting the grouped Fisher noncentral hypergeometric likelihood is an explicit extension. The protective comparator additionally caps the fitted log odds at zero; this direction-constrained version is an adaptation, not the published original procedure. All cells sharing a release time are predicted before their allocations enter the fit.", "")
text <- c(text, "At n=500 with interval observations, single-weight powers were:", "",
  table_md(dcast(s[n==500 & scheme=="interval" & truth=="alternative" & method %in%
    c("LR_fixed_original","GB_fixed","PP_fixed","FH_fixed","GT_fixed")], scenario~method,value.var="rejection_probability")), "")
for (scheme0 in c("interval", "exact")) {
  x <- dcast(s[n == 500 & scheme == scheme0 & truth == "alternative" & method %in% primary], scenario ~ method, value.var = "rejection_probability")
  text <- c(text, paste0("### n=500, ", scheme0, " observations: power"), "", table_md(x), "")
  x <- dcast(s[n == 500 & scheme == scheme0 & truth == "alternative" & method %in% primary], scenario ~ method, value.var = "mean_stop_capped")
  text <- c(text, "Mean capped stopping time (years; noncrossers receive the horizon):", "", table_md(x), "")
}
text <- c(text, "## Monte Carlo uncertainty and interpretation", "",
  sprintf("Maximum Bernoulli Monte Carlo SE over all rows: %.5f. Full results include 95%% Wilson intervals, which remain meaningful at rejection estimates zero or one. Paired rejection and capped-stop differences use within-dataset differences and retain their Monte Carlo SEs.", max(s$mc_se)), "",
  sprintf("Null rejection estimates range from %.4f to %.4f across the reported methods and cells. Individual estimates above 0.05 are expected from Monte Carlo variation across many cells; the CSV reports their intervals rather than treating each excursion as proof of invalidity.", min(null$rejection_probability), max(null$rejection_probability)), "",
  "No procedure should be selected retrospectively from these results and described as a prespecified design. The fixed configurations examine robustness across patterns; they do not replace the earlier scenario-specific pilot comparisons, which remain targeted calibrations. Comparisons of different weights include their joint coefficient/scale specification.", "",
  "## Reproducibility files", "",
  "- `design.json`, `design_full_precision.json`, `frozen_scales.csv` and `.rds`: frozen settings, source hashes and full-precision parameters. The original JSON displays rounded numbers; evaluation used exact expressions in the frozen source.",
  "- `summary.csv`: all rejection probabilities, Wilson intervals, MCSEs and stopping summaries.",
  "- `paired_differences.csv`: paired power and capped-stop differences with MCSEs.",
  "- `sim1_score_equivalence.csv`: score identity and variance/rejection differences.",
  "- `batches/*.rds`, `seed_ledger.csv`: replicate-level outcomes, independent batch seeds and pairing.",
  "- `session-info.txt`, `completion.json`: environment and completion record.",
  "- `qa_report.json`, `qa_null_calibration.csv/json`: independent output, replay and calibration audit.",
  "- `recovery_record.json`: terminal-interruption recovery after 943 saved batches, reusing the same frozen design and seeds. Completion timing measures the resumed stage.", "",
  "## Primary references", "",
  "- Fay MP and Shaw PA (2010), *Exact and Asymptotic Weighted Logrank Tests for Interval Censored Data: The interval R Package*. https://doi.org/10.18637/jss.v036.i02", 
  "- Official interval manual: https://cran.r-project.org/web/packages/interval/interval.pdf", 
  "- ter Schure et al. (2024), *The Anytime-Valid Logrank Test: Error Control under Continuous Monitoring with Unlimited Horizon*, §§3.2–3.3 and equations (3.7)–(3.8). https://nejsds.nestat.org/journal/NEJSDS/article/73/text", "")
writeLines(text, file.path(run, "REPORT.md"))

plot_data <- s[truth == "alternative" & method %in% primary]
plot_data[, label := factor(labels[method], levels = labels[primary])]
theme_set(theme_minimal(base_size = 10) + theme(
  plot.background = element_rect(fill = "white", colour = NA),
  panel.background = element_rect(fill = "white", colour = NA),
  panel.grid.major = element_line(colour = "#D9D9D9"), panel.grid.minor = element_blank(),
  text = element_text(colour = "#222222"), axis.text = element_text(colour = "#333333")))
palette <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00")
g <- ggplot(plot_data, aes(n, rejection_probability, color = label)) +
  geom_line(linewidth = .6) + geom_point(size = 1) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high), width = 0, alpha = .3) +
  facet_grid(scheme ~ scenario) + scale_x_log10(breaks = sort(unique(s$n))) +
  scale_y_continuous(limits = c(0,1)) + labs(x = "Total enrolled sample size", y = "Power by horizon", color = NULL) +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(file.path(run, "power_by_scenario.pdf"), g, width = 13, height = 6.5, bg = "white")
ggsave(file.path(run, "power_by_scenario.png"), g, width = 13, height = 6.5, dpi = 160, bg = "white")
g <- ggplot(plot_data, aes(n, mean_stop_capped / horizon, color = label)) +
  geom_line(linewidth = .6) + geom_point(size = 1) + facet_grid(scheme ~ scenario) +
  scale_x_log10(breaks = sort(unique(s$n))) +
  labs(x = "Total enrolled sample size", y = "Mean capped stop / horizon", color = NULL) +
  scale_color_manual(values = palette) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(file.path(run, "capped_stop_by_scenario.pdf"), g, width = 13, height = 6.5, bg = "white")
ggsave(file.path(run, "capped_stop_by_scenario.png"), g, width = 13, height = 6.5, dpi = 160, bg = "white")
single <- c("LR_fixed_original","GB_fixed","PP_fixed","FH_fixed","GT_fixed")
plot_data <- s[truth == "alternative" & method %in% single]
plot_data[, label := factor(labels[method], levels = labels[single])]
g <- ggplot(plot_data, aes(n, rejection_probability, color = label)) +
  geom_line(linewidth = .6) + geom_point(size = 1) +
  geom_errorbar(aes(ymin = ci95_low, ymax = ci95_high), width = 0, alpha = .3) +
  facet_grid(scheme ~ scenario) + scale_x_log10(breaks = sort(unique(s$n))) +
  scale_y_continuous(limits = c(0,1)) + scale_color_manual(values = palette) +
  labs(x = "Total enrolled sample size", y = "Power by horizon", color = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))
ggsave(file.path(run, "single_weights_power.pdf"), g, width = 13, height = 6.5, bg = "white")
ggsave(file.path(run, "single_weights_power.png"), g, width = 13, height = 6.5, dpi = 160, bg = "white")
cat("Report and figures saved in", run, "\n")
