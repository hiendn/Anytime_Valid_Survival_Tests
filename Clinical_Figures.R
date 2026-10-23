#!/usr/bin/env Rscript
# Use the manuscript plotting functions without executing diagnostic calculations.
suppressPackageStartupMessages(library(data.table))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
ROOT <- dirname(normalizePath(script, mustWork = TRUE))
OUT <- file.path(ROOT, "build/figure_previews")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
archive <- fread(file.path(ROOT, "reference/liming_revision/data/application_weight_eprocess_paths.csv"))
codes <- c("LR", "GB", "PP", "FH", "GT")
parsed <- parse(file.path(ROOT, "scripts/referee_clinical_diagnostics.R"), keep.source = FALSE)
definition <- function(e, name) is.call(e) && identical(e[[1L]], as.name("<-")) &&
  identical(e[[2L]], as.name(name))
for (name in c("transform_y", "plot_one")) {
  selected <- which(vapply(parsed, definition, logical(1), name = name))
  stopifnot(length(selected) == 1L)
  eval(parsed[[selected]], envir = environment())
}
calls <- which(vapply(parsed, function(e)
  is.call(e) && identical(e[[1L]], as.name("plot_one")), logical(1)))
stopifnot(length(calls) == 2L)
for (i in calls) eval(parsed[[i]], envir = environment())
cat("Clinical figures recreated from saved evidence paths.\n")
