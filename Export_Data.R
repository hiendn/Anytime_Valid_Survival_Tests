#!/usr/bin/env Rscript
# Export unmodified public example objects; no cohort selection or fitting.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- dirname(normalizePath(script, mustWork = TRUE))
out <- file.path(root, "data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
examples <- new.env()
data("msmdata", package = "msm", envir = examples)
for (name in c("cav", "psor")) {
  write.csv(get(name, envir = examples), file.path(out, paste0(name, ".csv")),
            row.names = FALSE, na = "")
}
cat("Public msm example data exported from package", as.character(packageVersion("msm")), "\n")
