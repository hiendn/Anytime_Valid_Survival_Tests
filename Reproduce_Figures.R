#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- dirname(normalizePath(script, mustWork = TRUE))
setwd(root)
source(file.path(root, "scripts/referee_sim_figures.R"), local = new.env(parent = globalenv()))
source(file.path(root, "Clinical_Figures.R"), local = new.env(parent = globalenv()))
