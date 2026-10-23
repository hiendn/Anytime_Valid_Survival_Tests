# Common design, shared by every active simulation driver.
# n is the TOTAL number of enrolled participants, not per treatment arm.
SAMPLE_SIZES <- c(100L, 200L, 500L, 1000L, 2000L)
REFERENCE_N <- 500L
PILOT_REPS <- 1000L
EVALUATION_REPS <- 5000L
ETA_GRID_SIZE <- 50L

validate_sample_sizes <- function(x) {
  x <- as.integer(x)
  if (!length(x) || anyNA(x) || any(x < 2L) || anyDuplicated(x))
    stop("Sample sizes must be distinct positive integers >= 2")
  sort(x)
}
