# Source-safe Fay--Shaw/Sun fixed-analysis interval-censored comparator.
# Primary method: interval::ictest(scores="logrank1", method="pclt").
# Source: https://doi.org/10.18637/jss.v036.i02
# Manual: https://cran.r-project.org/web/packages/interval/interval.pdf
# No simulations, installation, directory changes or manuscript writes on source().

sim1_pvalue <- function(z, alternative = c("less", "greater", "two.sided")) {
  alternative <- match.arg(alternative)
  switch(alternative, less = stats::pnorm(z),
         greater = stats::pnorm(z, lower.tail = FALSE),
         two.sided = 2 * stats::pnorm(-abs(z)))
}

# Matrices use the existing engine's units: one unit per planned assessment.
# For panel observation T in (k-1,k] is observed as precisely that interval;
# T>C is right censored at the last negative assessment C. No latent failure
# time is supplied to the interval-censored benchmark.
sim1_interval_records <- function(data, sc, replicate = 1L,
                                  scheme = c("interval", "exact"),
                                  visit_width = 1) {
  scheme <- match.arg(scheme)
  required <- c("event_time", "last_visit", "group")
  if (!all(required %in% names(data))) stop("Missing latent-data fields")
  if (length(replicate) != 1L || !is.finite(replicate) ||
      replicate != as.integer(replicate) || replicate < 1L) {
    stop("replicate must be a positive integer")
  }
  if (!is.matrix(data$group) || replicate > nrow(data$group) ||
      !identical(dim(data$group), dim(data$event_time)) ||
      !identical(dim(data$group), dim(data$last_visit))) {
    stop("Latent-data fields must be identically sized matrices")
  }
  if (length(visit_width) != 1L || !is.finite(visit_width) || visit_width <= 0) {
    stop("visit_width must be finite and positive")
  }
  t <- as.numeric(data$event_time[replicate, ])
  c <- pmin(as.numeric(data$last_visit[replicate, ]), sc$horizon)
  g <- as.integer(data$group[replicate, ])
  if (anyNA(t) || anyNA(c) || any(!is.finite(c)) || any(t < 0) ||
      any(c < 0) || any(!g %in% 0:1)) stop("Invalid latent records")
  observed <- t <= c
  L <- c
  R <- rep(Inf, length(t))
  if (scheme == "interval") {
    if (any(abs(c - round(c)) > 1e-10)) {
      stop("Panel censoring must occur on the engine's integer assessment grid")
    }
    # The generating distribution has positive times; explicitly reject T=0
    # rather than invent a negative lower endpoint for the first interval.
    if (any(t[observed] <= 0)) stop("Panel event times must be positive")
    R[observed] <- ceiling(t[observed])
    L[observed] <- R[observed] - 1
  } else {
    L[observed] <- R[observed] <- t[observed]
  }
  out <- data.frame(L = L * visit_width, R = R * visit_width, group = g)
  attr(out, "observation_scheme") <- scheme
  out
}

sim1_validate_records <- function(records) {
  if (!is.data.frame(records) || !all(c("L", "R", "group") %in% names(records))) {
    stop("records must have L, R and group columns")
  }
  L <- records$L
  R <- records$R
  g <- records$group
  if (!is.numeric(L) || !is.numeric(R) || length(L) < 2L ||
      anyNA(L) || anyNA(R) || anyNA(g) || any(!is.finite(L)) ||
      any(L < 0) || any(R < L) || any(!g %in% 0:1)) {
    stop("Invalid interval records; group must be numeric 0/1")
  }
  invisible(TRUE)
}

# A general records adapter, also usable off the common assessment grid.
# `less` always means fewer/longer event times in treatment group 1, matching
# Fixed's negative observed-minus-expected score. ictest follows first-occurring
# group order, so order treatment records first explicitly (factor levels alone
# do not control this in interval 1.1-1.0).
sim1_fay_shaw <- function(records,
                         alternative = c("less", "greater", "two.sided"),
                         scores = "logrank1", method = c("pclt", "scoretest"),
                         alpha = 0.05, icontrol = NULL, keep_fit = FALSE) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)
  sim1_validate_records(records)
  if (!requireNamespace("interval", quietly = TRUE)) {
    stop("Load the workspace-local .local/R-library containing interval")
  }
  n1 <- sum(records$group == 1L)
  n0 <- sum(records$group == 0L)
  base_diag <- list(package = as.character(utils::packageVersion("interval")),
                    scores = scores, method = method, n1 = n1, n0 = n0,
                    direction = "negative score means better treatment survival",
                    permutation_assumption = "exchangeable group labels under the null, including observation mechanism")
  if (!n1 || !n0 || all(is.infinite(records$R))) {
    return(list(score = 0, variance = 0, statistic = 0, pvalue = 1,
                reject = FALSE,
                diagnostics = c(base_diag, list(degenerate = TRUE, converged = NA,
                                                warnings = character()))))
  }
  ord <- c(which(records$group == 1L), which(records$group == 0L))
  z <- records[ord, ]
  warnings <- character()
  args <- list(L = z$L, R = z$R, group = z$group, scores = scores,
               alternative = alternative, method = method,
               mcontrol = interval::mControl(setSEED = FALSE))
  if (!is.null(icontrol)) args$icontrol <- icontrol
  fit <- withCallingHandlers(do.call(interval::ictest, args), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  score <- unname(fit$U[["1"]])
  if (method == "pclt") {
    variance <- n1 * n0 / ((n1 + n0) * (n1 + n0 - 1)) *
      sum((fit$scores - mean(fit$scores))^2)
  } else {
    variance <- if (is.null(fit$V)) 0 else as.numeric(fit$V[1L, 1L])
  }
  statistic <- if (variance > 0) score / sqrt(variance) else 0
  out <- list(score = score, variance = variance, statistic = statistic,
              pvalue = as.numeric(fit$p.value), reject = as.numeric(fit$p.value) <= alpha,
              diagnostics = c(base_diag, list(degenerate = variance <= 0,
                converged = isTRUE(fit$fit$converge),
                npmle_error = fit$fit$error, npmle_iterations = fit$fit$numit,
                package_statistic = unname(fit$statistic),
                package_method = fit$method, warnings = warnings,
                score_sum = sum(fit$scores))))
  if (keep_fit) out$fit <- fit
  out
}

# Exact common-grid reduction of the *permutation* form of the Sun logrank
# benchmark, not the likelihood scoretest variance. The common grid has every
# subject assessed at each integer visit until event or independent censoring;
# no delayed entry, stratification or mixed exact/panel observation is allowed.
# Input updates are the existing all-panel build(data,sc,"interval") output.
# n1/n0 are ENROLLED counts, needed for the permutation distribution, even when
# some participants leave before the first event-bearing cell.
#
# At event-bearing bin j, r_j subjects remain, d_j fail, h_j=d_j/r_j and
# H_j=sum_{l<=j}h_l. Sun individual scores are 1-H_j for failures in bin j
# and -H_j for censors after bin j. Thus U=sum_j(x_j-r1_j*h_j).
# The number censored after j and before the next event bin is
# r_j-d_j-r_{j+1}; set r_{J+1}=0. Censors before the first bin have score 0.
# Let Q=sum_j[d_j(1-H_j)^2 + (r_j-d_j-r_{j+1})H_j^2]. Then
# V_perm=n1*n0/[N(N-1)]*(Q-(sum_i c_i)^2/N).
# This is interval's permTS/pclt variance after rescaling the difference in
# group mean scores to U. It is generally not Fixed's hypergeometric V.
sim1_fast_updates <- function(updates, n1, n0,
                             alternative = c("less", "greater", "two.sided"),
                             alpha = 0.05, validate = TRUE) {
  alternative <- match.arg(alternative)
  required <- c("n1", "n0", "d", "x", "count", "time")
  if (!all(required %in% names(updates))) stop("Missing packed-update fields")
  nr <- nrow(updates$n1)
  if (is.null(nr)) stop("Packed updates must use matrices")
  if (length(n1) == 1L) n1 <- rep(n1, nr)
  if (length(n0) == 1L) n0 <- rep(n0, nr)
  if (length(n1) != nr || length(n0) != nr || anyNA(n1) || anyNA(n0) ||
      any(n1 < 0) || any(n0 < 0) || any(n1 != round(n1)) || any(n0 != round(n0))) {
    stop("n1 and n0 must be enrolled counts for every replicate")
  }
  U <- V <- Q <- total_score <- fixed_variance <- numeric(nr)
  for (i in seq_len(nr)) {
    k <- as.integer(updates$count[i])
    if (k <= 0L || n1[i] == 0 || n0[i] == 0) next
    idx <- seq_len(k)
    a <- as.numeric(updates$n1[i, idx])
    b <- as.numeric(updates$n0[i, idx])
    d <- as.numeric(updates$d[i, idx])
    x <- as.numeric(updates$x[i, idx])
    r <- a + b
    time <- updates$time[i, idx]
    cens <- r - d - c(r[-1L], 0)
    if (validate && (any(!is.finite(c(a,b,d,x,time))) ||
        any(abs(time - round(time)) > 1e-10) || any(diff(time) <= 0) ||
        any(d <= 0) || any(d > r) || any(cens < 0) ||
        any(x < 0) || any(x > pmin(d,a)) || any(d-x > b) ||
        a[1] > n1[i] || b[1] > n0[i])) {
      stop("Updates do not describe compatible common-grid all-panel records at replicate ", i)
    }
    H <- cumsum(d / r)
    U[i] <- sum(x - a * d / r)
    total_score[i] <- sum(d * (1-H) - cens * H)
    Q[i] <- sum(d * (1-H)^2 + cens * H^2)
    N <- n1[i] + n0[i]
    V[i] <- n1[i] * n0[i] / (N*(N-1)) *
      max(0, Q[i] - total_score[i]^2/N)
    informative <- r > 1
    fixed_variance[i] <- sum(d[informative] * a[informative] * b[informative] *
      (r[informative] - d[informative]) /
      (r[informative]^2 * (r[informative] - 1)))
  }
  Z <- numeric(nr)
  positive <- V > 0
  Z[positive] <- U[positive] / sqrt(V[positive])
  p <- sim1_pvalue(Z, alternative)
  p[!positive] <- 1
  list(score = U, variance = V, statistic = Z, pvalue = p,
       reject = p <= alpha,
       diagnostics = list(method = "Fay-Shaw/Sun logrank1 pclt common-grid reduction",
         scope = "all-panel common assessment grid; no delayed entry or strata",
         individual_score_square_sum = Q, pooled_score_sum = total_score,
         fixed_lr_variance = fixed_variance,
         variance_minus_fixed_lr = V - fixed_variance,
         n1 = n1, n0 = n0))
}
