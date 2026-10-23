#!/usr/bin/env Rscript
# Deterministic fixtures and small seeded validation samples, never production.
run_sim1_checks <- function(root = getwd(), output_dir = file.path(root, "results/session/sim1_validation")) {
  root <- normalizePath(root)
  .libPaths(c(file.path(root, ".local/R-library"), .libPaths()))
  source(file.path(root, "scripts/sim1_benchmark.R"), local = TRUE)
  battery <- new.env(parent = globalenv())
  source(file.path(root, "simulations/code/run_practical_battery.R"), local = battery)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  compare <- function(data, sc, label) {
    updates <- battery$build(data, sc, "interval")
    fast <- sim1_fast_updates(updates, rowSums(data$group), sc$n - rowSums(data$group))
    fixed <- battery$score_arrays(updates, "logrank")
    for (i in seq_len(nrow(data$group))) {
      records <- sim1_interval_records(data, sc, i)
      elapsed <- system.time(oracle <- sim1_fay_shaw(records,
        icontrol = interval::icfitControl(epsilon = 1e-10, maxit = 10000)))[["elapsed"]]
      rows[[length(rows)+1L]] <<- data.frame(
        case = label, n = sc$n, replicate = i,
        score_fast = fast$score[i], score_package = oracle$score,
        variance_fast = fast$variance[i], variance_package = oracle$variance,
        pvalue_fast = fast$pvalue[i], pvalue_package = oracle$pvalue,
        score_abs_diff = abs(fast$score[i] - oracle$score),
        variance_abs_diff = abs(fast$variance[i] - oracle$variance),
        pvalue_abs_diff = abs(fast$pvalue[i] - oracle$pvalue),
        reject_fast = fast$reject[i], reject_package = oracle$reject,
        fixed_lr_score = fixed$U[i,max(1L,updates$count[i])],
        fixed_lr_variance = fixed$V[i,max(1L,updates$count[i])],
        fixed_lr_reject = fixed$V[i,max(1L,updates$count[i])] > 0 &&
          fixed$Z[i,max(1L,updates$count[i])] <= qnorm(.05),
        package_converged = oracle$diagnostics$converged,
        warnings = paste(oracle$diagnostics$warnings, collapse = "; "),
        package_seconds = elapsed)
    }
    invisible(list(updates = updates, fast = fast))
  }
  fixture <- function(t, c, g) {
    sc <- battery$ENGINE$new_mixed_scenario("fixture", "fixture", "logrank",
      length(g), 8L, .21, .15, .15, .7, .12)
    data <- list(event_time = matrix(t,1), last_visit = matrix(as.integer(c),1),
      group = matrix(as.integer(g),1), regime_uniform = matrix(.5,1,length(g)),
      exit_order = matrix(order(pmin(t,c), seq_along(t)),1))
    compare(data, sc, "deterministic_fixture")
    list(data=data, sc=sc)
  }
  # T exactly on a visit stays in (k-1,k]; censoring at the same visit follows
  # event ascertainment. Includes a participant censored before the first event.
  f <- fixture(c(1,1.8,2.4,3,4.7,5.1,6.4,7.9,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
               c(8,8,8,8,8,8,8,8,0,1,2,3,4,5,6,8), rep(0:1,8))
  r <- sim1_interval_records(f$data, f$sc)
  stopifnot(r$L[1] == 0, r$R[1] == 1, r$L[4] == 2, r$R[4] == 3,
            r$L[9] == 0, is.infinite(r$R[9]))
  direction_records <- data.frame(L=c(rep(3,8),rep(0,8)),
                                  R=c(rep(4,8),rep(1,8)),
                                  group=c(rep(1L,8),rep(0L,8)))
  benefit <- sim1_fay_shaw(direction_records)
  harm <- sim1_fay_shaw(transform(direction_records,group=1L-group))
  stopifnot(benefit$score < 0, benefit$pvalue < .05,
            harm$score > 0, harm$pvalue > .95,
            abs(benefit$pvalue + harm$pvalue - 1) < 1e-10)
  # A heavily censored, balanced no-effect fixture checks finite p-values.
  fixture(c(1.2,1.2,3.3,3.3,Inf,Inf,Inf,Inf), rep(8,8), rep(0:1,4))
  no_event <- f$data
  no_event$event_time[,] <- Inf
  no_event$exit_order[1,] <- order(no_event$last_visit[1,], seq_len(f$sc$n))
  compare(no_event, f$sc, "no_events")
  one_group <- f$data
  one_group$group[,] <- 1L
  compare(one_group, f$sc, "one_group")
  # Same S1--S5 generating patterns used by the practical battery, under both
  # null and alternative, at every production sample size. Three reps per cell
  # are numerical implementation checks, not estimates of type-I error/power.
  # horizon, base rate (per assessment unit), cure0, cure1, hr1, hr2,
  # change time (assessment units), dropout probability.
  spec <- list(S1=c(8,.21,.15,.15,.65,.65,Inf,.14),
               S2=c(10,.25,.15,.34,1,1,Inf,.12),
               S3=c(12,.10,0,0,.65,1,5,1-exp(-.035)),
               S4=c(16,.09,0,0,1,.55,6,1-exp(-.03)),
               S5=c(10,.23,.20,.30,.82,.82,Inf,.035))
  index <- 0L
  for (n in c(100L,200L,500L,1000L,2000L)) for (s in names(spec)) {
    v <- spec[[s]]
    for (truth in c("null", "alternative")) {
      index <- index + 1L
      sc <- battery$ENGINE$new_mixed_scenario(s,s,"logrank",n,as.integer(v[1]),v[2],v[3],
         if(truth=="null")v[3] else v[4],1,v[8])
      data <- battery$custom_latent_piecewise(sc,3L,202610150L+index,
         treatment_probability=.5, hr_early=if(truth=="null")1 else v[5],
         hr_late=if(truth=="null")1 else v[6], change_time=v[7])
      compare(data,sc,paste(s,truth,sep="_"))
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out,file.path(output_dir,"package_equivalence_checks.csv"),row.names=FALSE)
  # Absolute tolerance accounts for the package NPMLE's numerical EM fit; the
  # algebraic fast reduction itself has no fitted nuisance parameters here.
  stopifnot(max(out$score_abs_diff) < 1e-5,
            max(out$variance_abs_diff) < 1e-5,
            max(out$pvalue_abs_diff) < 1e-6,
            all(out$reject_fast == out$reject_package),
            max(abs(out$score_fast-out$fixed_lr_score)) < 1e-10,
            all(is.finite(out$pvalue_fast)),
            !any(nzchar(out$warnings)))
  summary <- list(checks = nrow(out), all_pass = TRUE,
    direction_fixture_pass = TRUE,
    max_abs_score_difference = max(out$score_abs_diff),
    max_abs_variance_difference = max(out$variance_abs_diff),
    max_abs_pvalue_difference = max(out$pvalue_abs_diff),
    rejection_disagreements = sum(out$reject_fast != out$reject_package),
    fixed_lr_score_difference_max = max(abs(out$score_fast-out$fixed_lr_score)),
    fixed_lr_variance_difference_max = max(abs(out$variance_fast-out$fixed_lr_variance)),
    fixed_lr_rejection_disagreements = sum(out$reject_fast != out$fixed_lr_reject),
    package_seconds_total = sum(out$package_seconds),
    R = R.version.string,
    versions = setNames(vapply(c("interval","perm","Icens","MLEcens"),
      function(p) as.character(utils::packageVersion(p)), ""),
      c("interval","perm","Icens","MLEcens")))
  jsonlite::write_json(summary,file.path(output_dir,"summary.json"),auto_unbox=TRUE,pretty=TRUE)
  print(summary)
  invisible(summary)
}

if (sys.nframe() == 0L) run_sim1_checks()
