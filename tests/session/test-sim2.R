#!/usr/bin/env Rscript
# Deterministic tests only: no Monte Carlo data generation or driver execution.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
root <- dirname(dirname(dirname(normalizePath(script, mustWork = TRUE))))
source(file.path(root, "simulations", "code", "simulation_core.R"))
source(file.path(root, "scripts", "sim2_comparators.R"))
Sys.setenv(SIMULATION_RCPP_CACHE_DIR = file.path(root, "build", "rcpp_cache"))
library(testthat)

packed <- function(n1,n0,d,x,time,count=length(x)) {
  list(n1=matrix(as.integer(n1),1), n0=matrix(as.integer(n0),1),
       d=matrix(as.integer(d),1), x=matrix(as.integer(x),1),
       time=matrix(as.numeric(time),1), count=as.integer(count))
}
logz <- function(eta,b,a,d) {
  y <- seq.int(max(0L,d-a),min(b,d))
  terms <- lchoose(b,y)+lchoose(a,d-y)+eta*y
  max(terms)+log(sum(exp(terms-max(terms))))
}
independent_eta <- function(u, prior_cells, a=10,b=10) {
  f <- function(eta) {
    # The literal paper smoothing, evaluated stably in log odds.
    base <- eta-log(exp(log(a+1))+exp(eta+log(b+1)))-log(exp(log(a))+exp(eta+log(b+1)))
    if (length(prior_cells)) for (j in prior_cells)
      base <- base+eta*u$x[j]-logz(eta,u$n1[j],u$n0[j],u$d[j])
    base
  }
  optimize(f,c(-15,15),maximum=TRUE,tol=1e-11)$maximum
}

test_that("fixed unit-weight comparator matches the existing exact table kernel", {
  u <- packed(c(10,9,8,5,4),c(10,8,7,5,3),c(1,2,4,3,2),c(1,1,2,1,1),1:5)
  out <- sim2_comparators(u,10,10,5,keep_path=TRUE)
  inc <- table_loge(u$n1,u$n0,u$d,u$x,matrix(log(.7),1,5))
  expect_equal(as.numeric(out$fixed$loge_path),cumsum(as.numeric(inc)),tolerance=1e-11)
  expect_equal(out$fixed$terminal_loge,sum(inc),tolerance=1e-11)
})

test_that("smoothing and all prequential estimates match independent likelihood optimization", {
  u <- packed(c(10,10,9,8,6),c(10,9,8,7,5),c(1,2,2,3,2),c(0,1,1,2,0),1:5)
  out <- sim2_comparators(u,10,10,5,keep_path=TRUE)
  eta <- as.numeric(out$prequential$eta_path)
  expect_equal(eta[1],log(sqrt(10*11)/11),tolerance=1e-12)
  for (k in seq_along(eta))
    expect_equal(eta[k],independent_eta(u,seq_len(k-1)),tolerance=1e-6)
  expect_equal(as.numeric(out$prequential_protective$eta_path),pmin(eta,0),tolerance=1e-12)
})

test_that("current and later allocations cannot change the current predictor", {
  u <- packed(c(10,10,9,8),c(10,9,8,7),c(1,2,2,3),c(0,1,1,2),1:4)
  v <- u;v$x[3:4] <- c(0L,3L)
  a <- sim2_comparators(u,10,10,4,keep_path=TRUE)$prequential
  b <- sim2_comparators(v,10,10,4,keep_path=TRUE)$prequential
  expect_equal(a$eta_path[1:3],b$eta_path[1:3],tolerance=1e-12)
})

test_that("simultaneous cells share a predictor and decisions occur after the release", {
  u <- packed(c(10,8,9,6),c(10,8,9,6),c(1,2,3,2),c(0,0,3,0),c(1,2,2,3))
  v <- lapply(u,function(z) if(is.matrix(z)) z[,c(1,3,2,4),drop=FALSE] else z)
  a <- sim2_comparators(u,10,10,3,keep_path=TRUE)
  b <- sim2_comparators(v,10,10,3,keep_path=TRUE)
  for (name in names(a)) {
    expect_equal(a[[name]]$terminal_loge,b[[name]]$terminal_loge,tolerance=1e-11)
    expect_equal(a[[name]]$maximum_loge,b[[name]]$maximum_loge,tolerance=1e-11)
    expect_equal(a[[name]]$stop_time,b[[name]]$stop_time)
    expect_equal(a[[name]]$eta_path[2],a[[name]]$eta_path[3],tolerance=1e-12)
    expect_equal(a[[name]]$loge_path[2],a[[name]]$loge_path[3],tolerance=1e-12)
  }
})

test_that("all degenerate margins and padded cells contribute one", {
  u <- packed(c(0,4,3,0),c(4,0,2,0),c(2,1,5,0),c(0,1,3,0),c(1,2,3,Inf),count=3)
  out <- sim2_comparators(u,4,4,3,keep_path=TRUE)
  for (p in out) {
    expect_identical(p$reject,FALSE)
    expect_equal(p$terminal_loge,0)
    expect_equal(p$maximum_loge,0)
    expect_equal(p$stop_time,3)
    expect_equal(p$loge_path,matrix(0,1,4))
  }
  empty <- packed(0,0,0,0,Inf,count=0)
  expect_equal(sim2_comparators(empty,0,0,3)$prequential$terminal_loge,0)
})

test_that("adaptive conditional factors have null expectation one on an enumerated tree", {
  total <- c(fixed=0,prequential=0,prequential_protective=0)
  # Four exact observations from two groups of three; enumerate every possible
  # allocation, evaluating the unconditional central-hypergeometric probability.
  for (mask in 0:15) {
    xs <- as.integer(intToBits(mask)[1:4]);a<-3L;b<-3L;prob<-1
    aa<-bb<-integer(4)
    for (j in 1:4) {
      aa[j]<-a;bb[j]<-b
      if (xs[j]>b || 1L-xs[j]>a) { prob<-0;break }
      prob<-prob*if(xs[j]==1L)b/(a+b) else a/(a+b)
      b<-b-xs[j];a<-a-(1L-xs[j])
    }
    if(prob==0)next
    out<-sim2_comparators(packed(bb,aa,rep(1L,4),xs,1:4),3,3,4)
    for(name in names(total)) total[name]<-total[name]+prob*exp(out[[name]]$terminal_loge)
  }
  expect_equal(unname(total),rep(1,3),tolerance=1e-11)
})

test_that("invalid margins and unordered release times are rejected", {
  u<-packed(c(3,2),c(3,3),c(1,2),c(0,1),c(2,1))
  expect_error(sim2_comparators(u,3,3,3),"Release times")
  u$time[,]<-c(1,2);u$x[2]<-3L
  expect_error(sim2_comparators(u,3,3,3),"Invalid fixed-margin")
})
cat("SIM-2 deterministic comparator tests completed. No simulations run.\n")
