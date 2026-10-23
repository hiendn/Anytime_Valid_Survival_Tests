#!/usr/bin/env Rscript
# Bounded null stress test of the actual CAV rounding and PsA last-visit binning
# functions. Synthetic data only; no archived results or patient files changed.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
argv <- commandArgs(trailingOnly=TRUE)
getarg <- function(flag,default) {i<-match(flag,argv);if(is.na(i)) default else argv[i+1L]}
ROOT <- normalizePath(getarg("--root","."))
REPS <- as.integer(getarg("--reps","2000"))
CORES <- as.integer(getarg("--cores","4"))
OUT <- file.path(ROOT,getarg("--output","results/session/referee_2026-10-22/clinical/stress"))
stopifnot(REPS>=1L,REPS<=5000L)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
data.table::setDTthreads(1)
source(file.path(ROOT,"simulations/code/simulation_core.R"))
application_source <- file.path(ROOT,"simulations/code/run_real_data_applications.R")
# Evaluate only the original function definitions. This avoids the historical
# driver's top-level path setup and main(), and leaves all bodies unmodified.
functions <- c("log_mean_exp","eta_text","mixture_path","add_application_weights",
  "cav_intervals","cav_analysis","psor_cohort","psor_analysis")
loaded <- character()
for (e in parse(application_source)) {
  if(is.call(e) && identical(e[[1L]],as.name("<-")) &&
     is.symbol(e[[2L]]) && as.character(e[[2L]]) %in% functions &&
     is.call(e[[3L]]) && identical(e[[3L]][[1L]],as.name("function"))) {
    eval(e);loaded<-c(loaded,as.character(e[[2L]]))
  }
}
stopifnot(setequal(loaded,functions))
ensure_simulation_kernels()
mechanisms <- c("common_grid","independent_jitter","group_differential_jitter")
seed_base <- 202610220L

synthetic <- function(application,replicate_id,mechanism) {
  cav <- application=="CAV"; n<-if(cav)583L else 185L
  n0<-if(cav)270L else 120L; steps<-if(cav)6L else 4L
  lambda<-if(cav).12 else .025; q<-if(cav).25 else 6
  set.seed(seed_base+replicate_id+if(cav)0L else 100000L)
  group<-c(rep(0L,n0),rep(1L,n-n0)); event_time<-rexp(n,lambda)
  if(cav) {
    visits<-CJ(id=seq_len(n),k=seq_len(steps))
    visits[,u:=runif(.N)]
    visits[,group:=group[id]]
    visits[,right:=k+switch(mechanism,common_grid=0,
      independent_jitter=.24*u-.12,
      group_differential_jitter=.12*(u-1+group))]
    visits[,left:=shift(right,fill=0),by=id]
    visits[,positive:=as.integer(event_time[id]<=right)]
    visits[,keep:=seq_len(.N)<=if(any(positive==1L))which(positive==1L)[1L] else .N,by=id]
    intervals<-visits[keep==TRUE,.(id,group,left,right,event=positive,
      sex=2L,recipient_age=50,donor_age=35)]
    subjects<-intervals[,.(group=group[1L],outcome=if(any(event==1L))"CAV" else "censored",
      completed_intervals=.N),by=id]
    baseline<-data.table(PTNUM=seq_len(n),sex=2L,age=50,dage=35)
    res<-cav_analysis(intervals,subjects,baseline,q,primary=TRUE)
    probability<-copy(intervals)
    probability[,`:=`(rounded_left=q*round(left/q),rounded_right=q*round(right/q),
      p=-expm1(-lambda*(right-left)))]
    probability<-probability[rounded_right>rounded_left]
    variation<-probability[,.(range=max(p)-min(p)),by=.(rounded_left,rounded_right)]$range
    if(replicate_id==1L && mechanism=="independent_jitter") {
      # Verify the vectorized construction of completed intervals against the
      # actual subject-level CAV preprocessing on the same synthetic records.
      raw<-rbind(data.table(PTNUM=seq_len(n),years=0,state=1L,pdiag=ifelse(group==1L,"IHD","IDC"),
        sex=2L,age=50,dage=35),
        visits[,.(PTNUM=id,years=right,state=ifelse(positive==1L,2L,1L),
          pdiag=ifelse(group==1L,"IHD","IDC"),sex=2L,age=50,dage=35)])
      reference<-cav_intervals(raw)$intervals
      stopifnot(isTRUE(all.equal(reference[,.(id,group,left,right,event)],
        intervals[,.(id,group,left,right,event)],check.attributes=FALSE)))
    }
  } else {
    visits<-CJ(ptnum=seq_len(n),k=seq_len(steps),visit=1:2)
    visits[,`:=`(u=runif(.N),group=group[ptnum])]
    visits[,offset:=switch(mechanism,common_grid=5+as.numeric(visit==2L),
      independent_jitter=.2+5.6*u,
      group_differential_jitter=.2+2.6*u+3*(1-group))]
    visits[,aligned_months:=6*(k-1)+offset]
    visits[,`:=`(baseline_effusion=0L,state=ifelse(event_time[ptnum]<=aligned_months,2L,1L))]
    setorder(visits,ptnum,aligned_months)
    records<-visits[,.(ptnum,aligned_months,group,baseline_effusion,state)]
    raw_intervals<-records[,{
      ix<-which(state>1L);last<-if(length(ix))ix[1L] else .N
      z<-copy(.SD[seq_len(last)]);z[,event:=as.integer(state>1L)];z
    },by=ptnum]
    subjects<-raw_intervals[,.(id=ptnum[1L],group=group[1L],baseline_effusion=0L,
      outcome=if(any(event==1L))"progression" else "censored",completed_intervals=.N),by=ptnum]
    baseline<-data.table(group=group,baseline_effusion=0L)
    res<-psor_analysis(records,raw_intervals,subjects,baseline,q,primary=TRUE)
    pp<-copy(res$data)
    pp[,p:=-expm1(-lambda*(aligned_months-previous_actual))]
    variation<-pp[,.(range=max(p)-min(p)),by=.(rounded_left,rounded_right)]$range
    if(replicate_id==1L && mechanism=="independent_jitter") {
      raw<-rbind(data.table(ptnum=seq_len(n),months=0,state=1L,
        ollwsdrt=ifelse(group==1L,1L,0L),hieffusn=0L),
        records[,.(ptnum,months=aligned_months,state,ollwsdrt=ifelse(group==1L,1L,0L),hieffusn=0L)])
      reference<-psor_cohort(raw)
      stopifnot(sum(reference$raw_intervals$event)==sum(raw_intervals$event),
        nrow(reference$raw_intervals)==nrow(raw_intervals))
    }
  }
  lr<-res$path
  gt<-mixture_path(res$cells,res$eta,unique(res$cells$analysis),
    unique(res$cells$display_name),res$cells$weight_GT)
  cc<-res$cells
  informative<-cc$n0>0 & cc$n1>0 & cc$d>0 & cc$d<(cc$n0+cc$n1)
  rbindlist(lapply(c("LR","GT"),function(w){p<-if(w=="LR")lr else gt
    data.table(application=application,replicate=replicate_id,mechanism=mechanism,weight=w,
      crossed20=any(p$e_value>=20),maximum_log_e=max(0,p$log_e_value),
      terminal_log_e=tail(p$log_e_value,1),events=sum(cc$d),comparisons=nrow(cc),
      informative_comparisons=sum(informative),informative_events=sum(cc$d[informative]),
      max_probability_range=max(variation),mean_probability_range=mean(variation))}))
}

results<-list();idx<-0L;started<-proc.time()[[3L]]
for(app in c("CAV","PsA")) {
  for(first in seq.int(1L,REPS,100L)) {
    rr<-seq.int(first,min(REPS,first+99L))
    batch<-parallel::mclapply(rr,function(r) rbindlist(lapply(mechanisms,
      function(m) synthetic(app,r,m))),mc.cores=min(CORES,length(rr)),mc.set.seed=FALSE)
    stopifnot(all(vapply(batch,is.data.frame,logical(1))))
    idx<-idx+1L;results[[idx]]<-rbindlist(batch)
    cat(app,max(rr),"/",REPS,"elapsed seconds",round(proc.time()[[3L]]-started),"\n");flush.console()
    fwrite(rbindlist(results),file.path(OUT,"replicate_results.csv"))
  }
}
result<-rbindlist(results)
wilson<-function(k,n) {
  z<-qnorm(.975);p<-k/n;den<-1+z*z/n;centre<-(p+z*z/(2*n))/den
  half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;c(centre-half,centre+half)
}
summary<-result[,{
  k<-sum(crossed20);ci<-wilson(k,.N)
  list(replicates=.N,crossings=k,rate=k/.N,mcse=sqrt((k/.N)*(1-k/.N)/.N),
    wilson_lower=ci[1L],wilson_upper=ci[2L],mean_events=mean(events),
    mean_informative_comparisons=mean(informative_comparisons),
    median_max_probability_range=median(max_probability_range),
    mean_within_comparison_probability_range=mean(mean_probability_range))
},by=.(application,mechanism,weight)]
fwrite(summary,file.path(OUT,"summary.csv"))
fwrite(CJ(application=c("CAV","PsA"),replicate=seq_len(REPS))[
  ,seed:=seed_base+replicate+ifelse(application=="CAV",0L,100000L)],file.path(OUT,"seed_ledger.csv"))
write_json(list(replicates_per_application_mechanism=REPS,
  actual_application_function_names=functions,
  actual_application_source_md5=unname(tools::md5sum(application_source)),
  data_model="Independent exponential event times with a common rate in both groups; visits independent of event times; no loss to follow-up; one homogeneous prognostic stratum",
  cav=list(n=583,n0=270,n1=313,rate_per_year=.12,assessment_years=1:6,q=.25,
    jitter="Independent: Uniform(-.12,.12); differential: group 0 Uniform(-.12,0), group 1 Uniform(0,.12); added to each annual visit"),
  psa=list(n=185,n0=120,n1=65,rate_per_month=.025,windows=4,q=6,
    visits_per_window=2,
    jitter="Within each six-month window, common grid offsets 5 and 6; independent offsets Uniform(.2,5.8); differential offsets group 0 Uniform(3.2,5.8), group 1 Uniform(.2,2.8); last visit retained by actual function"),
  probability_diagnostic="1-exp(-lambda*(actual_right-actual_left)); maximum minus minimum among intervals assigned to each comparison; summarize over replications",
  paired_design="Event-time draws and underlying uniform visit draws are reused across mechanisms within each application/replicate",
  statistical_scope="A mechanism stress test under equal latent event hazards, not a validation of either clinical cohort or its observation assumptions",
  no_archived_results_modified=TRUE,elapsed_seconds=proc.time()[[3L]]-started,
  session_info=capture.output(sessionInfo())),file.path(OUT,"design_and_validation.json"),
  auto_unbox=TRUE,pretty=TRUE,digits=NA)
print(summary)
