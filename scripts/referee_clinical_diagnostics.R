#!/usr/bin/env Rscript
# Recompute the archived clinical paths, informative-table counts, coefficient
# ranges and the eight reference PsA contributions. Original files are read-only.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
ROOT <- normalizePath(if (length(args)) args[[1L]] else ".")
OUT <- file.path(ROOT, "results/session/referee_2026-10-22/clinical")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
datadir <- file.path(ROOT, "reference/liming_revision/data")
cells <- fread(file.path(datadir, "application_fixed_margin_updates.csv"))
archive <- fread(file.path(datadir, "application_weight_eprocess_paths.csv"))
summary <- fread(file.path(datadir, "application_weight_summary.csv"))
codes <- c("LR", "GB", "PP", "FH", "GT")
lse <- function(x) {a <- max(x); a + log(sum(exp(x-a)))}
cells[, n := n0 + n1]
cells[, variance := fifelse(n > 1, d * n1 * n0 * (n-d) / (n^2 * pmax(n-1,1)), 0)]
counts <- cells[, .(dataset=dataset[1L], q=q[1L], unit=q_unit[1L],
  intervals=sum(interval_records), comparisons=.N, looks=uniqueN(look_time),
  events=sum(d), informative_comparisons=sum(variance > 0),
  informative_events=sum(d[variance > 0]), zero_variance_events=sum(d[variance == 0])),
  by=analysis]
fwrite(counts, file.path(OUT, "comparison_counts.csv"))
ranges <- list(); contributions <- list(); recomputed <- list(); i <- 0L; cidx <- 0L; pidx <- 0L
for (a in unique(cells$analysis)) {
  cc <- copy(cells[analysis == a]); setorder(cc, look_time, update_id)
  eta <- if(cc$dataset[1L] == "cav") c(-log(2),-log(1.5),log(1.5),log(2)) else log(c(.5,.67,.8))
  for (w in codes) {
    weight <- cc[[paste0("weight_",w)]]
    i <- i+1L
    informative <- cc$variance > 0
    ranges[[i]] <- data.table(analysis=a,weight=w,weight_min=min(weight),weight_max=max(weight),
      eta_weight_min=min(outer(weight,eta)),eta_weight_max=max(outer(weight,eta)),
      informative_weight_min=min(weight[informative]),informative_weight_max=max(weight[informative]),
      informative_eta_weight_min=min(outer(weight[informative],eta)),
      informative_eta_weight_max=max(outer(weight[informative],eta)))
    components <- rep(0,length(eta)); old_log <- 0
    for (j in seq_len(nrow(cc))) {
      r <- cc[j]; yy <- seq.int(max(0,r$d-r$n0), min(r$d,r$n1))
      lp <- dhyper(yy,r$n1,r$n0,r$d,log=TRUE)
      inc <- vapply(eta*weight[j],function(u) u*r$x-lse(lp+u*yy),numeric(1))
      components <- components + inc
      current_log <- lse(components)-log(length(eta))
      if (a == "psor_q6_primary") {
        cidx <- cidx+1L
        contributions[[cidx]] <- data.table(analysis=a,weight=w,update_id=r$update_id,
          left=r$rounded_left,right=r$rounded_right,effusion=r$effusion_band,
          n0=r$n0,n1=r$n1,events0=r$d-r$x,events1=r$x,total_events=r$d,
          variance=r$variance,weight_value=weight[j],
          component_log_increment=paste(format(inc,digits=17),collapse=";"),
          mixture_multiplier=exp(current_log-old_log),
          cumulative_evidence=exp(current_log),log_mixture_increment=current_log-old_log)
      }
      if (j == nrow(cc) || cc$look_time[j+1L] != r$look_time) {
        pidx <- pidx+1L
        recomputed[[pidx]] <- data.table(analysis=a,weight_code=w,look_time=r$look_time,log_e_value=current_log)
      }
      old_log <- current_log
    }
  }
}
ranges <- rbindlist(ranges); contributions <- rbindlist(contributions)
recomputed <- rbindlist(recomputed)
mix <- merge(recomputed[weight_code=="LR"],recomputed[weight_code=="GT"],
  by=c("analysis","look_time"),suffixes=c("_lr","_gt"))
mix[, log_e_value := mapply(function(x,y) lse(c(x,y))-log(2),log_e_value_lr,log_e_value_gt)]
mix[,weight_code:="LR--GT"]
recomputed <- rbind(recomputed,mix[,.(analysis,weight_code,look_time,log_e_value)])
checks <- merge(recomputed,archive[,.(analysis,weight_code,look_time,archived_log_e_value=log_e_value)],
  by=c("analysis","weight_code","look_time"),all=TRUE)
stopifnot(nrow(checks)==nrow(archive),!anyNA(checks),max(abs(checks$log_e_value-checks$archived_log_e_value))<1e-8)
fwrite(ranges,file.path(OUT,"coefficient_weight_ranges.csv"))
fwrite(contributions,file.path(OUT,"psa_reference_table_contributions.csv"))
fwrite(checks,file.path(OUT,"independent_path_checks.csv"))
fwrite(cells[,.(analysis,update_id,variance)],file.path(OUT,"conditional_variances.csv"))
write_json(list(archived_path_values=nrow(archive),recomputed_path_values=nrow(recomputed),
  maximum_absolute_log_error=max(abs(checks$log_e_value-checks$archived_log_e_value)),
  informative_definition="positive conditional hypergeometric variance: n0*n1*d*(n0+n1-d)>0",
  contribution_order="saved update_id within look_time; intermediate mixtures are algebraic decompositions, not additional monitoring looks",
  original_data_files_unchanged=TRUE,patient_level_preprocessing_rerun=FALSE),
  file.path(OUT,"diagnostics_verification.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)

# Use the same mixed linear/log display as the original figures, with an
# explicit reference-threshold label. Coordinates are transformed, not data.
transform_y <- function(x) ifelse(x <= 1,x/3,1/3+log10(pmax(x,1)))
plot_one <- function(a,filename,title,xlab,ticks,upper) {
  d <- archive[analysis == a]
  draw <- function() {
    par(mar=c(4.2,4.4,4.7,1),mgp=c(2.5,.7,0),las=1,family="sans",cex=.95)
    plot(NA,xlim=c(0,max(d$look_time)),ylim=c(0,transform_y(upper)),
      xlab=xlab,ylab="Accumulated e-value",axes=FALSE)
    axis(1,at=ticks); yt<-c(0,1,5,10,20,if(upper>100)c(100,500))
    axis(2,at=transform_y(yt),labels=yt); box(bty="l")
    abline(h=transform_y(1),lty=3,col="grey60")
    abline(h=transform_y(20),lty=2,col="grey25")
    text(.02*max(d$look_time),transform_y(20),"Reference threshold (20)",adj=c(0,-.4),cex=.85)
    cols<-c("#185778","#BB6631","#368578","#9263A2","#AE3D47","#777777")
    allcodes<-c(codes,"LR--GT")
    for(k in seq_along(allcodes)) {
      dd<-d[weight_code==allcodes[k]][order(look_index)]
      lines(c(0,dd$look_time),transform_y(c(1,dd$e_value)),type="s",col=cols[k],lty=k,lwd=2.8)
      if(startsWith(a,"psor")) points(dd$look_time,transform_y(dd$e_value),pch=16,cex=.4,col=cols[k])
    }
    legend("topleft",inset=c(0,-.20),legend=sub("--","-",allcodes),col=cols,lty=seq_along(allcodes),
      lwd=2.8,ncol=3,bty="n",xpd=NA,cex=.9)
    mtext(title,side=3,line=3.6,adj=0,cex=1.05)
  }
  grDevices::pdf(file.path(ROOT,"figures",paste0(filename,".pdf")),width=7.6,height=4.8,useDingbats=FALSE)
  draw();dev.off()
  grDevices::png(file.path(OUT,paste0(filename,".png")),width=1368,height=864,res=180)
  draw();dev.off()
}
plot_one("cav_q0.25_primary","referee_cav_evidence","CAV: assessment width 0.25 years", "Aligned follow-up after transplantation (years)",seq(0,18,3),500)
plot_one("psor_q6_primary","referee_psa_evidence","PsA: assessment width 6 months", "Nominal follow-up after baseline (months)",seq(0,24,6),40)
print(counts)
print(ranges[analysis %in% c("cav_q0.25_primary","psor_q6_primary")])
print(contributions[weight %in% c("LR","GT"),.(weight,update_id,left,right,effusion,n0,n1,events0,events1,variance,mixture_multiplier,cumulative_evidence)])
