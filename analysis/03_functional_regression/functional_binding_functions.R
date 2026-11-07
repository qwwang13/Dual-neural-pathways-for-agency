suppressPackageStartupMessages({library(dplyr);library(patchwork)})

binding_settings <- function(binding, rules) {
  binding <- match.arg(binding,c("action","outcome"))
  cfg <- c(list(binding=binding,time_windows=c(10,20,40),lower=0,upper=1000,
    condition=if(binding=="action")"vaa_iaa"else"vss_iss",
    binding_col=if(binding=="action")"vaa_iaa_binding"else"vss_iss_binding"),rules)
  stopifnot(setequal(names(rules),c("conti_p","gap_ms","min_span_ms","density_cap","outside_runs_cap")),
    setequal(names(cfg$conti_p),c("10","20","40")),setequal(names(cfg$gap_ms),c("10","20","40")),
    all(is.finite(cfg$gap_ms)),all(cfg$gap_ms>=0),all(cfg$gap_ms%%10==0),
    length(cfg$min_span_ms)==1,cfg$min_span_ms>=10,cfg$min_span_ms%%10==0,
    length(cfg$density_cap)==1,cfg$density_cap>0,cfg$density_cap<=1,
    length(cfg$outside_runs_cap)==1,cfg$outside_runs_cap>=0,cfg$outside_runs_cap%%1==0)
  for(tw in cfg$time_windows) {
    p <- cfg$conti_p[[as.character(tw)]]
    stopifnot(length(p)>0,all(is.finite(p)),all(p>=1),all(p%%1==0),!anyDuplicated(p))
    cfg$conti_p[[as.character(tw)]] <- sort(p)
  }
  cfg
}

binding_strict_run <- function(pos) {
  if(!length(pos))return(0L)
  pos<-sort(unique(pos));max(tabulate(cumsum(c(TRUE,diff(pos)!=1))))
}

binding_intervals <- function(mask,grid,gap_ms=0) {
  idx<-which(mask)
  if(!length(idx))return(tibble::tibble(lower=numeric(),upper=numeric()))
  step<-diff(grid)[1]
  groups<-cumsum(c(TRUE,diff(idx)*step>step+gap_ms))
  dplyr::bind_rows(lapply(split(idx,groups),function(z)tibble::tibble(lower=min(grid[z]),upper=max(grid[z])+step)))
}

binding_scan_models <- function(models,cfg) {
  stopifnot(all(models$intercept))
  models<-models%>%dplyr::mutate(path_index=dplyr::row_number(),strict_run=purrr::map_int(pos,binding_strict_run))
  status<-list();coefficients<-list()
  for(rg in sort(unique(models$region)))for(tw in cfg$time_windows)for(p in cfg$conti_p[[as.character(tw)]]) {
    pool<-models%>%dplyr::filter(region==rg,time_window==tw,strict_run>=p)
    chosen<-pool%>%dplyr::arrange(n_selected,dplyr::desc(r2),alpha,path_index)%>%dplyr::slice_head(n=1)
    found<-nrow(chosen)==1;pos<-if(found)chosen$pos[[1]]else integer()
    beta<-if(found)chosen$coef[[1]]else numeric();total<-(cfg$upper-cfg$lower)/tw
    status[[length(status)+1]]<-tibble::tibble(region=rg,time_window=tw,conti_p=p,
      candidate_found=found,n_eligible=nrow(pool),n_selected=length(pos),total_bins=total,
      covered_ms=length(pos)*tw,selected_fraction=length(pos)/total,
      alpha=if(found)chosen$alpha else NA_real_,lambda=if(found)chosen$lambda else NA_real_,
      r2=if(found)chosen$r2 else NA_real_,intercept=if(found)TRUE else NA)
    if(found) {
      stopifnot(length(pos)==chosen$n_selected,length(beta)==length(pos),binding_strict_run(pos)>=p)
      coefficients[[length(coefficients)+1]]<-tibble::tibble(region=rg,time_window=tw,conti_p=p,
        alpha=chosen$alpha,lambda=chosen$lambda,r2=chosen$r2,intercept=TRUE,consecutive_mode="strict",
        p_number=pos,lower=cfg$lower+(pos-1)*tw,upper=cfg$lower+pos*tw,coef=beta,
        coef_display=beta/(tw/10))
    }
  }
  list(summary=dplyr::bind_rows(status),coefficients=dplyr::bind_rows(coefficients))
}

binding_select_windows <- function(scan,cfg) {
  grid<-seq(cfg$lower,cfg$upper-10,10);all_combinations<-list();choices<-list();chosen_coef<-list()
  for(rg in sort(unique(scan$summary$region))) {
    sets<-list()
    for(tw in cfg$time_windows)for(p in cfg$conti_p[[as.character(tw)]]) {
      s<-scan$summary%>%dplyr::filter(region==rg,time_window==tw,conti_p==p)
      cf<-scan$coefficients%>%dplyr::filter(region==rg,time_window==tw,conti_p==p,coef!=0)
      raw<-vapply(grid,function(t)any(cf$lower<=t & cf$upper>t),logical(1))
      env<-binding_intervals(raw,grid,cfg$gap_ms[[as.character(tw)]])
      envelope<-vapply(grid,function(t)any(env$lower<=t & env$upper>t),logical(1))
      nextrow<-scan$summary%>%dplyr::filter(region==rg,time_window==tw,conti_p==p+1)
      jump<-if(nrow(nextrow)==1 && isTRUE(nextrow$candidate_found) && isTRUE(s$candidate_found))100*(nextrow$selected_fraction-s$selected_fraction)else NA_real_
      sets[[paste(tw,p,sep="|")]]<-list(raw=raw,envelope=envelope,summary=s,jump=jump)
    }
    combinations<-tidyr::crossing(p10=cfg$conti_p$`10`,p20=cfg$conti_p$`20`,p40=cfg$conti_p$`40`)
    rows<-vector("list",nrow(combinations))
    for(i in seq_len(nrow(combinations))) {
      ps<-as.integer(unlist(combinations[i,],use.names=FALSE))
      entries<-Map(function(tw,p)sets[[paste(tw,p,sep="|")]],cfg$time_windows,ps)
      found<-all(vapply(entries,function(e)isTRUE(e$summary$candidate_found),logical(1)))
      actual<-Reduce(`&`,lapply(entries,`[[`,"raw"))
      common<-Reduce(`&`,lapply(entries,`[[`,"envelope"))
      intervals<-binding_intervals(common,grid)
      admissible<-intervals%>%dplyr::filter(upper-lower>=cfg$min_span_ms)
      if(!found)admissible<-admissible[0,]
      candidates<-if(nrow(admissible))admissible else intervals
      primary_lower<-primary_upper<-NA_real_;span<-support<-0
      outside_runs<-nrow(binding_intervals(entries[[1]]$raw,grid));outside_bins<-sum(entries[[1]]$raw)
      if(nrow(candidates)) {
        local<-candidates%>%dplyr::mutate(span_ms=upper-lower,
          actual_ms=vapply(seq_len(nrow(candidates)),function(j)sum(actual & grid>=candidates$lower[j] & grid<candidates$upper[j])*10,numeric(1)),
          off_runs=vapply(seq_len(nrow(candidates)),function(j){
            z<-entries[[1]]$raw & (grid<candidates$lower[j] | grid>=candidates$upper[j]);sum(z & !c(FALSE,head(z,-1)))
          },numeric(1)),
          off_bins=vapply(seq_len(nrow(candidates)),function(j)sum(entries[[1]]$raw & (grid<candidates$lower[j] | grid>=candidates$upper[j])),numeric(1)))
        localized<-local%>%dplyr::filter(off_runs<=cfg$outside_runs_cap)
        if(!nrow(localized))localized<-local
        primary<-localized%>%dplyr::arrange(dplyr::desc(span_ms),dplyr::desc(actual_ms),off_runs,off_bins,lower)%>%dplyr::slice_head(n=1)
        primary_lower<-primary$lower;primary_upper<-primary$upper;span<-primary$span_ms;support<-primary$actual_ms
        outside_runs<-primary$off_runs;outside_bins<-primary$off_bins
      }
      counts<-vapply(entries,function(e)e$summary$n_selected,numeric(1))
      fraction<-vapply(entries,function(e)e$summary$selected_fraction,numeric(1))
      jumps<-vapply(entries,`[[`,numeric(1),"jump")
      rows[[i]]<-tibble::tibble(region=rg,combination_id=i,p10=ps[1],p20=ps[2],p40=ps[3],
        all_candidates_found=found,qualifies=nrow(admissible)>0,n10=counts[1],n20=counts[2],n40=counts[3],
        fraction10=fraction[1],fraction20=fraction[2],fraction40=fraction[3],max_nonzero_fraction=max(fraction),
        mean_nonzero_fraction=if(found)mean(fraction)else NA_real_,
        min_jump_pp=if(all(is.finite(jumps)))min(jumps)else NA_real_,
        primary_lower_ms=primary_lower,primary_upper_ms=primary_upper,primary_span_ms=span,
        primary_actual_overlap_ms=support,off_target_runs10=outside_runs,off_target_bins10=outside_bins,
        primary_window=if(is.finite(primary_lower))paste0(primary_lower,"–",primary_upper)else"none",
        max_overlap_span_ms=if(nrow(intervals))max(intervals$upper-intervals$lower)else 0,
        localization_pass=outside_runs<=cfg$outside_runs_cap)
    }
    tab<-dplyr::bind_rows(rows)%>%dplyr::mutate(auto_passed=all_candidates_found & qualifies & localization_pass & max_nonzero_fraction<=cfg$density_cap+1e-12)
    passed<-tab%>%dplyr::filter(auto_passed)
    if(nrow(passed)) {
      choice<-passed%>%dplyr::arrange(dplyr::desc(primary_span_ms),dplyr::desc(primary_actual_overlap_ms),
        off_target_runs10,off_target_bins10,max_nonzero_fraction,dplyr::desc(min_jump_pp),mean_nonzero_fraction,p10,p20,p40)%>%dplyr::slice_head(n=1)
    } else {
      pool<-tab%>%dplyr::filter(all_candidates_found,max_nonzero_fraction<=cfg$density_cap+1e-12)
      short<-pool%>%dplyr::filter(max_overlap_span_ms>0,max_overlap_span_ms<30)
      if(nrow(short))pool<-short
      if(!nrow(pool))pool<-tab%>%dplyr::filter(all_candidates_found)
      if(!nrow(pool))stop("No three-resolution model combination for ",rg)
      choice<-pool%>%dplyr::arrange(max_nonzero_fraction,mean_nonzero_fraction,dplyr::desc(max_overlap_span_ms),dplyr::desc(min_jump_pp),p10,p20,p40)%>%dplyr::slice_head(n=1)
    }
    cf<-dplyr::bind_rows(lapply(cfg$time_windows,function(tw){
      p<-choice[[paste0("p",tw)]];scan$coefficients%>%dplyr::filter(region==rg,time_window==tw,conti_p==p)
    }))
    for(tw in cfg$time_windows){rr<-unique(cf$r2[cf$time_window==tw]);choice[[paste0("r2_",tw)]]<-if(length(rr))rr else NA_real_}
    choices[[rg]]<-choice;all_combinations[[rg]]<-tab;chosen_coef[[rg]]<-cf
  }
  list(choices=dplyr::bind_rows(choices),combinations=dplyr::bind_rows(all_combinations),coefficients=dplyr::bind_rows(chosen_coef))
}
