# Pure helpers: finite-tolerance endpoint-conditioned CTMC uniformization.
# Sourcing this file does not create an output, change RNG, or run a map.
bgb_bridge_logsum <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x-m)))
}

bgb_bridge_prepare_q <- function(Qmat, ctx) {
  if(!is.null(ctx$last_Q) && identical(Qmat,ctx$last_Q)) return(ctx$last_prepared_Q)
  qhash <- digest::digest(Qmat, algo="sha256", serialize=TRUE)
  if (exists(qhash, ctx$qcache, inherits=FALSE)) {
    z<-get(qhash,ctx$qcache);ctx$last_Q<-Qmat;ctx$last_prepared_Q<-z;return(z)
  }
  Q <- as.matrix(Qmat)
  if (!is.numeric(Q) || length(dim(Q))!=2L || nrow(Q)!=ncol(Q) || !nrow(Q) || any(!is.finite(Q))) stop("Bridge requires finite square Q")
  off <- Q; diag(off) <- 0
  if (any(off<0) || any(diag(Q)>0)) stop("Bridge Q has invalid rate sign; Q is not repaired")
  rs <- rowSums(Q)
  if (max(rs)>1e-12) stop("Bridge Q has positive row mass exceeding roundoff tolerance")
  mu <- max(-diag(Q)); n <- nrow(Q)
  if (mu==0 && any(off!=0)) stop("Zero diagonal with nonzero outgoing rates")
  R <- if(mu==0) Matrix::Diagonal(n) else Matrix::Diagonal(n)+Matrix::Matrix(Q/mu,sparse=TRUE)
  if(any(R@x<0)) stop("Bridge R has negative entries; R is not repaired")
  # rho also covers positive row-sum floating point error conservatively.
  rho <- max(1, as.numeric(Matrix::rowSums(R)))
  z <- list(Q_hash=qhash, R=R, mu=mu, n=n, rho=rho,
            min_Q_row_sum=min(rs), max_Q_row_sum=max(rs), ladders=new.env(parent=emptyenv()))
  assign(qhash,z,ctx$qcache); ctx$q_count <- ctx$q_count+1L
  ctx$last_Q<-Qmat;ctx$last_prepared_Q<-z
  z
}

bgb_native_effective_generator <- function(Qraw,abs_tol=1e-7,rel_tol=4*2^-23) {
  # Reproduce native waiting rates AND row-normalized destinations explicitly.
  # Reject material missing/killing mass: this is not generic Q repair.
  Q<-as.matrix(Qraw)
  if(!is.numeric(Q) || nrow(Q)!=ncol(Q) || !nrow(Q) || any(!is.finite(Q))) stop("Native effective Q requires finite square input")
  off<-Q; diag(off)<-0
  rate<--diag(Q); total<-rowSums(off)
  if(any(off<0) || any(rate<0) || any(rate==0 & total>0) || any(rate>0 & total==0)) stop("Invalid native waiting-rate/destination combination")
  residual<-total-rate
  relative_residual<-ifelse(pmax(total,rate)>0,abs(residual)/pmax(total,rate),0)
  if(max(abs(residual))>abs_tol || max(relative_residual)>rel_tol) stop("Qraw discrepancy exceeds allowed single-rounding tolerance; do not normalize killing away")
  factor<-rep(1,nrow(Q)); hit<-total>0; factor[hit]<-rate[hit]/total[hit]
  effective<-off*factor; diag(effective)<--rate
  relative_change<-if(any(hit))max(abs(factor[hit]-1)) else 0
  list(Q=effective,metadata=list(Qraw_sha256=digest::digest(Qraw,algo="sha256",serialize=TRUE),
    Qeffective_sha256=digest::digest(effective,algo="sha256",serialize=TRUE),
    native_Q_max_abs_row_residual=max(abs(residual)),native_Q_max_relative_row_residual=max(relative_residual),
    native_Q_max_abs_change=max(abs(effective-Q)),native_Q_max_relative_offdiag_change=relative_change,
    native_Q_conversion="explicit native waiting-rate plus row-normalized destination law",
    native_Q_abs_tolerance=abs_tol,native_Q_relative_tolerance=rel_tol))
}

bgb_ctmc_bridge <- function(Qmat,t,a,b,ctx) {
  start_clock <- proc.time()[["elapsed"]]
  z <- bgb_bridge_prepare_q(Qmat,ctx)
  if(length(t)!=1L || !is.finite(t) || t<=0 || length(a)!=1L || length(b)!=1L ||
     !is.finite(a) || !is.finite(b) || a!=as.integer(a) || b!=as.integer(b) ||
     a<1 || b<1 || a>z$n || b>z$n) stop("Bridge invalid duration or endpoint")
  a<-as.integer(a); b<-as.integer(b); lambda <- z$mu*t
  if(!is.finite(lambda) || lambda>ctx$max_lambda) stop("Bridge lambda exceeds declared computational bound")
  key <- as.character(b)
  if(!exists(key,z$ladders,inherits=FALSE)) {
    # Reverse reachability uses exactly the positive entries of R; self loops
    # do not create spurious paths. No endpoint probabilities are thresholded.
    reach <- rep(FALSE,z$n); reach[b]<-TRUE
    repeat {
      next_reach <- reach | as.numeric(z$R %*% as.numeric(reach))>0
      if(identical(next_reach,reach)) break
      reach <- next_reach
    }
    ladder <- new.env(parent=emptyenv())
    v0 <- numeric(z$n); v0[b]<-1
    ladder$v <- list(v0); ladder$logscale <- 0; ladder$reachable<-reach
    assign(key,ladder,z$ladders)
    ctx$vector_count<-ctx$vector_count+1L
  }
  ladder <- get(key,z$ladders,inherits=FALSE)
  if(!ladder$reachable[a]) stop("Bridge endpoint is structurally unreachable under original Q")
  logw <- numeric(); logB <- -Inf; logtail <- Inf; m<-0L
  repeat {
    idx <- m+1L
    if(length(ladder$v)<idx) {
      nxt <- as.numeric(z$R %*% ladder$v[[idx-1L]])
      if(any(!is.finite(nxt)) || any(nxt<0)) stop("Bridge vector recurrence invalid")
      scale <- max(nxt)
      if(scale==0) {
        ladder$v[[idx]] <- nxt; ladder$logscale[idx] <- -Inf
      } else {
        ladder$v[[idx]] <- nxt/scale
        ladder$logscale[idx] <- ladder$logscale[idx-1L]+log(scale)
      }
      ctx$vector_count<-ctx$vector_count+1L
      if(ctx$vector_count*z$n*8>ctx$max_vector_bytes) stop("Bridge deterministic vector cache reached declared memory bound")
    }
    val <- ladder$v[[idx]][a]
    logw[idx] <- if(val==0) -Inf else stats::dpois(m,lambda,log=TRUE)+log(val)+ladder$logscale[idx]
    logB <- bgb_bridge_logsum(c(logB,logw[idx]))
    logtail <- lambda*(z$rho-1)+stats::ppois(m,lambda*z$rho,lower.tail=FALSE,log.p=TRUE)
    if(is.finite(logB) && logtail-logB<=log(ctx$rel_tol)) break
    if(m>=ctx$max_terms) stop("Bridge failed endpoint-relative Poisson tail bound within maximum terms")
    m <- m+1L
  }
  norm <- max(logw)
  N <- sample.int(length(logw),size=1L,prob=exp(logw-norm))-1L
  virtual_states <- integer(N+1L); virtual_states[1L]<-a
  if(N>0L) for(k in seq_len(N)) {
    remaining<-N-k; cur<-virtual_states[k]
    row <- as.numeric(z$R[cur,]); back <- ladder$v[[remaining+1L]]
    lp <- log(row)+log(back); maximum<-max(lp)
    if(!is.finite(maximum)) stop("Bridge conditional intermediate state denominator vanished")
    virtual_states[k+1L] <- sample.int(z$n,size=1L,prob=exp(lp-maximum))
  }
  if(virtual_states[N+1L]!=b) stop("Bridge endpoint invariant failed")
  event_times <- if(N) sort(stats::runif(N,min=0,max=t)) else numeric()
  if(length(event_times) && (any(event_times<=0) || any(event_times>=t) || anyDuplicated(event_times))) stop("Bridge event-time floating point degeneracy")
  real <- if(N) which(diff(virtual_states)!=0L) else integer()
  states <- c(a,virtual_states[real+1L]); times<-event_times[real]
  if(length(real)) {
    supported <- vapply(seq_along(real),function(k) as.numeric(z$R[states[k],states[k+1L]])>0,logical(1))
    if(!all(supported)) stop("Bridge generated unsupported original-Q transition")
  }
  diag <- list(Q_sha256=z$Q_hash,t=t,a=a,b=b,n_states=z$n,mu=z$mu,lambda=lambda,
    log_endpoint_probability=logB,endpoint_probability=exp(logB),
    relative_tail_bound=exp(logtail-logB),last_series_term=m,virtual_events=N,
    real_jumps=length(real),self_loops=N-length(real),
    min_Q_row_sum=z$min_Q_row_sum,max_Q_row_sum=z$max_Q_row_sum,
    elapsed_seconds=proc.time()[["elapsed"]]-start_clock)
  list(states=states,times=times,diagnostics=diag)
}

bgb_bridge_branch_metadata <- function(nodenum_at_top_of_branch,trtable,stratified) {
  # Preserves installed native length extraction, including its special
  # subtree-root correction. Independent preflight must reconcile endpoint P.
  i<-nodenum_at_top_of_branch
  if(length(i)!=1L || !is.finite(i) || i!=as.integer(i) || i<1 || i>nrow(trtable)) stop("Invalid branch table row")
  if(!("SUBedge.length" %in% names(trtable))) {
    brlen<-trtable$edge.length[i]
  } else {
    subedge<-trtable$SUBedge.length[i]
    if(is.na(subedge)) subedge<-0
    unresolved<-FALSE
    if(subedge>trtable$reltimept[i]) {
      unresolved<-TRUE
      if(trtable$piececlass[i]=="subbranch") {
        if(is.na(trtable$fossils[i]) || !trtable$fossils[i]) brlen<-trtable$reltimept[i]
        else brlen<-trtable$time_bot[i]-trtable$time_bp[i]
        unresolved<-FALSE
      }
    } else brlen<-subedge
    if(trtable$SUBnode.type[i]=="root" && trtable$piececlass[i]=="subtree") {
      brlen<-trtable$time_bot[i]-trtable$time_bp[i]; unresolved<-FALSE
    }
    if(unresolved) stop("Native subtree branch length is ambiguous; bridge refuses to guess")
  }
  if(length(brlen)!=1L || !is.finite(brlen) || brlen<=0) stop("Bridge branch length is not finite positive")
  older_age <- if(stratified) trtable$time_top[i]+trtable$SUBtime_bp[i]+brlen else trtable$time_bp[i]+brlen
  if(length(older_age)!=1L || !is.finite(older_age)) stop("Bridge older endpoint age is invalid")
  list(t=brlen,older_age=older_age,a=trtable$sampled_states_AT_brbots[i],b=trtable$sampled_states_AT_nodes[i])
}

bgb_bridge_to_native_events <- function(path,nodenum,trtable,meta,state_indices_0based,ranges_list,areas) {
  if(length(path$times)==0L) return(NA)
  cols<-c("nodenum_at_top_of_branch","trynum","brlen","current_rangenum_1based","new_rangenum_1based",
    "current_rangetxt","new_rangetxt","abs_event_time","event_time","event_type","event_txt",
    "new_area_num_1based","lost_area_num_1based","dispersal_to","extirpation_from")
  rows<-lapply(seq_along(path$times),function(k) {
    old<-path$states[k]; new<-path$states[k+1L]
    from<-state_indices_0based[[old]]; to<-state_indices_0based[[new]]
    from<-from[!is.na(from)]; to<-to[!is.na(to)]
    added<-setdiff(to,from); lost<-setdiff(from,to)
    if(length(added)==1L && !length(lost)) {
      type<-"d"; newarea<-added+1L; lostarea<-"-"; dispersal<-areas[newarea]; extirpation<-"-"
    } else if(length(lost)==1L && !length(added)) {
      type<-"e"; newarea<-"-"; lostarea<-lost+1L; dispersal<-"-"; extirpation<-areas[lostarea]
    } else stop("Bridge event is not a one-area DEC gain or loss")
    c(nodenum,1L,meta$t,old,new,ranges_list[[old]],ranges_list[[new]],meta$older_age-path$times[k],
      path$times[k],type,paste0(ranges_list[[old]],"->",ranges_list[[new]]),newarea,lostarea,dispersal,extirpation)
  })
  ans<-as.data.frame(do.call(rbind,rows),stringsAsFactors=FALSE); names(ans)<-cols; rownames(ans)<-NULL
  ans
}

new_bsm_bridge_context <- function(rel_tol=1e-12,max_terms=10000L,max_lambda=1000,max_vector_bytes=512*1024^2,
                                   native_q_abs_tol=1e-7,native_q_rel_tol=4*2^-23) {
  if(length(rel_tol)!=1L || !is.finite(rel_tol) || rel_tol<=0 || rel_tol>=1) stop("Invalid bridge tolerance")
  if(length(max_terms)!=1L || !is.finite(max_terms) || max_terms<1 || max_terms!=as.integer(max_terms) ||
    length(max_lambda)!=1L || !is.finite(max_lambda) || max_lambda<=0 ||
    length(max_vector_bytes)!=1L || !is.finite(max_vector_bytes) || max_vector_bytes<1024) stop("Invalid bridge computation bound")
  ctx<-new.env(parent=emptyenv()); ctx$rel_tol<-rel_tol; ctx$max_terms<-as.integer(max_terms)
  ctx$max_lambda<-max_lambda; ctx$max_vector_bytes<-max_vector_bytes
  ctx$qcache<-new.env(parent=emptyenv()); ctx$native_qcache<-new.env(parent=emptyenv())
  ctx$last_Q<-NULL;ctx$last_prepared_Q<-NULL;ctx$last_Qraw<-NULL;ctx$last_native_Q<-NULL
  ctx$vector_count<-0L; ctx$q_count<-0L; ctx$audit<-list()
  branch <- function(nodenum_at_top_of_branch,trtable,Qmat,state_indices_0based,ranges_list,areas,
                     single_branch=FALSE,stratified=FALSE,maxtries=40000,manual_history_for_difficult_branches=TRUE) {
    meta<-bgb_bridge_branch_metadata(nodenum_at_top_of_branch,trtable,stratified)
    if(length(state_indices_0based)!=nrow(Qmat) || length(ranges_list)!=nrow(Qmat)) stop("Bridge state dictionary does not match Q")
    if(!is.null(ctx$last_Qraw) && identical(Qmat,ctx$last_Qraw)) native_q<-ctx$last_native_Q
    else {
      raw_hash<-digest::digest(Qmat,algo="sha256",serialize=TRUE)
      if(!exists(raw_hash,ctx$native_qcache,inherits=FALSE)) assign(raw_hash,
        bgb_native_effective_generator(Qmat,abs_tol=native_q_abs_tol,rel_tol=native_q_rel_tol),ctx$native_qcache)
      native_q<-get(raw_hash,ctx$native_qcache,inherits=FALSE)
      ctx$last_Qraw<-Qmat;ctx$last_native_Q<-native_q
    }
    path<-bgb_ctmc_bridge(native_q$Q,meta$t,meta$a,meta$b,ctx)
    events<-bgb_bridge_to_native_events(path,nodenum_at_top_of_branch,trtable,meta,state_indices_0based,ranges_list,areas)
    d<-c(path$diagnostics,native_q$metadata)
    d$call_id<-length(ctx$audit)+1L; d$nodenum_at_top_of_branch<-nodenum_at_top_of_branch
    d$master_node<-if("node"%in%names(trtable)) trtable$node[nodenum_at_top_of_branch] else NA_integer_
    d$older_endpoint_age<-meta$older_age; d$younger_endpoint_age<-meta$older_age-meta$t
    d$start_range<-ranges_list[[meta$a]]; d$end_range<-ranges_list[[meta$b]]
    d$method<-"uniformization_endpoint_bridge_relative_tail_control"
    ctx$audit[[length(ctx$audit)+1L]]<-d
    events
  }
  list(branch=branch,
    reset_audit=function() {ctx$audit<-list(); invisible(NULL)},
    get_audit=function() {if(!length(ctx$audit)) return(data.frame()); do.call(rbind,lapply(ctx$audit,function(x)as.data.frame(x,stringsAsFactors=FALSE)))},
    cache_stats=function() list(Q_count=ctx$q_count,vector_count=ctx$vector_count,max_vector_bytes=ctx$max_vector_bytes),
    kernel=function(Qmat,t,a,b) bgb_ctmc_bridge(Qmat,t,a,b,ctx))
}
