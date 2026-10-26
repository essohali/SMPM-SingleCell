# E6: localized phenotype-associated signals; base R only.
# Run from repository root via source('run_E6.R').
# Same 3-group between-group operator for all three representations.
# E6 compares representation choices, not independent supervised theories.
source('R/00_simulation_core.R')

B <- as.integer(Sys.getenv('SMPM_E6_B','500'))
n <- as.integer(Sys.getenv('SMPM_E6_N','500'))
M <- as.integer(Sys.getenv('SMPM_E6_M','201'))
seed <- as.integer(Sys.getenv('SMPM_E6_SEED','20261015'))
K_sine <- as.integer(Sys.getenv('SMPM_E6_SINE','3'))
stopifnot(B>=2L,n>=30L,M>=101L,K_sine>=1L)
set.seed(seed)
q <- 5L
nbins <- as.integer(Sys.getenv('SMPM_E6_WINDOWS', '20'))
stopifnot(nbins >= 2L, nbins <= M-1L, (M-1L) %% nbins == 0L)
levels <- c(50L,200L,1000L)
scenarios <- c('broad','narrow','adjacent')

# Smooth, compactly supported bump; exactly zero outside (center-width,center+width).
bump <- function(t, center, width) {
  z <- (t-center)/width
  ans <- numeric(length(t))
  ok <- abs(z)<1
  ans[ok] <- exp(1-1/(1-z[ok]^2))
  ans
}
# Each contrast is a positive-negative localized density perturbation.
# Integrals cancel by normalizing both bumps on the integration grid.
trap_cum <- function(y,t) {
  d <- diff(t)
  c(0,cumsum((head(y,-1L)+tail(y,-1L))*d/2))
}
contrast <- function(t,centers,width) {
  a <- bump(t,centers[1],width); b <- bump(t,centers[2],width)
  ia <- tail(trap_cum(a,t),1); ib <- tail(trap_cum(b,t),1)
  stopifnot(ia>0,ib>0)
  h <- a/ia-b/ib
  h/max(abs(h))  # bounded amplitude; numerical integral remains zero
}
scenario_params <- function(scenario) {
  switch(scenario,
    broad=list(centers=rbind(c(.15,.30),c(.66,.82)),width=.055),
    narrow=list(centers=rbind(c(.15,.30),c(.66,.82)),width=.018),
    adjacent=list(centers=rbind(c(.35,.46),c(.43,.54)),width=.035),
    stop('Unknown scenario'))
}

localized_config <- function(scenario,M=201L,inverse_M=5001L,
                             alpha=.4,beta=.3) {
  t <- seq(0,1,length.out=M)
  ti <- seq(0,1,length.out=inverse_M)
  w <- rep(1/(M-1),M); w[c(1,M)] <- w[c(1,M)]/2
  sw <- rep(sqrt(w),times=q)
  a <- cbind(c(1,1,0,0,0)/sqrt(2),c(0,0,1,-1,0)/sqrt(2))
  theta <- rbind(c(alpha,beta),c(-alpha,beta),c(0,-2*beta))
  prm <- scenario_params(scenario)
  hi <- vapply(1:2,function(k) contrast(ti,prm$centers[k,],prm$width),numeric(inverse_M))
  # Compute CDFs from high-resolution numerical integration, then interpolate.
  Hi <- apply(hi,2,trap_cum,t=ti)
  H <- sapply(1:2,function(k) approx(ti,Hi[,k],xout=t,rule=2)$y)
  h <- sapply(1:2,function(k) approx(ti,hi[,k],xout=t,rule=2)$y)
  inverse <- matrix(vector('list',3L*q),3L,q)
  means <- matrix(0,3L,q*M)
  for(g in 1:3) for(r in 1:q) {
    coeff <- theta[g,]*a[r,]
    dens <- 1+as.vector(hi%*%coeff)
    if(min(dens)<0.05) stop('Nonpositive density')
    Fi <- ti+as.vector(Hi%*%coeff)
    Fi <- Fi/Fi[length(Fi)] # numerical endpoint correction
    if(any(diff(Fi)<=0)) stop('Nonmonotone CDF')
    inverse[[g,r]] <- approxfun(Fi,ti,rule=2,ties='ordered')
    idx <- ((r-1L)*M+1L):(r*M)
    means[g,idx] <- approx(ti,Fi,xout=t,rule=2)$y*sqrt(w)
  }
  p <- rep(1/3,3)
  z <- colSums(means*p)
  A <- sweep(sweep(means,2,z,'-'),1,sqrt(p),'*')
  ee <- eigen(tcrossprod(A),symmetric=TRUE)
  lam <- ee$values[1:2]
  stopifnot(all(lam>1e-12),lam[1]>lam[2]*1.01)
  phi <- sweep(t(A)%*%ee$vectors[,1:2,drop=FALSE],2,sqrt(lam),'/')
  list(grid=t,w=w,sw=sw,p=p,means=means,phi=phi,lambda=lam,
       inverse=inverse,scenario=scenario,h=h,H=H,theta=theta,a=a,
       params=prm)
}

# Weighted projection, applied separately to each component of group means.
make_smoother <- function(t,w,Kbasis=3L) {
  basis <- cbind(t,sapply(seq_len(Kbasis),function(j) sin(2*pi*j*t)))
  Wb <- sweep(basis,1,w,'*')
  basis %*% solve(crossprod(basis,Wb),t(Wb))
}
make_binner <- function(t,nb=10L) {
  m <- length(t); knots <- seq(0,1,length.out=nb+1L)
  ix <- vapply(knots,function(z) which.min(abs(t-z)),integer(1))
  P <- matrix(0,m,m)
  for(j in seq_len(m)) {
    a <- max(which(knots<=t[j])); a <- min(a,nb)
    v <- (t[j]-knots[a])/(knots[a+1]-knots[a])
    P[j,ix[a]] <- P[j,ix[a]]+1-v
    P[j,ix[a+1]] <- P[j,ix[a+1]]+v
  }
  P
}
transform_means <- function(group_means,map,w) {
  out <- group_means; m <- length(w)
  for(g in seq_len(nrow(out))) for(r in 1:q) {
    ix <- ((r-1L)*m+1L):(r*m)
    out[g,ix] <- as.vector(map%*%(out[g,ix]/sqrt(w)))*sqrt(w)
  }
  out
}
fit_means <- function(means,counts) {
  p <- counts/sum(counts)
  z <- colSums(means*p)
  A <- sweep(sweep(means,2,z,'-'),1,sqrt(p),'*')
  e <- eigen(tcrossprod(A),symmetric=TRUE)
  val <- e$values[1:2]
  if(any(val<1e-12)) stop('Degenerate supervised eigenspace')
  sweep(t(A)%*%e$vectors[,1:2,drop=FALSE],2,sqrt(val),'/')
}

# Vector-measure total variation of the two-direction increment vector.
# Invariant under any orthogonal change of basis of the estimated subspace.

window_tv <- function(phi, w, nb = 20L) {
  m <- length(w)
  stopifnot(
    nrow(phi) == q * m,
    ncol(phi) == 2L
  )
  
  stopifnot((m-1L) %% nb == 0L)
  edges <- round(seq(1, m, length.out = nb + 1L))
  scores <- numeric(q * nb)
  
  for (r in seq_len(q)) {
    ix <- ((r - 1L) * m + 1L):(r * m)
    
    # Return to unweighted functional coordinates
    f <- sweep(
      phi[ix, , drop = FALSE],
      1, sqrt(w), "/"
    )
    
    for (j in seq_len(nb)) {
      steps <- seq.int(edges[j], edges[j + 1L])
      increments <- diff(f[steps, , drop = FALSE])
      
      # Rotation-invariant vector-measure variation
      scores[(r - 1L) * nb + j] <-
        sum(sqrt(rowSums(increments^2)))
    }
  }
  
  scores
}

# Oracle active windows are determined from density-contrast support, not
# by taking an arbitrary top fraction of estimated scores.
true_windows <- function(cfg,nb=20L) {
  # Analytic support: each bump is nonzero on (center-width, center+width).
  # The positive and negative bumps have disjoint supports in these scenarios.
  t <- cfg$grid; m <- length(t)
  stopifnot((m-1L) %% nb == 0L)
  edges <- round(seq(1, m, length.out=nb+1L))
  out <- logical(q*nb)
  for (r in seq_len(q)) for (j in seq_len(nb)) {
    lo <- t[edges[j]]; hi <- t[edges[j+1L]]
    relevant <- which(abs(cfg$a[r,]) > 0)
    active <- FALSE
    for (k in relevant) for (center in cfg$params$centers[k,]) {
      support_lo <- center - cfg$params$width
      support_hi <- center + cfg$params$width
      # Strict positive-length intersection; touching endpoints is not active.
      if (max(lo, support_lo) < min(hi, support_hi)) active <- TRUE
    }
    out[(r-1L)*nb+j] <- active
  }
  out
}

# Average precision (AP) for binary window labels; higher is better.
avg_precision <- function(score,label) {
  ord <- order(score,decreasing=TRUE)
  yy <- as.integer(label[ord]); k <- sum(yy)
  if(k==0L) return(NA_real_)
  sum(cumsum(yy)[yy==1L]/which(yy==1L))/k
}
# Top-K F1, K = number of oracle active windows. This is an oracle-budget
# evaluation, NOT a data-driven region-selection procedure.
topk_f1 <- function(score,label) {
  k <- sum(label)
  if(k==0L) return(NA_real_)
  selected <- order(score,decreasing=TRUE)[seq_len(k)]
  sum(label[selected])/k
}

one_replication <- function(cfg,depth,S,P,replication) {
  dat <- smpm_simulate(cfg,n=n,depth=depth,store_profiles=FALSE)
  candidates <- list(
    SMPM=dat$group_means,
    Smoothed_supervised=transform_means(dat$group_means,S,cfg$w),
    Binned_supervised=transform_means(dat$group_means,P,cfg$w)
  )
  truth <- true_windows(cfg,nbins)
  truth_tv <- window_tv(cfg$phi,cfg$w,nbins)
  do.call(rbind,lapply(names(candidates),function(method) {
    est <- fit_means(candidates[[method]],dat$counts)
    sv <- pmin(1,pmax(0,svd(crossprod(cfg$phi,est),nu=0,nv=0)$d))
    score <- window_tv(est,cfg$w,nbins)
    data.frame(scenario=cfg$scenario,depth=depth,replication=replication,
      method=method,subspace_min_cos=min(sv),
      subspace_sin_max=sqrt(max(0,1-min(sv)^2)),
      window_average_precision=avg_precision(score,truth),
      window_topK_f1=topk_f1(score,truth),
      window_spearman=suppressWarnings(cor(score,truth_tv,method='spearman')),
      active_windows=sum(truth))
  }))
}

# Self-tests: exact partition and orthogonal invariance on the truth.
for (scenario in scenarios) {
  cfg_test <- localized_config(scenario,M=M)
  Q <- matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),2L,2L)
  stopifnot(max(abs(window_tv(cfg_test$phi,cfg_test$w,nbins) -
                    window_tv(cfg_test$phi %*% Q,cfg_test$w,nbins))) < 1e-9)
  stopifnot(any(true_windows(cfg_test,nbins)), !all(true_windows(cfg_test,nbins)))
}

out <- vector('list',length(scenarios)*length(levels))
ix <- 0L
for(scenario in scenarios) {
  cfg <- localized_config(scenario,M=M)
  S <- make_smoother(cfg$grid,cfg$w,K_sine)
  P <- make_binner(cfg$grid,10L)
  cat(sprintf('E6 scenario=%s: eigenvalues %.6g, %.6g; active windows=%d/%d\n',
              scenario,cfg$lambda[1],cfg$lambda[2],
              sum(true_windows(cfg,nbins)),q*nbins))
  for(depth in levels) {
    ix <- ix+1L
    out[[ix]] <- do.call(rbind,lapply(seq_len(B),function(b)
      one_replication(cfg,depth,S,P,b)))
    cat(sprintf('E6 %s depth=%d completed (%d replications)\n',scenario,depth,B))
  }
}
raw <- do.call(rbind,out)
row.names(raw) <- NULL
metrics <- c('subspace_min_cos','subspace_sin_max','window_average_precision',
             'window_topK_f1','window_spearman')
groups <- split(raw,list(raw$scenario,raw$depth,raw$method),drop=TRUE)
summary <- do.call(rbind,lapply(groups,function(z) {
  do.call(rbind,lapply(metrics,function(metric) {
    x <- z[[metric]]; good <- is.finite(x)
    data.frame(scenario=z$scenario[1],depth=z$depth[1],method=z$method[1],
               metric=metric,mean=if(any(good)) mean(x[good]) else NA_real_,
               sd=if(sum(good)>1) sd(x[good]) else NA_real_,
               mc_se=if(sum(good)>1) sd(x[good])/sqrt(sum(good)) else NA_real_,
               n_valid=sum(good),active_windows=z$active_windows[1])
  }))
}))
# Paired differences because all methods use the same replicate within
# each scenario/depth; positive difference favors SMPM for these metrics.
paired <- do.call(rbind,lapply(split(raw,list(raw$scenario,raw$depth),drop=TRUE),function(z) {
  do.call(rbind,lapply(c('Smoothed_supervised','Binned_supervised'),function(competitor) {
    aa <- z[z$method=='SMPM',]; bb <- z[z$method==competitor,]
    aa <- aa[order(aa$replication),]; bb <- bb[order(bb$replication),]
    stopifnot(identical(aa$replication,bb$replication))
    do.call(rbind,lapply(c('subspace_min_cos','window_average_precision','window_topK_f1'),function(metric) {
      delta <- aa[[metric]]-bb[[metric]]
      data.frame(scenario=aa$scenario[1],depth=aa$depth[1],competitor=competitor,
                 metric=metric,mean_difference=mean(delta),
                 mc_se_difference=sd(delta)/sqrt(length(delta)))
    }))
  }))
}))
dir.create('results/tables',recursive=TRUE,showWarnings=FALSE)
# Distinct filenames preserve the earlier E6 results and the 20/40 comparison.
suffix <- sprintf('_analytic_W%d',nbins)
write.csv(raw,paste0('results/tables/E6_replications',suffix,'.csv'),row.names=FALSE)
write.csv(summary,paste0('results/tables/E6_summary',suffix,'.csv'),row.names=FALSE)
write.csv(paired,paste0('results/tables/E6_paired_differences',suffix,'.csv'),row.names=FALSE)
saveRDS(list(settings=list(B=B,n=n,M=M,seed=seed,depths=levels,
                           scenarios=scenarios,sine_basis=K_sine,windows=nbins,
                           truth_definition='analytic positive-length support overlap'),
             raw=raw,summary=summary,paired=paired),
        paste0('results/tables/E6_results',suffix,'.rds'))
cat(sprintf('E6 complete, %d windows. Files: results/tables/*%s*\n',nbins,suffix))
