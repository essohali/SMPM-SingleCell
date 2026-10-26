# E5: transparent, reproducible comparison. Run from repository root.
source('R/00_simulation_core.R')
B <- as.integer(Sys.getenv('SMPM_E5_B','500'))
n <- as.integer(Sys.getenv('SMPM_E5_N','500'))
M <- as.integer(Sys.getenv('SMPM_E5_M','201'))
D <- as.integer(Sys.getenv('SMPM_E5_DEPTH','200'))
seed <- as.integer(Sys.getenv('SMPM_E5_SEED','20261015'))
set.seed(seed)
stopifnot(B>=1,n>=10,M>=21,D>=1)
cfg <- smpm_config(M=M)
q <- 5L; K <- 2L; m <- length(cfg$grid); w <- cfg$w
# Smooth cumulative profiles by projecting on a low-dimensional sine basis,
# retaining the linear trend. This is a regularized supervised competitor,
# NOT a verbatim implementation of Nie et al. (2018).
make_smoother <- function(Kbasis=3L) {
  t <- cfg$grid
  basis <- cbind(t, sapply(seq_len(Kbasis),function(j) sin(2*pi*j*t)))
  # Weighted least-squares projection, applied separately to each mark.
  W <- diag(w)
  solve(crossprod(basis,W%*%basis),crossprod(basis,W)) -> coef
  basis%*%coef
}
S <- make_smoother(3L)
# Coarse-bin competitor: retain only CDF values on bin boundaries and
# reconstruct by linear interpolation. This is a discretized supervised
# contrast method, not a direct published sFPCA implementation.
make_binner <- function(nbins=10L) {
  t <- cfg$grid; knots <- seq(0,1,length.out=nbins+1L)
  ix <- vapply(knots,function(z) which.min(abs(t-z)),integer(1))
  P <- matrix(0,m,m)
  for (j in seq_len(m)) {
    a <- max(which(knots<=t[j])); if(a>=length(knots)) a<-length(knots)-1L
    v <- (t[j]-knots[a])/(knots[a+1L]-knots[a])
    P[j,ix[a]] <- P[j,ix[a]]+(1-v)
    P[j,ix[a+1L]] <- P[j,ix[a+1L]]+v
  }
  P
}
P <- make_binner(10L)
transform_group_means <- function(group_means, map) {
  out <- group_means
  for (g in seq_len(nrow(out))) for (r in 1:q) {
    ix <- ((r-1L)*m+1L):(r*m)
    raw <- out[g,ix]/sqrt(w)
    out[g,ix] <- as.vector(map%*%raw)*sqrt(w)
  }
  out
}
fit_means <- function(means, counts) {
  pp <- counts/sum(counts)
  z <- colSums(means*pp)
  A <- sweep(sweep(means,2,z,'-'),1,sqrt(pp),'*')
  e <- eigen(tcrossprod(A),symmetric=TRUE)
  val <- pmax(e$values[1:2],0)
  if(any(val<1e-14)) stop('Degenerate fit')
  phi <- t(A)%*%e$vectors[,1:2,drop=FALSE]
  sweep(phi,2,sqrt(val),'/')
}
mincos <- function(A,B) min(pmin(1,svd(crossprod(A,B),nu=0,nv=0)$d))
# Window scores computed from signed-measure increments. Smooth using
# fixed windows; no pointwise derivative of noisy empirical CDFs.
window_scores <- function(phi,nbins=20L) {
  t <- cfg$grid; knots <- seq(0,1,length.out=nbins+1L)
  ix <- vapply(knots,function(z) which.min(abs(t-z)),integer(1))
  scores <- numeric(q*nbins)
  for(r in 1:q) {
    rr <- ((r-1L)*m+1L):(r*m)
    f <- phi[rr]/sqrt(w)
    scores[((r-1L)*nbins+1L):(r*nbins)] <- abs(diff(f[ix]))
  }
  scores
}
# Compare ranking of windows, avoiding arbitrary binary thresholding.
# Spearman correlations undefined for all-tie vectors are set to NA.
region_cor <- function(est,truth) {
  a <- window_scores(est); b <- window_scores(truth)
  suppressWarnings(cor(a,b,method='spearman'))
}
# Fixed true positive windows = top 20% by oracle signed-measure mass.
region_f1 <- function(est,truth) {
  a <- window_scores(est); b <- window_scores(truth)
  k <- max(1L,ceiling(length(b)*.20))
  ia <- order(a,decreasing=TRUE)[seq_len(k)]
  ib <- order(b,decreasing=TRUE)[seq_len(k)]
  length(intersect(ia,ib))/k  # precision=recall=F1 at fixed top-k
}
measure <- function(depth) {
  dat <- smpm_simulate(cfg,n=n,depth=depth,store_profiles=FALSE)
  candidates <- list(
    SMPM=dat$group_means,
    Smoothed_supervised=transform_group_means(dat$group_means,S),
    Binned_supervised=transform_group_means(dat$group_means,P)
  )
  do.call(rbind,lapply(names(candidates),function(method) {
    est <- fit_means(candidates[[method]],dat$counts)
    # For nearly equal eigenvalues, subspace is more stable than components.
    data.frame(method=method,depth=depth,
      subspace_min_cos=mincos(cfg$phi,est),
      subspace_sin_max=sqrt(max(0,1-mincos(cfg$phi,est)^2)),
      phi1_abs_cos=abs(sum(cfg$phi[,1]*est[,1])),
      phi2_abs_cos=abs(sum(cfg$phi[,2]*est[,2])),
      window_rank_cor1=region_cor(est[,1],cfg$phi[,1]),
      window_rank_cor2=region_cor(est[,2],cfg$phi[,2]),
      window_top20_f1_1=region_f1(est[,1],cfg$phi[,1]),
      window_top20_f1_2=region_f1(est[,2],cfg$phi[,2]))
  }))
}
levels <- c(50L,200L,1000L)
raw <- do.call(rbind,lapply(levels,function(d) {
  out <- lapply(seq_len(B),function(b) {
    x <- measure(d); x$replication <- b; x
  })
  cat(sprintf('E5 depth=%d completed (%d replications)\n',d,B))
  do.call(rbind,out)
}))
row.names(raw) <- NULL
metric_names <- setdiff(names(raw),c('method','depth','replication'))
summary <- do.call(rbind,lapply(split(raw,list(raw$method,raw$depth),drop=TRUE),function(z) {
  do.call(rbind,lapply(metric_names,function(metric) {
    x <- z[[metric]]
    data.frame(method=z$method[1],depth=z$depth[1],metric=metric,
      mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),
      mc_se=sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x))))
  }))
}))
dir.create('results/tables',recursive=TRUE,showWarnings=FALSE)
write.csv(raw,'results/tables/E5_replications.csv',row.names=FALSE)
write.csv(summary,'results/tables/E5_summary.csv',row.names=FALSE)
saveRDS(list(settings=list(B=B,n=n,M=M,seed=seed,depths=levels),raw=raw,summary=summary),
        'results/tables/E5_results.rds')
cat('E5 complete. Compare means and Monte Carlo SE in results/tables/E5_summary.csv\n')
