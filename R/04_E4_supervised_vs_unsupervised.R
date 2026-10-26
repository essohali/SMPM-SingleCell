# E4 | Shared SMPM simulation engine and weighted L2 geometry
source('R/00_simulation_core.R')
B <- as.integer(Sys.getenv('SMPM_E4_B','500'))
n <- as.integer(Sys.getenv('SMPM_E4_N','500'))
depth <- as.integer(Sys.getenv('SMPM_E4_DEPTH','1000'))
M <- as.integer(Sys.getenv('SMPM_E4_M','201'))
rhos <- c(0,.30,.50,.70)
seed <- as.integer(Sys.getenv('SMPM_E4_SEED','20261010'))
stopifnot(B>=1,n>=10,depth>=1,M>=11,all(rhos*sqrt(2)<1))
set.seed(seed)
cfg <- smpm_config(M=M, p=rep(1/3,3))
# All profiles are in weighted coordinates: x_(r,j)=sqrt(w_j)*F_r(t_j).
# v0(t)=e5 H1(t), nuisance affects only component 5.
t <- cfg$grid
H1 <- function(x) sqrt(2)*sin(2*pi*x)/(2*pi)
u1 <- function(x) sqrt(2)*cos(2*pi*x)
q <- 5L; m <- length(t)
v0 <- numeric(q*m)
v0[(4L*m+1L):(5L*m)] <- H1(t)*sqrt(cfg$w)
phi0 <- v0/sqrt(sum(v0^2))
stopifnot(max(abs(crossprod(cfg$phi,phi0)))<1e-10)
# Reuse EXACT inverse-CDF construction of the shared engine (5001 grid),
# modifying only mark 5, whose baseline density is uniform in all groups.
invgrid <- seq(0,1,length.out=5001L)
make_inverse_5 <- function(rho,xi) {
  F <- invgrid + rho*xi*H1(invgrid)
  stopifnot(all(diff(F)>0), min(1+rho*xi*u1(invgrid))>0)
  approxfun(F,invgrid,ties='ordered',rule=2)
}
# The shared simulator draws Y, group-specific events and stores X;
# we replace component-five observations with the nuisance-perturbed draw.
# Other components, normalization and trapezoidal weights are unchanged.
simulate_E4 <- function(rho) {
  dat <- smpm_simulate(cfg,n=n,depth=depth,store_profiles=TRUE)
  xi <- sample(c(-1,1),n,replace=TRUE)
  inv <- list(make_inverse_5(rho,-1),make_inverse_5(rho,1))
  idx <- (4L*m+1L):(5L*m)
  for (i in seq_len(n)) {
    events <- sort(inv[[if(xi[i]<0) 1L else 2L]](runif(depth)))
    dat$X[i,idx] <- findInterval(t,events,rightmost.closed=TRUE)/depth*cfg$sw[idx]
  }
  # Recompute group means from the modified, weighted cell profiles.
  dat$group_means <- do.call(rbind,lapply(1:3,function(g)
    colMeans(dat$X[dat$y==g,,drop=FALSE])))
  dat
}
# Top 3 right singular vectors of centered weighted profiles.
# RSpectra strongly recommended for full runs; base R fallback uses
# a n x n Gram eigen-decomposition rather than huge feature covariance.
fit_unsupervised <- function(X) {
  Xc <- sweep(X,2,colMeans(X),'-')
  if (requireNamespace('RSpectra',quietly=TRUE)) {
    s <- RSpectra::svds(Xc,k=3,nu=0,nv=3)
    ord <- order(s$d,decreasing=TRUE)
    V <- s$v[,ord,drop=FALSE]
    values <- s$d[ord]^2/(nrow(X)-1)
  } else {
    K <- tcrossprod(Xc)
    e <- eigen(K,symmetric=TRUE)
    values <- e$values[1:3]/(nrow(X)-1)
    V <- crossprod(Xc,e$vectors[,1:3,drop=FALSE])
    V <- sweep(V,2,sqrt(pmax(e$values[1:3],.Machine$double.eps)),'/')
  }
  V <- sweep(V,2,sqrt(colSums(V^2)),'/')
  list(phi=V,lambda=values)
}
min_cos <- function(A,B) min(pmin(1,svd(crossprod(A,B),nu=0,nv=0)$d))
measure <- function(dat) {
  S <- smpm_fit(dat)
  U <- fit_unsupervised(dat$X)
  ref <- cfg$phi
  c(S_phi1=abs(drop(crossprod(ref[,1],S$phi[,1]))),
    U_phi1=abs(drop(crossprod(ref[,1],U$phi[,1]))),
    S_subspace_min=min_cos(ref,S$phi[,1:2,drop=FALSE]),
    U_top2_subspace_min=min_cos(ref,U$phi[,1:2,drop=FALSE]),
    U_top3_subspace_min=min_cos(ref,U$phi[,1:3,drop=FALSE]),
    U_first_nuisance=abs(drop(crossprod(phi0,U$phi[,1]))),
    S_lambda1=S$lambda[1],S_lambda2=S$lambda[2],
    U_lambda1=U$lambda[1])
}
dir.create('results/tables',recursive=TRUE,showWarnings=FALSE)
dir.create('results/figures',recursive=TRUE,showWarnings=FALSE)
raw_list <- vector('list',length(rhos))
for (j in seq_along(rhos)) {
  rho <- rhos[j]
  vals <- replicate(B,measure(simulate_E4(rho)))
  raw_list[[j]] <- data.frame(rho=rho,replication=seq_len(B),t(vals),check.names=FALSE)
  cat(sprintf('E4 rho=%.2f: %d replications completed\n',rho,B))
}
raw <- do.call(rbind,raw_list)
rownames(raw) <- NULL
write.csv(raw,'results/tables/E4_replications.csv',row.names=FALSE)
metrics <- setdiff(names(raw),c('rho','replication'))
summary_rows <- lapply(split(raw, raw$rho), function(z) {
  do.call(rbind, lapply(metrics, function(k) {
    data.frame(
      rho = z$rho[1], metric = k,
      mean = mean(z[[k]]), sd = sd(z[[k]]),
      q025 = unname(quantile(z[[k]], 0.025)),
      q975 = unname(quantile(z[[k]], 0.975))
    )
  }))
})
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
write.csv(summary,'results/tables/E4_summary.csv',row.names=FALSE)
pop <- data.frame(rho=rhos,lambda_pheno1=cfg$lambda[1],
  lambda_pheno2=cfg$lambda[2],lambda_nuisance=rhos^2/(4*pi^2))
write.csv(pop,'results/tables/E4_population_eigenvalues.csv',row.names=FALSE)
saveRDS(list(config=list(B=B,n=n,depth=depth,M=M,rhos=rhos,seed=seed),
             raw=raw,summary=summary,population=pop),
        'results/tables/E4_raw_results.rds')
# Figure: same scientific panels as the original E4.
getmean <- function(metric_name) {
  z <- summary[summary$metric == metric_name, , drop = FALSE]
  z <- z[order(z$rho), , drop = FALSE]
  
  if (nrow(z) != length(rhos)) {
    stop("Nombre de valeurs incorrect pour : ", metric_name)
  }
  
  z$mean
}

pdf('results/figures/Fig_E4_supervision.pdf',width=9,height=4.2,
    family='Helvetica',useDingbats=FALSE)
par(mfrow=c(1,2),mar=c(4.2,4.3,1.5,.6),mgp=c(2.5,.65,0),las=1)
plot(rhos,getmean('S_phi1'),type='b',pch=16,lwd=2,col='#176B87',
     ylim=c(0,1.03),xaxt='n',xlab=expression('Nuisance amplitude '*rho),
     ylab='Alignment with first phenotype direction')
axis(1,at=rhos); lines(rhos,getmean('U_phi1'),type='b',pch=17,lwd=2,col='#C05B39')
abline(v=sqrt(4*pi^2*cfg$lambda[1]),lty=3,col='gray50')
legend('bottomleft',c('SMPM','Unsupervised PC1'),lty=1,pch=c(16,17),
       col=c('#176B87','#C05B39'),bty='n',cex=.83)
mtext('(a)',side=3,adj=0,font=2)
plot(rhos,getmean('S_subspace_min'),type='b',pch=16,lwd=2,col='#176B87',
     ylim=c(0,1.03),xaxt='n',xlab=expression('Nuisance amplitude '*rho),
     ylab='Minimum canonical alignment')
axis(1,at=rhos)
lines(rhos,getmean('U_top2_subspace_min'),type='b',pch=17,lwd=2,col='#C05B39')
lines(rhos,getmean('U_top3_subspace_min'),type='b',pch=15,lwd=2,lty=2,col='#666666')
abline(v=sqrt(4*pi^2*cfg$lambda[2]),lty=3,col='gray50')
legend('bottomleft',c('SMPM top 2','Unsupervised top 2','Unsupervised top 3'),
       lty=c(1,1,2),pch=c(16,17,15),col=c('#176B87','#C05B39','#666666'),bty='n',cex=.78)
mtext('(b)',side=3,adj=0,font=2)
dev.off()
cat('E4 complete. Outputs in results/tables and results/figures.\n')
