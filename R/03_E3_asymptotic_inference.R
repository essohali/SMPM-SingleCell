# E3: feasible influence-function inference for the two positive eigenvalues.
# Run from repository root; base R only.
if (!exists('smpm_config', mode='function')) source('R/00_simulation_core.R')
set.seed(as.integer(Sys.getenv('SMPM_E3_SEED', '20261010')))
B <- as.integer(Sys.getenv('SMPM_E3_B', '500'))
stopifnot(B >= 2L)
ns <- c(200L, 500L, 1000L)
depths <- c(50L, 200L, 1000L)
alpha_ci <- .05
zcrit <- qnorm(1-alpha_ci/2)
nb_size <- 10
mark_factors <- c(.60,.85,1.00,1.20,1.35)
group_factors_raw <- c(.80,1.00,1.20)
p_bal <- rep(1/3,3)
p_imb <- c(325,384,1306)/2015
# E3 heterogeneity sensitivity: balanced n=500 at all three mean depths.
# Set SMPM_E3_IMBALANCED=0 to omit the imbalanced phenotype scenario.
include_imbalanced <- Sys.getenv('SMPM_E3_IMBALANCED','1') == '1'

simulate_e3_hetero <- function(cfg,n,mean_depth) {
  y <- sample.int(3L,n,replace=TRUE,prob=cfg$p)
  counts <- tabulate(y,nbins=3L)
  if (any(counts==0L)) stop('Empty phenotype group')
  m <- length(cfg$grid); d <- 5L*m
  X <- matrix(0,n,d)
  sums <- matrix(0,3L,d)
  gf <- group_factors_raw/sum(cfg$p*group_factors_raw)
  mf <- mark_factors/mean(mark_factors)
  target <- mean_depth*outer(gf,mf)
  if (any(target <= 1)) stop('Invalid depth mean')
  dd <- matrix(NA_integer_,n,5L)
  for (i in seq_len(n)) for (r in seq_len(5L)) {
    D <- 1L+rnbinom(1L,size=nb_size,mu=target[y[i],r]-1)
    dd[i,r] <- D
    events <- sort(cfg$inverse[[y[i],r]](runif(D)))
    idx <- ((r-1L)*m+1L):(r*m)
    X[i,idx] <- findInterval(cfg$grid,events,rightmost.closed=TRUE)/D*cfg$sw[idx]
  }
  for (g in seq_len(3L)) sums[g,] <- colSums(X[y==g,,drop=FALSE])
  list(y=y,counts=counts,group_means=sweep(sums,1,counts,'/'),
       X=X,depths=dd)
}

scenarios <- expand.grid(n=ns,mean_depth=depths,
                         scenario='fixed_balanced',stringsAsFactors=FALSE)
scenarios <- rbind(scenarios,
  data.frame(n=500L,mean_depth=depths,scenario='hetero_balanced'))
if (include_imbalanced) scenarios <- rbind(scenarios,
  data.frame(n=500L,mean_depth=depths,scenario='hetero_imbalanced'))
rownames(scenarios) <- NULL
outdir <- 'results/tables'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
raw <- vector('list',nrow(scenarios))
summary_rows <- vector('list',nrow(scenarios)*2L)
for (j in seq_len(nrow(scenarios))) {
  sc <- scenarios[j,]
  cfg <- smpm_config(p=if(sc$scenario=='hetero_imbalanced') p_imb else p_bal)
  z <- matrix(NA_real_,nrow=B,ncol=15L)
  colnames(z) <- c('lambda1','lambda2','sigma2_1','sigma2_2',
                   'ci_lo1','ci_lo2','ci_hi1','ci_hi2',
                   'cover1','cover2','length1','length2',
                   'eigengap','relative_gap','realized_mean_depth')
  for (b in seq_len(B)) {
    sim <- if(sc$scenario=='fixed_balanced') {
      smpm_simulate(cfg,sc$n,depth=sc$mean_depth,store_profiles=TRUE)
    } else simulate_e3_hetero(cfg,sc$n,sc$mean_depth)
    fit <- smpm_fit(sim,reference=NULL)
    if (is.null(fit$sigma2) || any(!is.finite(fit$sigma2)) ||
        any(fit$sigma2 < 0)) stop('Invalid influence-function variance')
    se <- sqrt(fit$sigma2/sc$n)
    lo <- fit$lambda-zcrit*se; hi <- fit$lambda+zcrit*se
    cover <- as.numeric(lo <= cfg$lambda & cfg$lambda <= hi)
    z[b,] <- c(fit$lambda,fit$sigma2,lo,hi,cover,hi-lo,
              diff(rev(fit$lambda)),
              (fit$lambda[1]-fit$lambda[2])/fit$lambda[1],
              if (is.null(sim$depths)) sc$mean_depth else mean(sim$depths))
    if (b %% 25L==0L || b==B) cat(sprintf('E3 %s n=%d D=%d %d/%d\n',
          sc$scenario,sc$n,sc$mean_depth,b,B))
  }
  raw[[j]] <- data.frame(z,replication=seq_len(B),
     scenario=sc$scenario,n=sc$n,mean_depth=sc$mean_depth)
  for (k in 1:2) {
    lambda_hat <- z[,paste0('lambda',k)]
    sig2 <- z[,paste0('sigma2_',k)]
    mc_scaled <- sc$n*var(lambda_hat)
    plugin <- mean(sig2)
    coverage <- mean(z[,paste0('cover',k)])
    summary_rows[[(j-1L)*2L+k]] <- data.frame(
      scenario=sc$scenario,n=sc$n,mean_depth=sc$mean_depth,k=k,
      B=B,lambda_true=cfg$lambda[k],
      bias=mean(lambda_hat)-cfg$lambda[k],
      rmse=sqrt(mean((lambda_hat-cfg$lambda[k])^2)),
      n_mc_variance=mc_scaled,mean_plugin_variance=plugin,
      mc_over_plugin=if(plugin>0) mc_scaled/plugin else NA_real_,
      coverage_95=coverage,
      coverage_mc_se=sqrt(coverage*(1-coverage)/B),
      mean_ci_length=mean(z[,paste0('length',k)]),
      mean_relative_eigengap=mean(z[,'relative_gap']),
      mean_depth_realized=mean(z[,'realized_mean_depth']),
      near_tie_fraction=mean(z[,'relative_gap'] < .10))
  }
}
summary <- do.call(rbind,summary_rows)
rownames(summary) <- NULL
names(raw) <- paste(scenarios$scenario,scenarios$n,scenarios$mean_depth,sep='_')
write.csv(summary,file.path(outdir,'E3_summary.csv'),row.names=FALSE)
saveRDS(raw,file.path(outdir,'E3_raw_results.rds'))
writeLines(capture.output(sessionInfo()),'results/sessionInfo_E3.txt')
print(summary,digits=6,row.names=FALSE)
cat('E3 complete. CIs use only estimated eigenpairs and cell-level profiles.\n')
