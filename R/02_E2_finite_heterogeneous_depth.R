# E2: fixed and heterogeneous depth; base R; source from project root.
if (!exists('smpm_config', mode='function')) source('R/00_simulation_core.R')
set.seed(20261009)
B <- as.integer(Sys.getenv('SMPM_B','500'))
stopifnot(B >= 2L)
n <- 500L
mean_depths <- c(50L,200L,1000L)
# The imbalanced proportions mirror the relative sizes of three biological groups.
p_balanced <- rep(1/3,3)
p_imbalanced <- c(325,384,1306)/2015
# Positive NB depth: 1 + NB(mu=target-1,size=dispersion).
# E[D|g,r] = target * group_factor[g] * mark_factor[r] (exact).
nb_size <- 10
mark_factors <- c(.60,.85,1.00,1.20,1.35)
group_factors_raw <- c(.80,1.00,1.20)

# Draw one replicate, retaining group means only. The estimand is unchanged
# because event sampling is conditionally unbiased, even with varying depths.
simulate_E2 <- function(cfg, n, mean_depth, design, nb_size=10) {
  y <- sample.int(3L,n,replace=TRUE,prob=cfg$p)
  counts <- tabulate(y,nbins=3L)
  if (any(counts==0L)) stop('Empty group')
  m <- length(cfg$grid)
  d <- 5L*m
  sums <- matrix(0,3L,d)
  depths <- matrix(NA_integer_,n,5L)
  if (design == 'fixed') {
    depth_means <- matrix(mean_depth,3L,5L)
  } else if (design == 'heterogeneous') {
    gf <- group_factors_raw / sum(cfg$p * group_factors_raw)
    mf <- mark_factors / mean(mark_factors)
    depth_means <- mean_depth * outer(gf,mf)
  } else stop('Unknown design')
  for (i in seq_len(n)) {
    row <- numeric(d)
    for (r in seq_len(5L)) {
      D <- if (design=='fixed') as.integer(mean_depth) else
        1L + rnbinom(1L,size=nb_size,mu=depth_means[y[i],r]-1)
      depths[i,r] <- D
      events <- sort(cfg$inverse[[y[i],r]](runif(D)))
      idx <- ((r-1L)*m+1L):(r*m)
      row[idx] <- findInterval(cfg$grid,events,rightmost.closed=TRUE)/D * cfg$sw[idx]
    }
    sums[y[i],] <- sums[y[i],]+row
  }
  list(y=y,counts=counts,group_means=sweep(sums,1,counts,'/'),
       X=NULL,depths=depths)
}

outdir <- 'results/tables'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
scenarios <- expand.grid(mean_depth=mean_depths,
                         scenario=c('A_fixed_balanced','B_hetero_balanced','B_hetero_imbalanced'),
                         stringsAsFactors=FALSE)
raw <- vector('list',nrow(scenarios))
for (j in seq_len(nrow(scenarios))) {
  s <- scenarios[j,]
  cfg <- smpm_config(p=if (s$scenario=='B_hetero_imbalanced') p_imbalanced else p_balanced)
  design <- if (s$scenario=='A_fixed_balanced') 'fixed' else 'heterogeneous'
  z <- matrix(NA_real_,B,10L,
    dimnames=list(NULL,c('lambda1_hat','lambda2_hat','phi1_error','phi2_error',
                         'angle1','angle2','hs_error','mean_depth_realized',
                         'harmonic_depth_realized','depth_cv_realized')))
  for (b in seq_len(B)) {
    sim <- simulate_E2(cfg,n,s$mean_depth,design,nb_size)
    fit <- smpm_fit(sim,reference=cfg)
    dd <- as.vector(sim$depths)
    z[b,] <- c(fit$lambda,fit$phi_error,fit$angle,fit$hs_error,
              mean(dd),1/mean(1/dd),sd(dd)/mean(dd))
    if (b %% 25L==0L || b==B)
      cat(sprintf('E2 %s meanN=%d replicate %d/%d\n',s$scenario,s$mean_depth,b,B))
  }
  raw[[j]] <- as.data.frame(z)
  raw[[j]]$scenario <- s$scenario
  raw[[j]]$mean_depth <- s$mean_depth
}
names(raw) <- paste(scenarios$scenario,scenarios$mean_depth,sep='_N')
summary <- do.call(rbind,lapply(seq_len(nrow(scenarios)),function(j) {
  s <- scenarios[j,]; x <- raw[[j]]
  cfg <- smpm_config(p=if (s$scenario=='B_hetero_imbalanced') p_imbalanced else p_balanced)
  do.call(rbind,lapply(1:2,function(k) {
    vals <- x[[paste0('lambda',k,'_hat')]]
    data.frame(scenario=s$scenario,mean_depth=s$mean_depth,n=n,k=k,
       lambda_true=cfg$lambda[k],bias=mean(vals)-cfg$lambda[k],
       rmse=sqrt(mean((vals-cfg$lambda[k])^2)),mc_sd=sd(vals),
       mc_variance=var(vals),mean_phi_error=mean(x[[paste0('phi',k,'_error')]]),
       mean_angle_deg=mean(x[[paste0('angle',k)]])*180/pi,
       mean_HS_error=mean(x$hs_error),
       realized_mean_depth=mean(x$mean_depth_realized),
       realized_harmonic_depth=mean(x$harmonic_depth_realized),
       realized_depth_cv=mean(x$depth_cv_realized))
  }))
}))
rownames(summary) <- NULL
# Comparisons of variance only within the SAME phenotype design.
summary$variance_ratio_to_fixed <- NA_real_
for (j in seq_len(nrow(summary))) {
  ss <- summary[j,]
  if (ss$scenario=='B_hetero_balanced') {
    ref <- subset(summary,scenario=='A_fixed_balanced' & mean_depth==ss$mean_depth & k==ss$k)$mc_variance
    summary$variance_ratio_to_fixed[j] <- ss$mc_variance/ref
  }
}
write.csv(summary,file.path(outdir,'E2_summary.csv'),row.names=FALSE)
saveRDS(raw,file.path(outdir,'E2_raw_results.rds'))
writeLines(capture.output(sessionInfo()),'results/sessionInfo_E2.txt')
print(summary,digits=6,row.names=FALSE)
cat('E2 complete. Note: imbalanced-group population eigenpairs are recalculated.\n')
