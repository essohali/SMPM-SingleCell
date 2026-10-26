# E1 consistency: run from project root with source('R/01_E1_consistency.R')
if (!exists('smpm_config', mode='function')) source('R/00_simulation_core.R')
set.seed(20261005)
cfg <- smpm_config()
stopifnot(cfg$min_density >= 0)
n_values <- c(200L,500L,1000L)
B <- as.integer(Sys.getenv('SMPM_B', '500'))
stopifnot(B>=2)
depth <- 1000L
outdir <- 'results/tables'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
raw <- vector('list',length(n_values))
for (j in seq_along(n_values)) {
  n <- n_values[j]
  ans <- matrix(NA_real_, B, 7L,
                dimnames=list(NULL,c('lambda1_hat','lambda2_hat','phi1_error',
                                     'phi2_error','angle1','angle2','hs_error')))
  for (b in seq_len(B)) {
    fit <- smpm_fit(smpm_simulate(cfg,n,depth),reference=cfg)
    ans[b,] <- c(fit$lambda,fit$phi_error,fit$angle,fit$hs_error)
    if (b %% 25L==0L || b==B) cat(sprintf('E1 n=%d replication %d/%d\n',n,b,B))
  }
  raw[[j]] <- as.data.frame(ans)
  raw[[j]]$n <- n
}
names(raw) <- paste0('n',n_values)
summary <- do.call(rbind,lapply(seq_along(n_values),function(j) {
  n <- n_values[j]; x <- raw[[j]]
  do.call(rbind,lapply(1:2,function(k) {
    vals <- x[[paste0('lambda',k,'_hat')]]
    data.frame(n=n,component=k,lambda_true=cfg$lambda[k],
               mean_lambda_hat=mean(vals),bias=mean(vals)-cfg$lambda[k],
               rmse=sqrt(mean((vals-cfg$lambda[k])^2)),
               mean_phi_error=mean(x[[paste0('phi',k,'_error')]]),
               mean_angle_rad=mean(x[[paste0('angle',k)]]),
               mean_angle_deg=mean(x[[paste0('angle',k)]])*180/pi,
               mean_HS_error=mean(x$hs_error),
               sqrt_n_HS_error=sqrt(n)*mean(x$hs_error))
  }))
}))
rownames(summary) <- NULL
write.csv(summary,file.path(outdir,'E1_summary.csv'),row.names=FALSE)
saveRDS(raw,file.path(outdir,'E1_raw_results.rds'))
writeLines(capture.output(sessionInfo()),'results/sessionInfo.txt')
print(summary,digits=6,row.names=FALSE)
cat('E1 finished. Monte Carlo results are newly generated, not copied from old tables.\n')
