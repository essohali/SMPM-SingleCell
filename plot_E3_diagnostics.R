# Diagnostics E3 : aucune simulation n'est relancee.
raw <- readRDS('results/tables/E3_raw_results.rds')
stopifnot(is.list(raw), length(raw) > 0L)

lambda_bal <- c(0.4^2/(6*pi^2), 0.3^2/(8*pi^2))
lambda_imb <- c(0.001449346, 0.001143000)

build_one <- function(z) {
  required <- c('lambda1','lambda2','sigma2_1','sigma2_2','scenario','n','mean_depth')
  stopifnot(is.data.frame(z), all(required %in% names(z)))
  scenario <- unique(z$scenario)
  stopifnot(length(scenario) == 1L)
  target <- if (scenario == 'hetero_imbalanced') lambda_imb else lambda_bal
  data.frame(
    scenario = rep(z$scenario, each = 2L),
    n = rep(z$n, each = 2L),
    mean_depth = rep(z$mean_depth, each = 2L),
    replication = rep(if ('replication' %in% names(z)) z$replication else seq_len(nrow(z)), each = 2L),
    k = rep(1:2, times = nrow(z)),
    lambda_hat = as.vector(t(as.matrix(z[c('lambda1','lambda2')]))),
    sigma2 = as.vector(t(as.matrix(z[c('sigma2_1','sigma2_2')]))),
    lambda_true = rep(target, times = nrow(z))
  )
}
all_results <- do.call(rbind, lapply(raw, build_one))
all_results$T <- with(all_results, sqrt(n)*(lambda_hat-lambda_true)/sqrt(sigma2))
if (any(!is.finite(all_results$T))) stop('Statistiques studentisees non finies : verifier sigma2.')

dir.create('results/figures', recursive = TRUE, showWarnings = FALSE)

select_T <- function(scenario, n, depth, k) {
  # Eviter la capture de variables dans subset() : selection explicite.
  idx <- all_results$scenario == scenario & all_results$n == n &
    all_results$mean_depth == depth & all_results$k == k
  t <- all_results$T[idx]
  if (length(t) < 30L) stop('Moins de 30 observations pour ',scenario,' n=',n,' D=',depth,' k=',k)
  t
}

# 4 panneaux : cas difficile et cas favorable, chacun pour k=1,2
configs <- data.frame(scenario = 'fixed_balanced',
                      n = c(200,200,1000,1000),
                      depth = c(50,50,1000,1000),
                      k = c(1,2,1,2))

pdf('results/figures/E3_density.pdf', width=9.5, height=7.5, family='Helvetica')
par(mfrow=c(2,2), mar=c(4.2,4.2,3,1), oma=c(0,0,0,0))
for (i in seq_len(nrow(configs))) {
  cc <- configs[i,]
  tt <- select_T(cc$scenario, cc$n, cc$depth, cc$k)
  dd <- density(tt, from=-4.5, to=4.5)
  plot(dd, xlim=c(-4.5,4.5), ylim=c(0,max(dd$y,dnorm(0))*1.12),
       lwd=2.2, col='#2166AC', xlab='Studentized statistic', ylab='Density',
       main=sprintf('(%s) n = %d, D = %d, k = %d',letters[i],cc$n,cc$depth,cc$k))
  curve(dnorm(x), from=-4.5,to=4.5,add=TRUE,col='#B2182B',lwd=2,lty=2)
  abline(v=0,col='gray65',lty=3)
  if (i==1) legend('topright',c('Monte Carlo KDE','Standard normal'),
                   col=c('#2166AC','#B2182B'),lty=c(1,2),lwd=2,bty='n',cex=.78)
}
dev.off()

pdf('results/figures/E3_qqplots.pdf', width=9.5, height=7.5, family='Helvetica')
par(mfrow=c(2,2),mar=c(4.2,4.2,3,1))
for (i in seq_len(nrow(configs))) {
  cc <- configs[i,]
  tt <- select_T(cc$scenario,cc$n,cc$depth,cc$k)
  qqnorm(tt, main=sprintf('(%s) n = %d, D = %d, k = %d',letters[i],cc$n,cc$depth,cc$k),
         xlab='Theoretical N(0,1) quantiles', ylab='Observed studentized quantiles',
         pch=16,cex=.42,col='#2166AC')
  abline(a=0,b=1,col='#B2182B',lty=2,lwd=2)
}
dev.off()

split_index <- split(seq_len(nrow(all_results)),
                     interaction(all_results$scenario,all_results$n,all_results$mean_depth,
                                 all_results$k,drop=TRUE))
diag <- do.call(rbind,lapply(split_index,function(ii) {
  d <- all_results[ii,]; t <- d$T
  data.frame(scenario=d$scenario[1],n=d$n[1],mean_depth=d$mean_depth[1],k=d$k[1],
             B=length(t),mean_T=mean(t),sd_T=sd(t),
             tail_probability=mean(abs(t)>qnorm(.975)),
             coverage_from_T=mean(abs(t)<=qnorm(.975)),
             skewness=mean((t-mean(t))^3)/mean((t-mean(t))^2)^(3/2))
}))
row.names(diag) <- NULL
write.csv(diag,'results/figures/E3_diagnostics_summary.csv',row.names=FALSE)
cat('Figures creees dans results/figures/\n')
print(diag[,c('scenario','n','mean_depth','k','mean_T','sd_T','coverage_from_T')],row.names=FALSE)
