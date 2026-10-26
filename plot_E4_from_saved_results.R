# E4 figure only: reuse saved Monte Carlo results; no simulations.
graphics.off()
e4 <- readRDS('results/tables/E4_raw_results.rds')
summary_e4 <- e4$summary
rhos <- sort(unique(summary_e4$rho))
stopifnot(length(rhos) == 4L)
getmean <- function(metric_name) {
  z <- summary_e4[summary_e4$metric == metric_name, , drop = FALSE]
  z <- z[order(z$rho), , drop = FALSE]
  if (nrow(z) != length(rhos) || !isTRUE(all.equal(z$rho, rhos))) {
    stop('Missing or duplicated metric/rho combinations: ', metric_name)
  }
  z$mean
}
sp1 <- getmean('S_phi1'); un1 <- getmean('U_phi1')
sp2 <- getmean('S_subspace_min')
un2 <- getmean('U_top2_subspace_min')
un3 <- getmean('U_top3_subspace_min')
rho1_star <- 0.4 * sqrt(2/3)
rho2_star <- 0.3 / sqrt(2)
blue <- '#176B87'; orange <- '#C05B39'; gray <- '#666666'
dir.create('results/figures', recursive = TRUE, showWarnings = FALSE)
outfile <- 'results/figures/Fig_E4_supervision_corrected.pdf'
pdf(outfile, width = 12.2, height = 5.4, family = 'Helvetica', useDingbats = FALSE)
par(mfrow = c(1, 2), mar = c(7.1, 4.7, 3.0, 0.9),
    mgp = c(2.6, 0.65, 0), las = 1, xpd = FALSE)
# Panel (a)
plot(rhos, sp1, type = 'b', pch = 16, lwd = 2, col = blue,
     ylim = c(0, 1.04), xlim = c(-0.025, 0.73), xaxt = 'n',
     xlab = expression('Nuisance amplitude ' * rho),
     ylab = 'Alignment with first phenotype direction')
axis(1, at = rhos)
lines(rhos, un1, type = 'b', pch = 17, lwd = 2, col = orange)
abline(v = rho1_star, lty = 3, col = 'gray50', lwd = 1.2)
title(main = '(a) Leading phenotype direction', adj = 0, font.main = 2, cex.main = 0.96)
par(xpd = NA)
legend(x = par('usr')[1], y = par('usr')[3] - 0.18,
       legend = c('SMPM', 'Unsupervised FPCA (first component)'),
       lty = 1, pch = c(16, 17), col = c(blue, orange),
       bty = 'n', cex = 0.77, xjust = 0, yjust = 1)
par(xpd = FALSE)
# Panel (b)
plot(rhos, sp2, type = 'b', pch = 16, lwd = 2, col = blue,
     ylim = c(0, 1.04), xlim = c(-0.025, 0.73), xaxt = 'n',
     xlab = expression('Nuisance amplitude ' * rho),
     ylab = 'Minimum canonical alignment')
axis(1, at = rhos)
lines(rhos, un2, type = 'b', pch = 17, lwd = 2, col = orange)
lines(rhos, un3, type = 'b', pch = 15, lwd = 2, lty = 2, col = gray)
abline(v = rho2_star, lty = 3, col = 'gray50', lwd = 1.2)
title(main = '(b) Phenotype subspace recovery', adj = 0, font.main = 2, cex.main = 0.96)
par(xpd = NA)
legend(x = par('usr')[1], y = par('usr')[3] - 0.18,
       legend = c('SMPM (2 components)', 'Unsupervised FPCA (2 components)',
                  'Unsupervised FPCA (3 components)'),
       lty = c(1, 1, 2), pch = c(16, 17, 15),
       col = c(blue, orange, gray), bty = 'n', cex = 0.73,
       xjust = 0, yjust = 1)
par(xpd = FALSE)
dev.off()
cat('Figure saved:', outfile, '\n')
