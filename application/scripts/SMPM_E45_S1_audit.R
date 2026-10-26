#!/usr/bin/env Rscript
# Audit of the SMPM E4.5 S1 analysis. Run in the folder containing the input files.
# Rscript SMPM_E45_S1_audit.R /path/to/input /path/to/audit_output
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>=1) args[1] else 'data'
out <- if(length(args)>=2) args[2] else 'validation/audit_rerun'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
original <- file.path('scripts','SMPM_E45_S1_application.R')
if(!file.exists(original)) stop('Original script missing: ',original)
# Run the original estimation in its own environment, using the existing input and results.
# The original script reads commandArgs; isolate its outputs to a separate subfolder
# by overriding the script's out assignment in memory, without editing the original.
code <- readLines(original,warn=FALSE)
line <- grep("^out <- if\\(length\\(args\\)>=2\\)",code)
stopifnot(length(line)==1L)
code[line] <- sprintf('out <- %s',deparse(file.path(out,'original_rerun')))
env <- new.env(parent=globalenv())
env$input <- input
# commandArgs still has original values; inject input assignment as well
li <- grep("^input <- if\\(length\\(args\\)>=1\\)",code)
stopifnot(length(li)==1L)
code[li] <- sprintf('input <- %s',deparse(input))
eval(parse(text=code),envir=env)
X <- env$X; meta <- env$meta; marks <- env$marks; classes <- env$classes
M <- env$M; Lchr <- env$Lchr; w <- env$w; sw <- env$sw; full <- env$full
n <- nrow(X); q <- length(marks); p <- as.numeric(table(factor(meta$celltype,levels=classes)))/n
# Correct grid: each increment is assigned to (t_{j-1}, t_j] (first grid node 0 is a point mass).
inc <- list()
for(k in 1:2) for(r in seq_along(marks)) {
  v <- full$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
  left <- c(0,(0:(M-1L))/M*Lchr)
  right <- (0:M)/M*Lchr
  inc[[length(inc)+1L]] <- data.frame(component=k,mark=marks[r],
     left_bp=round(left),right_bp=round(right),signed_mass=c(v[1],diff(v)),
     interpretation=c('point_mass_at_zero',rep('grid_increment_interval',M)))
}
inc <- do.call(rbind,inc)
write.csv(inc,file.path(out,'corrected_grid_increments.csv'),row.names=FALSE)
# The first row (at zero) is a point mass; others must have positive interval widths.
stopifnot(all(inc$right_bp[inc$interpretation=='grid_increment_interval'] >
              inc$left_bp[inc$interpretation=='grid_increment_interval']))

# The spectral estimator based on arbitrary group means; this is the same weighted Gram approach.
spectral_means <- function(gm,p,mask=rep(TRUE,ncol(gm))) {
  overall <- colSums(gm*p)
  A <- sweep(gm,2,overall,'-')
  Z <- sweep(A[,mask,drop=FALSE],2,sw[mask],'*')*sqrt(p)
  e <- eigen(tcrossprod(Z),symmetric=TRUE)
  val <- pmax(e$values[1:2],0)
  phi <- matrix(0,ncol(X),2)
  for(k in 1:2) if(val[k]>1e-16) {
    phi[mask,k] <- as.vector(crossprod(Z,e$vectors[,k]))/sqrt(val[k])/sw[mask]
  }
  list(lambda=val,phi=phi)
}
alignment <- function(a,b,mask=rep(TRUE,ncol(X))) {
  # Full and reduced subspaces are compared on their common mark coordinates.
  A <- sweep(a$phi[mask,,drop=FALSE],1,sw[mask],'*')
  B <- sweep(b$phi[mask,,drop=FALSE],1,sw[mask],'*')
  QA <- qr.Q(qr(A)); QB <- qr.Q(qr(B))
  svd(crossprod(QA,QB),nu=0,nv=0)$d
}
base_means <- full$gm
no27 <- !rep(marks=='H3K27me3',each=M+1L)
reduced <- spectral_means(base_means,p,no27)
# Compare like with like: re-estimate full analysis restricted to the 4 common marks.
base_common <- spectral_means(base_means,p,no27)
# Comparison of the 5-mark supervised subspace projected to 4 marks with 4-mark refit.
cc27 <- alignment(full,reduced,no27)
# Relative contribution of each mark to squared L2 norm of each eigenfunction.
mark_norm <- do.call(rbind,lapply(1:2,function(k) data.frame(
 component=k,mark=marks,L2_energy=vapply(seq_along(marks),function(r) {
   ii <- (r-1L)*(M+1L)+seq_len(M+1L)
   sum(w[ii]*full$phi[ii,k]^2)
 },numeric(1)))))
write.csv(mark_norm,file.path(out,'mark_energy.csv'),row.names=FALSE)

# Batch-adjusted mean curves, at a common batch composition (global batch frequencies).
batches <- sort(unique(as.character(meta$batch)))
if(length(batches)!=2L) stop('Expected two batches; found: ',paste(batches,collapse=', '))
batch_means <- lapply(classes,function(g) do.call(rbind,lapply(batches,function(b) {
  idx <- which(meta$celltype==g & meta$batch==b)
  if(length(idx)<5L) stop('Insufficient cells for ',g,' in ',b)
  colMeans(X[idx,,drop=FALSE])
})))
common_batch <- as.numeric(table(factor(meta$batch,levels=batches)))/n
gm_batch <- do.call(rbind,lapply(batch_means,function(z) colSums(z*common_batch)))
fit_batch <- spectral_means(gm_batch,p)
cc_batch <- alignment(full,fit_batch)

# Joint adjusted means from multivariate linear regression of profiles on group,
# centered log depth and batch. This is a sensitivity model, NOT proof of unbiasedness.
D <- as.matrix(meta[,paste0('D_',marks),drop=FALSE])
LD <- log1p(D)
LD <- sweep(LD,2,colMeans(LD),'-')
Z <- cbind(sapply(classes,function(g) as.numeric(meta$celltype==g)),
          LD,as.numeric(meta$batch==batches[2])-common_batch[2])
colnames(Z) <- c(paste0('group_',classes),paste0('logdepth_',marks),'batch')
if(qr(Z)$rank!=ncol(Z)) stop('Depth/batch regression design is rank deficient')
coef <- qr.coef(qr(Z),X)
gm_adj <- coef[seq_along(classes),,drop=FALSE]
fit_adj <- spectral_means(gm_adj,p)
cc_adj <- alignment(full,fit_adj)

summary <- rbind(
 data.frame(scenario='original_5_marks',lambda1=full$lambda[1],lambda2=full$lambda[2],cosine_min=1,cosine_max=1),
 data.frame(scenario='omit_H3K27me3_4_marks',lambda1=reduced$lambda[1],lambda2=reduced$lambda[2],cosine_min=min(cc27),cosine_max=max(cc27)),
 data.frame(scenario='common_batch_composition',lambda1=fit_batch$lambda[1],lambda2=fit_batch$lambda[2],cosine_min=min(cc_batch),cosine_max=max(cc_batch)),
 data.frame(scenario='adjust_batch_and_logdepth',lambda1=fit_adj$lambda[1],lambda2=fit_adj$lambda[2],cosine_min=min(cc_adj),cosine_max=max(cc_adj)))
write.csv(summary,file.path(out,'audit_spectral_sensitivity.csv'),row.names=FALSE)
write.csv(as.data.frame.matrix(table(meta$celltype,meta$batch)),file.path(out,'lineage_by_batch.csv'))
write.csv(data.frame(celltype=classes,do.call(rbind,lapply(classes,function(g)
 apply(D[meta$celltype==g,,drop=FALSE],2,median)))),file.path(out,'depth_medians_by_lineage.csv'),row.names=FALSE)

# Region-level changes in absolute signed contributions. Do not claim biological
# annotation: gene tracks, mappability and blacklists are not included in these inputs.
region_edges <- unique(c(seq(0,Lchr,by=1e7),Lchr))
region_table <- function(fit,scenario) {
  z <- list()
  for(k in 1:2) for(r in seq_along(marks)) {
    phi <- fit$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
    coord <- (0:M)/M*Lchr
    idx <- pmin(length(region_edges)-1L,findInterval(coord,region_edges,rightmost.closed=TRUE))
    totals <- numeric(length(region_edges)-1L)
    agg <- rowsum(c(phi[1],diff(phi)),idx,reorder=FALSE)
    totals[as.integer(rownames(agg))] <- agg[,1]
    z[[length(z)+1L]] <- data.frame(scenario=scenario,component=k,mark=marks[r],
      start_bp=head(region_edges,-1),end_bp=tail(region_edges,-1),signed_mass=totals)
  }
  do.call(rbind,z)
}
reg <- do.call(rbind,list(region_table(full,'original'),region_table(reduced,'omit_H3K27me3'),
  region_table(fit_batch,'batch_standardized'),region_table(fit_adj,'batch_depth_adjusted')))
write.csv(reg,file.path(out,'audit_regional_contributions.csv'),row.names=FALSE)
writeLines(c('AUDIT NOTES',
 'Corrected increments: first row represents a point at 0; all other rows represent grid increments on (left,right].',
 'SMPM directions are identified up to sign; signed regional masses across scenarios require sign alignment.',
 'Principal cosines compare 2D subspaces, not individual signed principal measures.',
 'Removing H3K27me3 changes the Hilbert space; alignment is computed on the 4 shared mark coordinates.',
 'Depth adjustment is an observational regression sensitivity, not validation of conditional unbiasedness.',
 'Batch standardization uses global batch proportions shared across lineages.',
 'No biological gene annotation or mappability/blacklist verification is performed.'),file.path(out,'README_audit.txt'))
print(summary)
cat('Audit written to ',normalizePath(out),'\n')
