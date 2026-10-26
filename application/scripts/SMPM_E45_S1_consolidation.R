#!/usr/bin/env Rscript
# E4.5 S1 SMPM consolidation: grid resolution, signed regional stability,
# spectral alignment, mark omission, batch and depth sensitivity.
# Run: Rscript SMPM_E45_S1_consolidation.R INPUT_DIR OUTPUT_DIR [GENE_BED]
# Requires only base R; genome annotation is optional and MUST match mm10/GRCm38.
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>=1) args[1] else 'data'
out <- if(length(args)>=2) args[2] else 'validation/consolidation_rerun'
annotation <- if(length(args)>=3) args[3] else ''
dir.create(out,recursive=TRUE,showWarnings=FALSE)
marks <- c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3')
classes <- c('EPI','PE','TE'); Lchr <- 195471971
meta <- read.csv(file.path(input,'E45_S1_cell_depths.csv'),check.names=FALSE)
stopifnot(all(c('barcode','celltype','batch',paste0('D_',marks)) %in% names(meta)))
meta <- meta[meta$celltype %in% classes,,drop=FALSE]
stopifnot(!anyDuplicated(meta$barcode),nrow(meta)>20)
n <- nrow(meta); q <- length(marks); y <- factor(meta$celltype,levels=classes)
p <- as.numeric(table(y))/n
D <- as.matrix(meta[,paste0('D_',marks),drop=FALSE]); stopifnot(all(is.finite(D)),all(D>0))
# Read once, retaining chromosome-1 atom coordinates, cell index and genome-wide depth.
atom <- vector('list',q)
for(r in seq_along(marks)) {
  path <- file.path(input,sprintf('E45_S1_%s_empirical_atoms.tsv.gz',marks[r]))
  con <- gzfile(path,'rt'); hdr <- strsplit(readLines(con,n=1L),'\t',fixed=TRUE)[[1]]
  needed <- c('barcode','chrom','midpoint_bp','multiplicity','D_ir')
  stopifnot(all(needed %in% hdr))
  ix <- match(needed,hdr); chunks <- list(); j <- 0L
  repeat {
    ln <- readLines(con,n=50000L,warn=FALSE); if(!length(ln)) break
    z <- matrix(unlist(strsplit(ln,'\t',fixed=TRUE),use.names=FALSE),ncol=length(hdr),byrow=TRUE)[,ix,drop=FALSE]
    keep <- z[,2]=='chr1'; if(!any(keep)) next
    z <- z[keep,,drop=FALSE]; ii <- match(z[,1],meta$barcode)
    keep <- !is.na(ii); if(!any(keep)) next
    z <- z[keep,,drop=FALSE]; ii <- ii[keep]
    pos <- as.numeric(z[,3]); mult <- as.numeric(z[,4]); depth <- as.numeric(z[,5])
    stopifnot(all(is.finite(pos)),all(pos>=0 & pos<=Lchr),all(mult>0),all(depth>0))
    stopifnot(all(abs(depth-D[ii,r])<1e-8))
    j <- j+1L; chunks[[j]] <- data.frame(i=ii,pos=pos,mass=mult/depth)
  }
  close(con); atom[[r]] <- do.call(rbind,chunks)
  message(marks[r],': ',nrow(atom[[r]]),' chr1 atoms')
}
make_X <- function(M) {
  X <- matrix(0,n,q*(M+1L))
  for(r in seq_along(marks)) {
    a <- atom[[r]]; first <- pmin(M+1L,as.integer(ceiling(a$pos/Lchr*M-1e-12))+1L)
    lin <- a$i+(first-1L)*n
    agg <- rowsum(a$mass,lin,reorder=FALSE)
    dif <- matrix(0,n,M+1L); dif[as.integer(rownames(agg))] <- as.numeric(agg[,1])
    f <- t(apply(dif,1,cumsum))
    X[,(r-1L)*(M+1L)+seq_len(M+1L)] <- f
  }
  X
}
fit_means <- function(gm,p,w,mask=rep(TRUE,ncol(gm))) {
  mu <- colSums(gm*p); A <- sweep(gm,2,mu,'-')
  Z <- sweep(A[,mask,drop=FALSE],2,sqrt(w[mask]),'*')*sqrt(p)
  ee <- eigen(tcrossprod(Z),symmetric=TRUE)
  lam <- pmax(ee$values[1:2],0)
  phi <- matrix(0,ncol(gm),2)
  for(k in 1:2) if(lam[k]>1e-16) phi[mask,k] <- as.vector(crossprod(Z,ee$vectors[,k]))/sqrt(lam[k])/sqrt(w[mask])
  list(lambda=lam,phi=phi)
}
normalize_cols <- function(A) {
  Q <- qr.Q(qr(A)); if(ncol(Q)<2L) stop('Rank below two')
  Q[,1:2,drop=FALSE]
}
subcos <- function(ref,alt,w,mask=rep(TRUE,length(w))) {
  A <- normalize_cols(sweep(ref$phi[mask,,drop=FALSE],1,sqrt(w[mask]),'*'))
  B <- normalize_cols(sweep(alt$phi[mask,,drop=FALSE],1,sqrt(w[mask]),'*'))
  svd(crossprod(A,B),nu=0,nv=0)$d
}
# Align each eigenvector to its best corresponding reference direction using
# weighted absolute correlations and a one-to-one permutation; no rotations
# are imposed when reporting signed individual measures.
align_phi <- function(ref,alt,w,mask=rep(TRUE,length(w))) {
  corr <- crossprod(ref$phi[mask,,drop=FALSE]*w[mask],alt$phi[mask,,drop=FALSE])
  denom <- outer(sqrt(colSums(ref$phi[mask,,drop=FALSE]^2*w[mask])),
                 sqrt(colSums(alt$phi[mask,,drop=FALSE]^2*w[mask])))
  corr <- corr/denom
  perm <- if(sum(abs(diag(corr))) >= sum(abs(corr[cbind(1:2,2:1)]))) 1:2 else 2:1
  signs <- sign(corr[cbind(1:2,perm)]); signs[signs==0] <- 1
  z <- alt; z$phi <- sweep(alt$phi[,perm,drop=FALSE],2,signs,'*')
  z$lambda <- alt$lambda[perm]
  list(fit=z,perm=perm,sign=signs,correlation=abs(corr[cbind(1:2,perm)]))
}
# Region signed mass = difference in cumulative function at boundaries.
# This avoids grid-boundary attribution artifacts from allocating jumps to nodes.
region_table <- function(fit,M,scenario,edges=unique(c(seq(0,Lchr,by=1e7),Lchr))) {
  outrows <- list(); idx <- 0L
  for(k in 1:2) for(r in seq_along(marks)) {
    v <- fit$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
    # Piecewise-constant right-continuous grid representation.
    # A jump assigned to node j belongs to (node[j-1],node[j]].
    # At boundaries evaluate cumulative value at grid node <= boundary.
    node <- pmin(M,as.integer(floor(edges/Lchr*M+1e-12)))+1L
    cum <- v[node]
    mass <- diff(cum)
    idx <- idx+1L
    outrows[[idx]] <- data.frame(scenario=scenario,grid_M=M,component=k,mark=marks[r],
      start_bp=head(edges,-1),end_bp=tail(edges,-1),signed_mass=mass,
      abs_mass=abs(mass))
  }
  do.call(rbind,outrows)
}
rows <- list(); regs <- list(); fits <- list(); xbase <- NULL
for(M in c(500L,1000L,2000L)) {
  message('Fitting grid M=',M)
  X <- make_X(M)
  w <- rep(c(.5,rep(1,M-1L),.5)/M,q)
  gm <- do.call(rbind,lapply(classes,function(g) colMeans(X[y==g,,drop=FALSE])))
  full <- fit_means(gm,p,w)
  if(M==1000L) {xbase <- X; wbase <- w; fullbase <- full; gmbase <- gm}
  fits[[as.character(M)]] <- list(M=M,fit=full,w=w)
  rows[[length(rows)+1L]] <- data.frame(scenario='original',grid_M=M,
    lambda1=full$lambda[1],lambda2=full$lambda[2],cosine_min=1,cosine_max=1,
    component1_abs_corr=1,component2_abs_corr=1,swapped=FALSE)
  regs[[length(regs)+1L]] <- region_table(full,M,'original')
  if(M!=1000L) rm(X)
}
# Compare grid 500 and 2000 to reference 1000 on a COMMON 2000-interval grid,
# interpolating cumulative functions (piecewise-linear) and using same L2 weights.
interpolate_phi <- function(fit,M,target=2000L) {
  z <- matrix(0,q*(target+1L),2)
  for(r in seq_along(marks)) for(k in 1:2) {
    old <- fit$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
    z[(r-1L)*(target+1L)+seq_len(target+1L),k] <-
      approx(seq(0,1,length.out=M+1L),old,xout=seq(0,1,length.out=target+1L),ties='ordered')$y
  }
  list(phi=z,lambda=fit$lambda)
}
wcommon <- rep(c(.5,rep(1,1999L),.5)/2000,q)
refgrid <- interpolate_phi(fits[['1000']]$fit,1000L)
for(M in c(500L,2000L)) {
  alt <- interpolate_phi(fits[[as.character(M)]]$fit,M)
  al <- align_phi(refgrid,alt,wcommon); cc <- subcos(refgrid,alt,wcommon)
  rows[[length(rows)+1L]] <- data.frame(scenario='grid_vs_1000',grid_M=M,
    lambda1=al$fit$lambda[1],lambda2=al$fit$lambda[2],
    cosine_min=min(cc),cosine_max=max(cc),component1_abs_corr=al$correlation[1],
    component2_abs_corr=al$correlation[2],swapped=!identical(al$perm,1:2))
}
# At M=1000 assess removal of low-depth H3K27me3, batch-standardization,
# and adjustment for centered log1p(depth) + batch.
M <- 1000L; X <- xbase; w <- wbase; full <- fullbase; gm <- gmbase
mask <- rep(marks!='H3K27me3',each=M+1L)
reduced <- fit_means(gm,p,w,mask)
batches <- sort(unique(as.character(meta$batch)))
stopifnot(length(batches)==2L)
for(g in classes) for(b in batches) stopifnot(sum(y==g & meta$batch==b)>=5L)
pb <- as.numeric(table(factor(meta$batch,levels=batches)))/n
gm_batch <- do.call(rbind,lapply(classes,function(g) {
  mm <- do.call(rbind,lapply(batches,function(b) colMeans(X[y==g & meta$batch==b,,drop=FALSE])))
  colSums(mm*pb)
}))
batchfit <- fit_means(gm_batch,p,w)
LD <- log1p(D); LD <- sweep(LD,2,colMeans(LD),'-')
Z <- cbind(
  sapply(classes, function(g) as.numeric(y == g)),
  LD,
  as.numeric(meta$batch == batches[2]) - pb[2]
)
stopifnot(qr(Z)$rank==ncol(Z))
coef <- qr.coef(qr(Z),X)
adjfit <- fit_means(coef[seq_along(classes),,drop=FALSE],p,w)
scenarios <- list(omit_H3K27me3=reduced,common_batch_composition=batchfit,
                  adjust_batch_and_logdepth=adjfit)
for(nm in names(scenarios)) {
  z <- scenarios[[nm]]; use <- if(nm=='omit_H3K27me3') mask else rep(TRUE,length(w))
  al <- align_phi(full,z,w,use); cc <- subcos(full,z,w,use)
  rows[[length(rows)+1L]] <- data.frame(scenario=nm,grid_M=M,
    lambda1=al$fit$lambda[1],lambda2=al$fit$lambda[2],
    cosine_min=min(cc),cosine_max=max(cc),component1_abs_corr=al$correlation[1],
    component2_abs_corr=al$correlation[2],swapped=!identical(al$perm,1:2))
  regs[[length(regs)+1L]] <- region_table(al$fit,M,nm)
}
summary <- do.call(rbind,rows)
write.csv(summary,file.path(out,'spectral_and_grid_sensitivity.csv'),row.names=FALSE)
reg <- do.call(rbind,regs)
write.csv(reg,file.path(out,'aligned_regional_signed_masses.csv'),row.names=FALSE)
# Compare regional vectors for each component, across all marks. Signed Pearson
# correlation and relative L2 distance are complementary (scale-sensitive vs invariant).
base <- reg[reg$scenario=='original' & reg$grid_M==1000L,,drop=FALSE]
comparisons <- list()
for(nm in setdiff(unique(reg$scenario),'original')) {
  alt <- reg[reg$scenario==nm & reg$grid_M==1000L,,drop=FALSE]
  for(k in 1:2) {
    b <- base[base$component==k,]; a <- alt[alt$component==k,]
    key <- function(z) paste(z$mark,z$start_bp,z$end_bp,sep='|')
    a <- a[match(key(b),key(a)),]; ok <- !is.na(a$signed_mass)
    if(nm=='omit_H3K27me3') ok <- ok & b$mark!='H3K27me3'
    v <- b$signed_mass[ok]; u <- a$signed_mass[ok]
    comparisons[[length(comparisons)+1L]] <- data.frame(scenario=nm,component=k,
      n_region_mark=length(v),signed_correlation=if(sd(v)>0 && sd(u)>0) cor(v,u) else NA_real_,
      relative_L2_distance=sqrt(sum((v-u)^2))/sqrt(sum(v^2)),
      top10_overlap={top <- function(z) order(abs(z),decreasing=TRUE)[seq_len(min(10,length(z)))];
        length(intersect(top(v),top(u)))/min(10,length(v))})
  }
}
# Grid regional comparisons (same fixed 10 Mb windows).
for(M in c(500L,2000L)) for(k in 1:2) {
  b <- base[base$component==k,]
  a <- reg[reg$scenario=='original' & reg$grid_M==M & reg$component==k,]
  key <- function(z) paste(z$mark,z$start_bp,z$end_bp,sep='|')
  a <- a[match(key(b),key(a)),]
  # Sign-align using interpolated functional eigenvectors.
  alt <- interpolate_phi(fits[[as.character(M)]]$fit,M)
  al <- align_phi(refgrid,alt,wcommon)
  sign <- al$sign[k]
  u <- a$signed_mass*sign; v <- b$signed_mass
  comparisons[[length(comparisons)+1L]] <- data.frame(scenario=paste0('grid_',M),component=k,
    n_region_mark=length(v),signed_correlation=cor(v,u),
    relative_L2_distance=sqrt(sum((v-u)^2))/sqrt(sum(v^2)),
    top10_overlap=length(intersect(order(abs(v),decreasing=TRUE)[1:10],
                                   order(abs(u),decreasing=TRUE)[1:10]))/10)
}
comp <- do.call(rbind,comparisons)
write.csv(comp,file.path(out,'regional_stability_metrics.csv'),row.names=FALSE)
# Rank 10-Mb regions by combined absolute contribution across marks, not a formal
# peak discovery procedure. Signed masses have already been direction-aligned.
region_rank <- aggregate(abs_mass~scenario+grid_M+component+start_bp+end_bp,reg,sum)
region_rank <- region_rank[order(region_rank$scenario,region_rank$grid_M,
  region_rank$component,-region_rank$abs_mass),]
write.csv(region_rank,file.path(out,'regional_rankings.csv'),row.names=FALSE)
# Optional independent gene annotation: BED 0-based half-open coordinates,
# columns chrom/start/end/gene. No online retrieval or genome-build guessing.
if(nzchar(annotation)) {
  stopifnot(file.exists(annotation))
  bed <- read.table(annotation,sep='\t',header=FALSE,quote='',comment.char='#',
                    stringsAsFactors=FALSE,fill=TRUE)
  stopifnot(ncol(bed)>=4L)
  bed <- bed[bed[[1]] %in% c('chr1','1'),1:4,drop=FALSE]
  names(bed) <- c('chrom','start0','end0','gene')
  stopifnot(all(is.finite(bed$start0)),all(is.finite(bed$end0)),
            all(bed$start0>=0),all(bed$end0>bed$start0),all(bed$end0<=Lchr))
  # A gene interval overlaps 10 Mb region when gene end > region start and
  # gene start < region end. These are overlaps, NOT causal gene assignments.
  rr <- subset(region_rank,scenario=='original' & grid_M==1000L)
  rr <- rr[order(rr$component,-rr$abs_mass),]
  rr <- do.call(rbind,lapply(split(rr,rr$component),head,10))
  rr$overlapping_genes <- vapply(seq_len(nrow(rr)),function(i)
    paste(unique(bed$gene[bed$end0>rr$start_bp[i] & bed$start0<rr$end_bp[i]]),collapse=';'),character(1))
  write.csv(rr,file.path(out,'candidate_regions_gene_overlaps.csv'),row.names=FALSE)
}
writeLines(c('SMPM CONSOLIDATION NOTES',
 'All numeric comparisons are descriptive and not independent biological validation.',
 'Grid sensitivity uses common continuous interpolation for functional subspace comparison.',
 'Signed regional masses use cumulative differences at fixed 10 Mb boundaries.',
 'Individual eigenvectors are aligned by permutation and sign only, not Procrustes rotation.',
 'When eigenvalues nearly coincide, compare subspaces rather than individual components.',
 'IMPORTANT: for exactly 3 phenotype groups and rank 2, removing one mark',
 'can produce canonical cosine 1 identically on the shared coordinates.',
 'Therefore the 1.000 result is algebraically expected and NOT evidence of robustness.',
 'Depth regression can remove biological signal; it is not a test of unbiased capture.',
 'The optional BED annotation must be independently sourced and verified mm10/GRCm38.',
 'A BED overlap is not proof of gene regulation or biological mechanism.',
 paste('R:',R.version.string)),file.path(out,'README.txt'))
print(summary); print(comp)
message('Completed: ',normalizePath(out))
