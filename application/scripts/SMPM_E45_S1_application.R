#!/usr/bin/env Rscript
# Reproducible E4.5 S1 supervised multivariate principal measures, chr1
# Input: E45_S1_cell_depths.csv, five E45_S1_*_empirical_atoms.tsv.gz,
# optional S1_multigene_scores_by_cell.csv. Base R only; irlba optional.
# Usage: Rscript SMPM_E45_S1_application.R /path/to/data /path/to/results
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>=1) args[1] else 'data'
out <- if(length(args)>=2) args[2] else 'results'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261011)
marks <- c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3')
classes <- c('EPI','PE','TE'); M <- 1000L; B <- 200L
# mm10 chromosome 1 length (bp); check the reference assembly of input fragments.
Lchr <- 195471971
meta <- read.csv(file.path(input,'E45_S1_cell_depths.csv'),check.names=FALSE)
stopifnot(all(c('barcode','celltype',paste0('D_',marks)) %in% names(meta)))
meta <- meta[meta$celltype %in% classes,,drop=FALSE]
stopifnot(!anyDuplicated(meta$barcode),all(is.finite(as.matrix(meta[,paste0('D_',marks)]))),
          all(as.matrix(meta[,paste0('D_',marks)])>0))
n <- nrow(meta); stopifnot(n>3L)
message('Cells: ',n,'; groups: ',paste(names(table(meta$celltype)),table(meta$celltype),collapse=', '))
# A right-continuous CDF on u_l=l/M, including l=0. Each atom contributes
# to all grid nodes with u_l >= midpoint/Lchr. Read in chunks, chromosome 1 only.
X <- matrix(0,nrow=n,ncol=length(marks)*(M+1L))
rownames(X) <- meta$barcode
names_by_mark <- vector('list',length(marks))
for(r in seq_along(marks)) {
  path <- file.path(input,sprintf('E45_S1_%s_empirical_atoms.tsv.gz',marks[r]))
  if(!file.exists(path)) stop('Missing atom file: ',path)
  con <- gzfile(path,'rt'); hdr <- strsplit(readLines(con,n=1L),'\t',fixed=TRUE)[[1]]
  needed <- c('barcode','chrom','midpoint_bp','multiplicity','D_ir')
  stopifnot(all(needed %in% hdr))
  ix <- match(needed,hdr); dif <- matrix(0,nrow=n,ncol=M+2L)
  count <- 0L
  repeat {
    ln <- readLines(con,n=40000L,warn=FALSE); if(!length(ln)) break
    spl <- strsplit(ln,'\t',fixed=TRUE)
    if(any(lengths(spl)!=length(hdr))) stop('Malformed atom file: ',path)
    z <- matrix(unlist(spl,use.names=FALSE),ncol=length(hdr),byrow=TRUE)[,ix,drop=FALSE]
    keep <- z[,2]=='chr1'; if(!any(keep)) next
    z <- z[keep,,drop=FALSE]; i <- match(z[,1],meta$barcode)
    ok <- !is.na(i); if(!any(ok)) next
    z <- z[ok,,drop=FALSE]; i <- i[ok]
    pos <- as.numeric(z[,3]); mult <- as.numeric(z[,4]); d <- as.numeric(z[,5]);
    stopifnot(all(is.finite(pos)),all(pos>=0 & pos<=Lchr),all(d>0),all(mult>0))
    # Verify atom depth matches independently provided per-cell genome-wide depth
    if(any(abs(d-meta[[paste0('D_',marks[r])]][i])>1e-8)) stop('Depth mismatch for ',marks[r])
    # first grid index 1..M+1 with u >= pos/Lchr (grid node l=0,...,M)
    first <- pmin(M+1L,as.integer(ceiling(pos/Lchr*M-1e-12))+1L)
    # accumulate signed first differences via grouped tabulation
    lin <- i+(first-1L)*n
    accum <- rowsum(mult/d,lin,reorder=FALSE)
    ids <- as.integer(rownames(accum))
    dif[ids] <- dif[ids]+as.numeric(accum[,1]); count <- count+length(i)
  }
  close(con)
  f <- t(apply(dif[,seq_len(M+1L),drop=FALSE],1,cumsum))
  stopifnot(all(f>=-1e-10),all(f<=1+1e-8))
  X[,(r-1L)*(M+1L)+seq_len(M+1L)] <- f
  names_by_mark[[r]] <- paste0(marks[r],':',0:M)
  message(marks[r],': ',count,' chr1 fragment records; max chr1 mass ',round(max(f[,M+1L]),4))
}
colnames(X) <- unlist(names_by_mark)
write.csv(data.frame(barcode=meta$barcode,celltype=meta$celltype,
  sapply(seq_along(marks),function(r) X[,r*(M+1L)])),
  file.path(out,'chr1_terminal_masses.csv'),row.names=FALSE)
# trapezoidal quadrature on [0,1], independently for each histone mark
w <- rep(c(0.5,rep(1,M-1L),0.5)/M,length(marks))
sw <- sqrt(w)
# spectral estimator via 3 x 3 Gram matrix, NOT a 5005 x 5005 matrix
fit <- function(ind) {
  yy <- factor(meta$celltype[ind],levels=classes)
  if(any(table(yy)<2L)) stop('Too few cells in a lineage')
  p <- as.numeric(table(yy))/length(ind)
  gm <- do.call(rbind,lapply(classes,function(g) colMeans(X[ind[yy==g],,drop=FALSE])))
  overall <- colSums(gm*p)
  A <- sweep(gm,2,overall,'-')
  Z <- sweep(A,2,sw,'*')*sqrt(p)
  ee <- eigen(tcrossprod(Z),symmetric=TRUE)
  vals <- pmax(ee$values[1:2],0)
  phi <- matrix(0,nrow=ncol(X),ncol=2L)
  for(k in 1:2) if(vals[k]>1e-16) phi[,k] <- as.vector(crossprod(Z,ee$vectors[,k]))/sqrt(vals[k])/sw
  # Deterministic sign: highest absolute coordinate is positive
  for(k in 1:2) if(phi[which.max(abs(phi[,k])),k]<0) phi[,k] <- -phi[,k]
  list(lambda=vals,phi=phi,p=p,gm=gm,mean=overall,ind=ind)
}
full <- fit(seq_len(n))
write.csv(data.frame(component=1:2,eigenvalue=full$lambda,
 explained_supervised=full$lambda/sum(full$lambda)),
 file.path(out,'supervised_eigenvalues.csv'),row.names=FALSE)
# Export the supervised eigenfunctions for exact figure reproduction.
# The ordering and sign convention are those of full$phi above.
eigenfunctions <- do.call(rbind,lapply(seq_along(marks),function(r) {
  idx <- (r-1L)*(M+1L)+seq_len(M+1L)
  data.frame(mark=marks[r],position_bp=(0:M)/M*Lchr,
             SMPC1=full$phi[idx,1],SMPC2=full$phi[idx,2])
}))
write.csv(eigenfunctions,file.path(out,'supervised_eigenfunctions.csv'),row.names=FALSE)

# Pointwise principal measure on the discretized grid: signed increments of phi.
# This is a discretization of the signed-measure representation, not raw genomic atoms.
windows <- list()
for(k in 1:2) for(r in seq_along(marks)) {
  phi <- full$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
  inc <- c(phi[1],diff(phi))
  # phi[1] represents a point mass at zero; subsequent differences
  # correspond to nondegenerate intervals (t_{j-1}, t_j].
  left <- c(0,(0:(M-1L))/M*Lchr)
  right <- (0:M)/M*Lchr
  windows[[length(windows)+1L]] <- data.frame(component=k,mark=marks[r],
    left_bp=round(left),right_bp=round(right),signed_mass=inc,
    interpretation=c('point_mass_at_zero',rep('grid_increment_interval',M)))
}
grid_increments <- do.call(rbind,windows)
stopifnot(nrow(grid_increments)==2L*length(marks)*(M+1L),
  all(grid_increments$right_bp[grid_increments$interpretation=='grid_increment_interval'] >
      grid_increments$left_bp[grid_increments$interpretation=='grid_increment_interval']),
  all(grid_increments$left_bp[grid_increments$interpretation=='point_mass_at_zero']==0),
  all(grid_increments$right_bp[grid_increments$interpretation=='point_mass_at_zero']==0))
write.csv(grid_increments,file.path(out,'principal_measure_grid_increments.csv'),row.names=FALSE)
# 10 Mb regional signed contributions: these are integrated signed increments,
# not inferential peak calls or tests of individual loci.
region_edges <- unique(c(seq(0,Lchr,by=10000000),Lchr))
regional <- list()
for(k in 1:2) for(r in seq_along(marks)) {
  phi <- full$phi[(r-1L)*(M+1L)+seq_len(M+1L),k]
  inc <- c(phi[1],diff(phi)); coord <- (0:M)/M*Lchr
  bin <- pmin(length(region_edges)-1L,findInterval(coord,region_edges,rightmost.closed=TRUE))
  v <- as.numeric(rowsum(inc,bin,reorder=FALSE)[,1]); bins <- as.integer(rownames(rowsum(inc,bin,reorder=FALSE)))
  regional[[length(regional)+1L]] <- data.frame(component=k,mark=marks[r],
    start_bp=region_edges[bins],end_bp=region_edges[bins+1L],signed_mass=v)
}
reg <- do.call(rbind,regional); reg$absolute_mass <- abs(reg$signed_mass)
write.csv(reg,file.path(out,'regional_signed_contributions_10Mb.csv'),row.names=FALSE)
# Scores in the supervised Hilbert directions (scores for display, not validation)
xc <- sweep(X,2,full$mean,'-')
scores <- (xc*rep(w,each=n))%*%full$phi
write.csv(data.frame(barcode=meta$barcode,celltype=meta$celltype,SMPC1=scores[,1],SMPC2=scores[,2]),
 file.path(out,'cell_scores.csv'),row.names=FALSE)
pdf(file.path(out,'SMPM_application_figures.pdf'),width=11,height=7)
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(scores,col=as.integer(factor(meta$celltype,levels=classes)),pch=16,cex=.45,
 xlab='SMPC1 score',ylab='SMPC2 score',main='Phenotype-directed cell scores')
legend('topright',legend=classes,col=1:3,pch=16,cex=.8)
barplot(full$lambda,names.arg=c('SMPC1','SMPC2'),main='Supervised eigenvalues',ylab='Eigenvalue')
for(k in 1:2) {
  plot(NA,xlim=c(0,Lchr/1e6),ylim=range(full$phi[,k]),xlab='chr1 position (Mb)',
       ylab='Eigenfunction',main=paste('SMPC',k,'signed-measure cumulative function'))
  for(r in seq_along(marks)) lines((0:M)/M*Lchr/1e6,
    full$phi[(r-1L)*(M+1L)+seq_len(M+1L),k],col=r,lwd=1.5)
  legend('topright',legend=marks,col=seq_along(marks),lty=1,cex=.65)
}
dev.off()
# Stratified bootstrap: uncertainty of eigenvalues and principal directions.
boot <- matrix(NA_real_,B,2L)
for(b in seq_len(B)) {
  ind <- unlist(lapply(classes,function(g) sample(which(meta$celltype==g),sum(meta$celltype==g),replace=TRUE)))
  boot[b,] <- fit(ind)$lambda
}
write.csv(data.frame(component=1:2,eigenvalue=full$lambda,
  ci_low=apply(boot,2,quantile,.025),ci_high=apply(boot,2,quantile,.975)),
  file.path(out,'bootstrap_eigenvalues.csv'),row.names=FALSE)
# Compare subspaces using principal cosines, insensitive to sign and rotations.
subspace_alignment <- function(a,b) {
  A <- sweep(a$phi,1,sw,'*'); C <- sweep(b$phi,1,sw,'*')
  svd(crossprod(A,C),nu=0,nv=0)$d
}
# Sensitivity: discordant cells removed, then uncertain cells removed.
score_path <- file.path(input,'S1_multigene_scores_by_cell.csv')
sens <- data.frame(analysis='all',n=n,lambda1=full$lambda[1],lambda2=full$lambda[2],
   min_canonical_cosine=1,median_canonical_cosine=1)
if(file.exists(score_path)) {
  m <- read.csv(score_path,stringsAsFactors=FALSE)
  stopifnot(all(c('cellid','group','discordant','uncertain') %in% names(m)))
  mi <- match(meta$barcode,m$cellid)
  if(anyNA(mi)) stop('Missing marker scores for some S1 cells')
  if(any(as.character(m$group[mi])!=as.character(meta$celltype))) stop('Marker group mismatch')
  scenarios <- list(no_discordance=!as.logical(m$discordant[mi]),
                    no_uncertainty=!as.logical(m$uncertain[mi]),
                    no_discordance_or_uncertainty=!(as.logical(m$discordant[mi])|as.logical(m$uncertain[mi])))
  for(nm in names(scenarios)) {
    ind <- which(scenarios[[nm]])
    if(any(table(factor(meta$celltype[ind],levels=classes))<2L)) next
    z <- fit(ind); cc <- subspace_alignment(full,z)
    sens <- rbind(sens,data.frame(analysis=nm,n=length(ind),lambda1=z$lambda[1],lambda2=z$lambda[2],
      min_canonical_cosine=min(cc),median_canonical_cosine=median(cc)))
  }
}
# Depth sensitivity: retain cells above mark-specific lower quartiles in ALL 5 marks
D <- as.matrix(meta[,paste0('D_',marks),drop=FALSE])
cut <- apply(D,2,quantile,.25)
ind <- which(rowSums(sweep(D,2,cut,'>='))==length(marks))
if(all(table(factor(meta$celltype[ind],levels=classes))>=2L)) {
  z <- fit(ind); cc <- subspace_alignment(full,z)
  sens <- rbind(sens,data.frame(analysis='all_depths_above_Q1',n=length(ind),
     lambda1=z$lambda[1],lambda2=z$lambda[2],min_canonical_cosine=min(cc),median_canonical_cosine=median(cc)))
}
write.csv(sens,file.path(out,'sensitivity_results.csv'),row.names=FALSE)
# Optional unsupervised comparison (leading two functional PCs); never use PCA
# score separation alone to establish better biological relevance.
if(requireNamespace('irlba',quietly=TRUE)) {
  Xw <- sweep(xc,2,sw,'*')
  up <- irlba::irlba(Xw,nv=2,nu=2,center=NULL)
  U <- sweep(full$phi,1,sw,'*')
  cs <- svd(crossprod(U,up$v),nu=0,nv=0)$d
  write.csv(data.frame(canonical_cosine=cs),file.path(out,'supervised_unsupervised_alignment.csv'),row.names=FALSE)
} else message('Optional comparison skipped: install.packages("irlba")')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
message('Completed: ',normalizePath(out))
