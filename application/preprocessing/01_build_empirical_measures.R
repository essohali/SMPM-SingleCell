#!/usr/bin/env Rscript
# Reproducible scMTR-seq E4.5 S1 empirical point measures (base R only).
# Usage: Rscript 01_build_empirical_measures.R INPUT_DIR ANNOTATIONS_CSV OUTPUT_DIR
# Annotation source: Wang et al. supplementary Table S2, E4.5 sheet,
# supplied as S2_E45_raw.csv (columns new.cellid, sub.lib, batch, celltype).
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) stop('Usage: Rscript 01_build_empirical_measures.R INPUT_DIR ANNOTATIONS_CSV OUTPUT_DIR')
in_dir <- args[1]; ann_file <- args[2]; out_dir <- args[3]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
marks <- c('H3K27ac','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K9me3')
ann <- read.csv(ann_file, stringsAsFactors=FALSE, check.names=FALSE)
need <- c('new.cellid','sub.lib','batch','celltype')
if (!all(need %in% names(ann))) stop('Missing annotation columns: ', paste(setdiff(need,names(ann)),collapse=', '))
ann <- ann[ann[['sub.lib']] == 'S1' & ann[['celltype']] %in% c('EPI','PE','TE'), need, drop=FALSE]
if (anyNA(ann) || anyDuplicated(ann[['new.cellid']])) stop('Missing or duplicate annotations')
if (nrow(ann) == 0L) stop('No S1 EPI/PE/TE annotations')
barcodes <- ann[['new.cellid']]; n <- length(barcodes)
celltype <- setNames(ann[['celltype']],barcodes); batch <- setNames(ann[['batch']],barcodes)
message('Annotated S1 cells: ',n,'; ',paste(names(table(ann$celltype)),as.integer(table(ann$celltype)),collapse=', '))
# Read gzip input in chunks; no extra R packages or whole-file load required.
chunk_size <- 50000L
parse_lines <- function(lines, filename) {
  fields <- strsplit(lines, '\t', fixed=TRUE)
  if (any(lengths(fields) != 5L)) stop('Invalid 5-column TSV row in ',filename)
  mat <- matrix(unlist(fields, use.names=FALSE), ncol=5L, byrow=TRUE)
  start <- suppressWarnings(as.numeric(mat[,2])); end <- suppressWarnings(as.numeric(mat[,3]))
  multiplicity <- suppressWarnings(as.numeric(mat[,5]))
  if (anyNA(start) || anyNA(end) || anyNA(multiplicity) ||
      any(end <= start) || any(start < 0) || any(multiplicity <= 0) ||
      any(multiplicity != floor(multiplicity))) stop('Invalid coordinates/multiplicities in ',filename)
  list(chrom=mat[,1], start=start, end=end, barcode=mat[,4], multiplicity=multiplicity)
}
input_path <- function(mark) file.path(in_dir,paste0('GSM8529560_E45_S1_',mark,'.fragments.tsv.gz'))
if (!all(file.exists(vapply(marks,input_path,character(1))))) stop('Missing files: ',paste(vapply(marks,input_path,character(1))[!file.exists(vapply(marks,input_path,character(1)))],collapse=', '))
depths <- matrix(0, nrow=n,ncol=length(marks),dimnames=list(barcodes,marks))
rowcounts <- depths
file_audit <- data.frame(mark=marks, total_rows=0, matched_rows=0, unmatched_rows=0, stringsAsFactors=FALSE)
# Pass 1: compute D_ir from all chromosomes, without chromosome-specific normalization.
for (k in seq_along(marks)) {
  path <- input_path(marks[k]); con <- gzfile(path,'rt'); rows <- 0; matched <- 0
  repeat {
    lines <- readLines(con,n=chunk_size,warn=FALSE)
    if (!length(lines)) break
    x <- parse_lines(lines,path); rows <- rows+length(lines)
    idx <- match(x$barcode,barcodes); keep <- !is.na(idx)
    if (any(keep)) {
      matched <- matched+sum(keep)
      sums <- rowsum(x$multiplicity[keep],idx[keep],reorder=FALSE)
      ii <- as.integer(rownames(sums))
      depths[ii,k] <- depths[ii,k]+as.numeric(sums[,1])
      cnt <- tabulate(idx[keep],nbins=n)
      rowcounts[,k] <- rowcounts[,k]+cnt
    }
  }
  close(con)
  file_audit$total_rows[k] <- rows; file_audit$matched_rows[k] <- matched
  file_audit$unmatched_rows[k] <- rows-matched
  message(marks[k], ': ', matched, ' annotated fragment rows, ',sum(depths[,k]>0),' cells with depth > 0')
}
# Pass 2: emit one weighted atom per original fragment row, only for annotated cells.
# GEO fragment BED intervals are zero-based half-open; midpoint is (start+end)/2.
# Multiplicities are preserved; atom_mass=multiplicity/D_ir.
for (k in seq_along(marks)) {
  path <- input_path(marks[k]); dest <- file.path(out_dir,paste0('E45_S1_',marks[k],'_empirical_atoms.tsv.gz'))
  con <- gzfile(path,'rt'); out <- gzfile(dest,'wt',compression=6)
  writeLines('barcode\tcelltype\tbatch\tchrom\tmidpoint_bp\tmultiplicity\tD_ir\tatom_mass',out)
  repeat {
    lines <- readLines(con,n=chunk_size,warn=FALSE)
    if (!length(lines)) break
    x <- parse_lines(lines,path); idx <- match(x$barcode,barcodes); keep <- !is.na(idx)
    if (!any(keep)) next
    idx <- idx[keep]; mult <- x$multiplicity[keep]; D <- depths[idx,k]
    if (any(D <= 0)) stop('Nonpositive depth for an observed fragment')
    # Use format(..., scientific=FALSE) for barcodes and genomic positions.
    tab <- data.frame(barcode=x$barcode[keep],celltype=unname(celltype[x$barcode[keep]]),
       batch=unname(batch[x$barcode[keep]]),chrom=x$chrom[keep],
       midpoint_bp=(x$start[keep]+x$end[keep])/2,
       multiplicity=mult,D_ir=D,atom_mass=mult/D,check.names=FALSE)
    write.table(
      tab, out,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE,
      na = "NA"
    )
  }
  close(con); close(out)
}
# Per-cell depths, number of rows, audit, and simple depth summary.
res <- data.frame(barcode=barcodes,celltype=unname(celltype[barcodes]),batch=unname(batch[barcodes]),check.names=FALSE)
for (k in seq_along(marks)) res[[paste0('D_',marks[k])]] <- depths[,k]
for (k in seq_along(marks)) res[[paste0('fragment_rows_',marks[k])]] <- rowcounts[,k]
res <- res[order(res$barcode),,drop=FALSE]
write.csv(res,file.path(out_dir,'E45_S1_cell_depths.csv'),row.names=FALSE,na='')
write.csv(file_audit,file.path(out_dir,'E45_S1_file_audit.csv'),row.names=FALSE)
summary <- data.frame(mark=marks, cells_with_positive_depth=colSums(depths>0),
   min_depth=apply(depths,2,min),median_depth=apply(depths,2,median),
   p10=apply(depths,2,quantile,probs=.1,names=FALSE),
   p90=apply(depths,2,quantile,probs=.9,names=FALSE),
   max_depth=apply(depths,2,max),stringsAsFactors=FALSE)
write.csv(summary,file.path(out_dir,'E45_S1_depth_summary.csv'),row.names=FALSE)
# Algebraic mass check: all nonzero cell/mark totals are normalized to 1.
# Additional rounding differences can arise from writing decimal weights to disk.
if (any(depths < 0) || any(rowcounts < 0)) stop('Negative depth/count')
if (any(colSums(rowcounts) != file_audit$matched_rows)) stop('Row count audit mismatch')
message('Completed. Output directory: ',normalizePath(out_dir))
