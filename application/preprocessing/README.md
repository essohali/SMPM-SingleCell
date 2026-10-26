# Preprocessing: scMTR-seq E4.5 S1 empirical point measures

This directory contains `01_build_empirical_measures.R`, the archived base-R
preprocessing script used to construct cell-level depths and empirical genomic
fragment atoms for the SMPM analysis. The repository includes its **prepared
outputs** in `../data/`; it does not include the upstream raw fragment files.

## Source and reference

Wang et al. (2025), *Science Advances*, 11(32), eadu3308,
DOI: https://doi.org/10.1126/sciadv.adu3308.
Source sequencing records: GEO [GSE277856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277856).
The source study's processing resources include Zenodo records
[10.5281/zenodo.15046901](https://doi.org/10.5281/zenodo.15046901)
and [10.5281/zenodo.15046989](https://doi.org/10.5281/zenodo.15046989).
These references identify upstream resources; they do not imply redistribution permission.

## Required local inputs (not bundled)

1. `S2_E45_raw.csv`: E4.5 annotations extracted from supplementary Table S2,
   including `new.cellid`, `sub.lib`, `batch`, `celltype`.
2. A local directory containing six gzip fragment TSV files named:
   `GSM8529560_E45_S1_<MARK>.fragments.tsv.gz`, with `<MARK>` one of
   `H3K27ac`, `H3K27me3`, `H3K36me3`, `H3K4me1`, `H3K4me3`, `H3K9me3`.
   Each row is five tab-separated fields: chromosome, start, end, barcode,
   multiplicity. Check the downloaded files and naming against the source
   archive before execution; they are not supplied here.

## Run from the repository root

```sh
Rscript preprocessing/01_build_empirical_measures.R /path/to/fragments /path/to/S2_E45_raw.csv preprocessing/output
```

Base R only. The script processes each fragment file in 50,000-row chunks
using two passes, and selects `sub.lib == "S1"` and cell types EPI/PE/TE.
For each cell and histone mark, the genome-wide depth `D_ir` is the sum of
fragment multiplicities across **all chromosomes**. Each fragment produces
an atom at `(start + end)/2` with mass `multiplicity / D_ir`. The main
SMPM analysis later restricts atoms to chromosome 1; it does not recompute
chromosome-specific depth. The script builds all six marks, while the
SMPM model uses five (excluding H3K9me3).

## Outputs

- `E45_S1_cell_depths.csv`: 2,015 cells × 15 columns, with lineage, batch,
  six mark depths and fragment-row counts.
- `E45_S1_<MARK>_empirical_atoms.tsv.gz`: one file per mark, six total.
- `E45_S1_file_audit.csv`: total, matched and unmatched fragment rows.
- `E45_S1_depth_summary.csv`: per-mark depth summaries.

## Archived-output provenance checks

The depth file in `../data/` was exactly identical, after R import, to the
archived preprocessing output (2,015 × 15; `identical(A,B) == TRUE`). All five
analysis-mark atom files were compared in full against archived preprocessing
outputs: 5,645,484 lines, same order and categorical fields; genomic
midpoints, multiplicities and depths agreed exactly after numeric parsing.
The largest absolute difference in serialized atom masses was
`3.330669e-16`, below the comparison tolerance of `1e-12`. Compressed file
hashes differed, so byte identity of the gzip files is **not** claimed.

**Limitations:** the preprocessing script has not been rerun independently
from the original GEO fragment files as part of this release audit. The
separate construction/provenance of `../data/S1_multigene_scores_by_cell.csv`
is not documented by this script. The upstream raw reads-to-fragments
pipeline is not reproduced by this script. Verify data redistribution rights
before public hosting.
