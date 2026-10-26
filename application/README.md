# Supervised multivariate principal measures — E4.5 single-cell application

Reproducibility materials for the chromosome 1 (mm10/GRCm38) E4.5
single-cell application in the manuscript **“Supervised multivariate principal
measures for single-cell point processes.”** This package contains four R
analysis/validation scripts, prepared inputs, archived numerical outputs,
Figures 2–3, and a recovered empirical-measure preprocessing script.

**Release status:** analysis-level numerical checks passed; **public release
is pending the author's license decision and confirmation of the right to
redistribute the processed data**. This is not a claim of end-to-end
reproduction from raw sequencing reads or a fresh independent rerun in every
environment.

## Study and provenance

- 2,015 S1 cells: EPI 325, PE 384, TE 1,306; batches E45_1 719 and E45_2 1,296.
- Five analysis marks: H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3.
- Chromosome 1 (mm10/GRCm38), length 195,471,971 bp.
- Upstream study: Wang et al. (2025), *Science Advances* 11(32), eadu3308,
  https://doi.org/10.1126/sciadv.adu3308.
- GEO accession: [GSE277856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277856).
- Related source-study Zenodo records: https://doi.org/10.5281/zenodo.15046901
  (processing resources) and https://doi.org/10.5281/zenodo.15046989
  (figure/code resources).

The `preprocessing/` directory documents the construction of empirical
fragment measures from annotations and six mark-specific fragment files.
The preprocessing script uses genome-wide per-cell depths and produces six
marks, including H3K9me3; the SMPM application uses the five marks above
and selects chromosome 1 at analysis time. Raw fragment inputs and the
upstream FASTQ-to-fragment pipeline are **not bundled**. See
[`preprocessing/README.md`](preprocessing/README.md) for required filenames,
command-line invocation and limitations.

`data/S1_multigene_scores_by_cell.csv` is an ancillary input. Its exact
upstream construction has not yet been documented by the recovered
empirical-measure preprocessing script; do not describe it as independently
reconstructed.

## Repository layout

```text
SMPM_E45_S1_GitHub_release/
├── README.md
├── LICENSE                         # license decision pending
├── .gitignore
├── SHA256SUMS.txt                  # SHA-256 of files, excluding itself
├── preprocessing/
│   ├── 01_build_empirical_measures.R
│   └── README.md
├── data/
│   ├── E45_S1_cell_depths.csv
│   ├── E45_S1_<mark>_empirical_atoms.tsv.gz  # five files
│   └── S1_multigene_scores_by_cell.csv
├── scripts/                        # four R scripts
├── results/                        # archived main numerical outputs
├── figures/                        # Figures 2 and 3 in PDF and PNG
└── validation/
    ├── audit_rerun/                # archived rerun and diagnostic outputs
    ├── consolidation_rerun/        # archived sensitivity summaries
    └── *.csv                       # figure and comparison diagnostics
```

Unvalidated `robustness_exploratory` outputs are deliberately excluded.
The validation files listed above are included and should remain tracked
by Git.

## Requirements

The primary application and preprocessing use base R; `irlba` is optional
for the unsupervised alignment comparison. The figures script requires
`ggplot2`, `dplyr`, `tidyr`, `readr`, and `gridExtra`; `ragg` is optional
for PNG output. The archived R environment is documented in
`results/sessionInfo.txt`.

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "gridExtra"))
```

## Reproduce the analysis

**Run from the repository root**, not from `scripts/`. To avoid overwriting
archived outputs, first work on a copy or choose new output directories.

```sh
Rscript scripts/SMPM_E45_S1_application.R data results
Rscript scripts/SMPM_E45_S1_figures.R
Rscript scripts/SMPM_E45_S1_audit.R data validation/audit_rerun
Rscript scripts/SMPM_E45_S1_consolidation.R data validation/consolidation_rerun
```

Equivalent interactive RStudio usage (working directory at the root):

```r
source("scripts/SMPM_E45_S1_application.R")
source("scripts/SMPM_E45_S1_figures.R")
source("scripts/SMPM_E45_S1_audit.R")
source("scripts/SMPM_E45_S1_consolidation.R")
```

The main script processes large gzip inputs and runs a stratified bootstrap
with `B = 200`. The audit reruns the application under
`validation/audit_rerun/original_rerun/`. The figure script reads
`results/cell_scores.csv` and `results/supervised_eigenfunctions.csv`.
Rerunning any script may replace tracked result files; check Git diffs
before committing a new release.

## Manuscript-to-file mapping

| Output | Reference file |
|---|---|
| Figure 2: supervised cell scores | `figures/Figure2_SMPM_scores.pdf` (PNG also included) |
| Figure 3: signed increments | `figures/Figure3_SMPM_signed_increments.pdf` (PNG also included) |
| Table 11: eigenvalues and explained variation | `results/supervised_eigenvalues.csv` |
| Table 11: bootstrap intervals | `results/bootstrap_eigenvalues.csv` |
| Table 11: mark energies | `validation/audit_rerun/mark_energy.csv` |
| Table 12: leading 10-Mb windows | aggregate `absolute_mass` by component/window across marks in `results/regional_signed_contributions_10Mb.csv` |
| Table 13: mark omission, batch/depth adjustment | `validation/audit_rerun/audit_spectral_sensitivity.csv` |
| Table 13: grid sensitivity | `validation/consolidation_rerun/spectral_and_grid_sensitivity.csv` |
| Regional stability diagnostics | `validation/consolidation_rerun/regional_stability_metrics.csv` |

Publication figures use manuscript numbering Figures **2** and **3**.
Their files were opened and visually inspected for readability; verify
content against the final accepted manuscript before tagging a release.

## Archived validation and numerical comparisons

- Table 11, Table 12 and Table 13 numerical values were checked against the
  reported manuscript precision in the earlier analysis audit.
- Three archived versions of `supervised_eigenfunctions.csv` (original,
  release and audit rerun) were identical after R import (5,005 × 4;
  numerical differences zero).
- The 11 shared non-PDF numerical results under `results/` and
  `validation/audit_rerun/original_rerun/` are byte-identical in the
  archived release.
- The `data/E45_S1_cell_depths.csv` file is exactly identical after R
  import to the archived preprocessing output (2,015 × 15).
- All five empirical-atom files were compared with archived preprocessing
  outputs in their entirety: **5,645,484 fragment rows**. Barcode, lineage,
  batch, chromosome, midpoint, multiplicity and depth agreed; the largest
  absolute `atom_mass` difference was **3.330669 × 10⁻¹⁶** (comparison
  tolerance 10⁻¹²). Compressed-file hashes were not equal because their
  serialized bytes differ; do not claim byte-level equality.

These are archived-output comparisons. The raw-fragment preprocessing was
**not** independently rerun as part of the verification described here.

## Sensitivity analysis interpretation

The reference batch-composition sensitivity averages lineage-by-batch
profiles using pooled batch weights `q = (719/2015, 1296/2015)`.
The joint batch/depth adjustment regresses profile coordinates on three
lineage indicators, five centered `log(1 + depth)` covariates and a centered
batch-2 indicator, without an additional intercept. Spectral subspaces
are compared with canonical cosines on a common grid; signed regional
increments are matched and sign-aligned before regional comparisons.
These observational adjustments are not causal estimates of technical effects.

The corrected grid-increment export treats the initial node as a point
mass at zero and subsequent rows as grid increments. Eigenfunction signs
are indeterminate, and near-degenerate eigenvalues favor subspace over
individual-direction comparisons. Signed regional contributions are
**descriptive**, not inferential genomic peak calls or validated gene
regulatory assignments. Phenotype labels construct the supervised
components; apparent group separation is **not** independent predictive
validation. The four-mark canonical cosine of 1 can arise algebraically
from the three-group, rank-two setting and is not, by itself, evidence
of robustness.

## File integrity

`SHA256SUMS.txt` records SHA-256 hashes of every distributed file other than
itself. Hashes refer to the extracted files, not to the ZIP archive. After
extracting and entering the repository root, check with:

```sh
sha256sum -c SHA256SUMS.txt
```

On Windows, use an equivalent SHA-256 verification tool. Regenerate the
manifest after **any** edit to a tracked file.

## Citation, license and public availability

Please cite the SMPM manuscript when using these materials (full
bibliographic details/DOI to be added after publication), and cite the
upstream Wang et al. (2025) dataset and GEO accession.

The bundled `LICENSE` is an **explicit no-license notice**, not an
open-source permission grant. The author must select a code license and
independently confirm rights to redistribute third-party-derived
processed inputs before making the repository public. A source study's
availability or another pipeline's GPL license does not automatically
license this R code or these processed data. **Do not treat this ZIP as
cleared for public redistribution until these decisions are resolved.**
