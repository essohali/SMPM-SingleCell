# Supervised multivariate principal measures for single-cell point processes

Code and archived Monte Carlo results for **“Supervised multivariate principal measures for single-cell point processes”** by Essoham Ali, Adewale Loukman and Ibukun A. Akinyemi.

This repository contains experiments **E1–E6**. The experiment previously labeled *E5-B* is now called **E6**. **No Monte Carlo replications were rerun to prepare this repository.**

## Read-only verification (no simulations)

From the repository root:

```bash
python validation/check_saved_results.py
```

Requires Python 3 only (standard library). It recalculates available archived CSV summary statistics for **E4, E5, and E6** and checks replicate counts and paired comparisons. It does not execute R or overwrite results.

For a *syntax-only* R check, use R/RStudio from the repository root:

```r
source("validation/check_R_syntax.R")
```

This calls `parse()` on source files and checks expected file presence. It does not execute any Monte Carlo experiment.

**Audit scope:** E1–E6 were reviewed individually before consolidation. These tools do not repeat that scientific audit; the CSV checker does not independently decode or recalculate E1–E3 raw RDS objects or prove full numerical reproducibility.

## Experiments

| Experiment | Focus | Archived outputs |
| --- | --- | --- |
| E1 | Consistency as sample size increases | `E1_summary.csv`, `E1_raw_results.rds` |
| E2 | Fixed/heterogeneous depth and group balance | `E2_summary.csv`, `E2_raw_results.rds` |
| E3 | Eigenvalue inference and interval diagnostics | `E3_summary.csv`, `E3_imbalanced_summary.csv`, `E3_raw_results.rds` |
| E4 | Supervised versus unsupervised directions | `E4_replications.csv`, `E4_summary.csv`, corrected figure |
| E5 | Direct, smoothed and binned representations | `E5_replications.csv`, `E5_summary.csv`, `E5_results.rds` |
| E6 | Localization for broad/narrow/adjacent signals at W=20/40 | `E6_*_analytic_W20/W40.csv` and RDS |

All listed tabular files reside under `results/tables/`. See [experiment index](docs/EXPERIMENTS.md) and [E6 details](docs/E6.md).

## Directory layout

```text
R/                           Shared engine and E1–E6 simulation implementations
run_E1.R ... run_E6.R        Simulation entry points; NOT needed for read-only checks
results/tables/              Archived CSV and RDS outputs
results/figures/             Archived figures and E3 diagnostics
results/legacy_E6/           Older E6 outputs, not used for analytic-support SI tables
validation/                  Read-only CSV and R source-syntax checks
docs/                        Experiment index, provenance and release notes
plot_E4_from_saved_results.R  E4 figure from archived RDS
plot_E3_diagnostics.R         E3 diagnostics from archived RDS
```

## Requirements and reproduction precautions

- CSV verification: **Python 3 standard library**.
- R syntax check and original simulation scripts: **base R**.
- The simulation entry points `run_E1.R` ... `run_E6.R` **write into `results/` and may overwrite the archived files**. Run them only in a separate working copy if a full rerun is desired; a rerun is *not required* to inspect the published results.
- Figure-generation scripts read archived RDS but may overwrite figure PDFs. They are not part of read-only verification.
- Work from the repository root for all commands. E1/E2 use `SMPM_B`, E3 uses `SMPM_E3_B` and `SMPM_E3_IMBALANCED`, and E4/E5/E6 have experiment-specific environment variables. E6 selects resolution using `SMPM_E6_WINDOWS` (20 or 40).

## Scientific and provenance notes

The E6 score is a **sum of Euclidean norms of consecutive increment vectors** (discrete vector total variation). The square root is **inside** the sum over increments, matching manuscript Eq. (36), corrected SI Eq. (S6), and `window_tv()`. E6 compares fixed smoothing and ten-bin binning competitors; the top-K F1 uses the oracle number of active windows. See [E6 documentation](docs/E6.md).

`docs/SOURCE_MANIFEST.csv` is the **historical original-source** manifest; its SHA-256 hashes describe files in the original archives, not the current edited scripts. `docs/ARCHIVED_RESULTS_SHA256.csv` gives hashes of the packaged archived results to facilitate byte-preservation checks. Packaging did not establish that every saved result was generated with the *current* script revision.

## Citation and license

Please cite the accompanying manuscript. A final manuscript DOI and an explicit software license were not supplied; no license is assumed. Authors should add the final citation and select a license before public release.
