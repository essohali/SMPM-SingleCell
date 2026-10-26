# E1–E6 experiment and results index

This file maps the **existing** archived results; no simulations were rerun.

| ID | Design / question | Results in `results/tables/` |
| --- | --- | --- |
| E1 | n=200,500,1000; consistency of eigenvalues, eigendirections and operator | `E1_summary.csv`, `E1_raw_results.rds` |
| E2 | Depths 50,200,1000; fixed and heterogeneous depths; balanced/imbalanced groups | `E2_summary.csv`, `E2_raw_results.rds` |
| E3 | Inference/coverage under sample size, depth and phenotype balance variations | `E3_summary.csv`, `E3_imbalanced_summary.csv`, `E3_raw_results.rds` |
| E4 | Four nuisance amplitudes, 500 saved replicates per amplitude | `E4_replications.csv`, `E4_summary.csv`, `E4_raw_results.rds` |
| E5 | Depths 50,200,1000; SMPM, smoothed, binned; 500 saved replicates per cell | `E5_replications.csv`, `E5_summary.csv`, `E5_results.rds` |
| E6 | Broad/narrow/adjacent; depths 50,200,1000; three methods; W=20/40; 500 saved replicates per cell | `E6_replications_analytic_W*.csv`, `E6_summary_analytic_W*.csv`, `E6_paired_differences_analytic_W*.csv`, `E6_results_analytic_W*.rds` |

The E6 W=20/W=40 summary files correspond to SI Tables S2/S3; paired differences correspond to S4/S5. See [E6 details](E6.md).

## Read-only checks

```bash
python validation/check_saved_results.py
```

The checker recalculates E4–E6 summaries and E6 paired comparisons from saved CSVs. It does **not** rerun the simulations, independently audit E1–E3 RDS files or compare manuscript tables.

In RStudio, the following parses R scripts without executing them:

```r
source("validation/check_R_syntax.R")
```

## Figures and legacy files

`plot_E4_from_saved_results.R` uses E4 archived RDS to produce `results/figures/Fig_E4_supervision_corrected.pdf`. `plot_E3_diagnostics.R` uses archived E3 RDS. These plotting scripts may replace PDFs and are not part of read-only checks. The `results/legacy_E6/` files are retained for historical traceability but do not substantiate the analytic-support SI tables.
