# Syntax-only check: parses R code; does NOT source simulation scripts.
# Run from repository root: source("validation/check_R_syntax.R")
scripts <- c(list.files("R", pattern="\\.[Rr]$", full.names=TRUE),
             list.files(".", pattern="\\.[Rr]$", full.names=TRUE))
stopifnot(length(scripts) >= 13L)
for (f in scripts) {
  tryCatch(parse(file=f), error=function(e) stop("Syntax error in ", f,
                 ": ", conditionMessage(e)))
}
required <- c(paste0("run_E",1:6,".R"),
              paste0("results/tables/E",1:5,"_summary.csv"),
              "results/tables/E6_summary_analytic_W20.csv",
              "results/tables/E6_summary_analytic_W40.csv",
              "results/figures/Fig_E4_supervision_corrected.pdf")
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing expected files: ", paste(missing, collapse=", "))
cat("PASS:", length(scripts), "R scripts parsed; expected files present. No simulation executed.\n")
