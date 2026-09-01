# =============================================================================
# scripts/accuracy_impact/ws1_reliability.R
#
# WS1a and WS1b. Measure the reliability of a district's survey prevalence, and
# compare the measurement against the analytic ceiling the project has been
# using. See R/reliability_empirical.R for what the split does and does not
# reproduce.
#
#   PROFILE=smoke   Ghana child_iron only, B = 50
#   WS1_B=400       override the number of splits
#
#   Rscript scripts/accuracy_impact/ws1_reliability.R
# -> results/tables/reliability_empirical.csv
# -> results/tables/reliability_analytic_vs_empirical.csv
# =============================================================================
suppressPackageStartupMessages({library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE   <- here("_targets_full")
SUF     <- if (PROFILE == "smoke") "_SMOKE" else ""
B       <- as.integer(Sys.getenv("WS1_B", if (PROFILE == "smoke") "50" else "200"))
SEED    <- 20260901L
TDIR    <- here("results", "tables")
cs <- if (PROFILE == "smoke") "Ghana" else NULL
os <- if (PROFILE == "smoke") "child_iron" else NULL
cat(sprintf("[ws1] profile=%s  B=%d  seed=%d\n", PROFILE, B, SEED))

rel <- build_reliability_empirical(STORE, B = B, seed = SEED,
                                   countries = cs, outcomes = os)
if (is.null(rel)) stop("No reliability rows produced.")

# The design effect that would reconcile the analytic formula with the
# measurement, per cell. This is the diagnostic WS1b turns on.
rel$implied_deff <- NA_real_
for (i in seq_len(nrow(rel))) {
  sfx <- paste0(tolower(gsub(" ", "", rel$country[i])), "_", rel$outcome[i])
  sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", sfx), store = STORE),
                 error = function(e) NULL)
  rel$implied_deff[i] <- round(implied_deff(sv, rel$lambda_emp[i]), 3)
}
readr::write_csv(rel, file.path(TDIR, sprintf("reliability_empirical%s.csv", SUF)))

W <- rel[rel$scheme == "within", , drop = FALSE]
cat("\n=== WS1a: empirical reliability, 'within' scheme ===\n")
print(as.data.frame(W[, c("country","outcome","n_areas_median","r_halfhalf",
                          "lambda_emp","r_max_emp","r_max_emp_lo","r_max_emp_hi",
                          "frac_neg_lambda")]), row.names = FALSE)

cat("\n=== WS1b: analytic against empirical ===\n")
cmp <- W[, c("country","outcome","r_max_analytic","r_max_emp",
             "lambda_analytic","lambda_emp","implied_deff","median_n_analytic")]
cmp$diff <- round(cmp$r_max_emp - cmp$r_max_analytic, 3)
cmp <- cmp[order(cmp$diff), ]
print(as.data.frame(cmp), row.names = FALSE)
readr::write_csv(cmp, file.path(TDIR,
  sprintf("reliability_analytic_vs_empirical%s.csv", SUF)))

fin <- is.finite(cmp$r_max_analytic) & is.finite(cmp$r_max_emp)
cat(sprintf("\ncells compared: %d\n", sum(fin)))
cat(sprintf("median r_max analytic: %.4f   median r_max empirical: %.4f\n",
            stats::median(cmp$r_max_analytic[fin]), stats::median(cmp$r_max_emp[fin])))
cat(sprintf("empirical exceeds analytic in %d of %d cells\n",
            sum(cmp$r_max_emp[fin] > cmp$r_max_analytic[fin]), sum(fin)))
cat(sprintf("analytic r_max == 0 exactly in %d cells; empirical == 0 in %d\n",
            sum(cmp$r_max_analytic[fin] == 0), sum(cmp$r_max_emp[fin] == 0)))
cat(sprintf("median implied design effect: %.3f (the formula assumes 1.5)\n",
            stats::median(cmp$implied_deff[is.finite(cmp$implied_deff)])))

if (any(rel$scheme == "cluster")) {
  C <- rel[rel$scheme == "cluster", , drop = FALSE]
  cat(sprintf("\n=== design-faithful 'cluster' scheme: %d cells computable ===\n", nrow(C)))
  print(as.data.frame(C[, c("country","outcome","n_areas_median","lambda_emp","r_max_emp")]),
        row.names = FALSE)
}
