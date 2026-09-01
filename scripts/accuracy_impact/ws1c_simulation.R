# =============================================================================
# scripts/accuracy_impact/ws1c_simulation.R
#
# WS1c. Bias of the ceiling estimators, measured against a known truth.
#
# The real survey structure is kept (actual districts, clusters and respondents
# per cluster) and only the outcome is simulated, so the sampling geometry is
# preserved while the quantity every estimator targets becomes observable.
# See R/reliability_simulation.R.
#
#   PROFILE=smoke   Ghana child_iron only, R = 20
#   Rscript scripts/accuracy_impact/ws1c_simulation.R
# -> results/tables/reliability_simulation.csv
# =============================================================================
suppressPackageStartupMessages({library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
R_REP <- as.integer(Sys.getenv("WS1C_R", if (PROFILE == "smoke") "20" else "100"))
SEED  <- 20260902L
SD_GRID <- c(0.2, 0.5, 1.0, 1.5)

# Four cells spanning the countries and both the well- and poorly-measured end
# of the reliability range, rather than all 24: the simulation asks about the
# ESTIMATOR, and the estimator does not change from cell to cell.
CELLS <- if (PROFILE == "smoke") list(c("Ghana", "child_iron")) else list(
  c("Ghana", "child_iron"), c("Gambia", "women_iron"),
  c("Malawi", "child_vitA"), c("Sierra Leone", "women_iron"))

cfgs <- get_country_configs()
key <- function(x) tolower(gsub(" ", "", x))
rows <- list()
for (cell in CELLS) {
  cn <- cell[1]; ocn <- cell[2]
  cc <- cfgs[[names(cfgs)[key(names(cfgs)) == key(cn)][1]]]
  oc <- cc$outcomes[[ocn]]; if (is.null(oc)) next
  od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", key(cn), "_", ocn),
                                       store = STORE), error = function(e) NULL)
  if (is.null(od)) { cat("[skip]", cn, ocn, "\n"); next }
  cat(sprintf("[sim] %s %s ...\n", cn, ocn))
  r <- simulate_ceiling_bias(od$data, cc$admin2_col, cc$cluster_id,
                             cc$weight_col, oc$binary,
                             sd_logit = SD_GRID, R = R_REP, deff = 1.5,
                             seed = SEED)
  if (is.null(r)) next
  r$country <- cc$country; r$outcome <- ocn
  rows[[length(rows) + 1L]] <- r
}
if (!length(rows)) stop("No simulation rows produced.")
sim <- do.call(rbind, rows)
front <- c("country", "outcome")
sim <- sim[, c(front, setdiff(names(sim), front))]
readr::write_csv(sim, here("results", "tables",
                           sprintf("reliability_simulation%s.csv", SUF)))

cat("\n=== WS1c: ceiling estimator bias against a known attainable correlation ===\n")
print(as.data.frame(sim[, c("country","outcome","sd_logit","icc","lambda_true",
                            "r_oracle","r_max_analytic","r_max_empirical",
                            "bias_analytic","bias_empirical","pct_analytic_zero")]),
      row.names = FALSE)

cat("\n--- pooled over cells and settings ---\n")
cat(sprintf("mean bias of the ANALYTIC ceiling  : %+.4f\n",
            mean(sim$bias_analytic, na.rm = TRUE)))
cat(sprintf("mean bias of the EMPIRICAL ceiling : %+.4f\n",
            mean(sim$bias_empirical, na.rm = TRUE)))
cat(sprintf("analytic returns exactly zero in %.1f%% of replicates on average\n",
            mean(sim$pct_analytic_zero, na.rm = TRUE)))
