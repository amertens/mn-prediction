# =============================================================================
# scripts/accuracy_impact/wsd_calibration_gate.R
#
# WS-D. Apply the calibration gate retroactively to every committed arm.
#
#   Rscript scripts/accuracy_impact/wsd_calibration_gate.R
# -> results/tables/calibration_gate_report.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")

# The implied-shift audit already recorded, per cell, the un-anchored
# population-weighted national aggregate and the design-based survey estimate.
# That is exactly what the gate needs, so the retroactive pass reads it rather
# than refitting.
sh <- read.csv(file.path(TDIR, "anchor_implied_shifts_B1.csv"), stringsAsFactors = FALSE)
sh <- sh[sh$arm == "national", , drop = FALSE]

rows <- list()
for (i in seq_len(nrow(sh))) {
  cn <- sh$country[i]; ocn <- sh$outcome[i]
  # the survey's own SE for the national estimate, for the CI-relative term
  se <- NA_real_
  od <- tryCatch(targets::tar_read_raw(
    paste0("outcome_data_", tolower(cn), "_", ocn), store = STORE),
    error = function(e) NULL)
  cfg <- get_country_configs()
  cck <- names(cfg)[tolower(gsub("[^a-z]", "", tolower(names(cfg)))) ==
                    tolower(gsub("[^a-z]", "", tolower(cn)))]
  if (length(cck) && !is.null(od)) {
    cc <- cfg[[cck[1]]]; oc <- cc$outcomes[[ocn]]
    nat <- tryCatch(national_design_based(od, cc, oc), error = function(e) NULL)
    if (!is.null(nat) && !is.null(nat$n) && is.finite(nat$prev) && nat$n > 1)
      se <- sqrt(nat$prev * (1 - nat$prev) / nat$n)   # design effect not applied
  }
  g <- calibration_gate(pred = sh$before[i], pop = NULL,
                        national_prev = sh$target[i], national_se = se)
  g$country <- cn; g$outcome <- ocn
  g$arm <- "admin-2 ridge (LORO), un-anchored"
  rows[[i]] <- g
}
rep <- dplyr::bind_rows(rows)
front <- c("country", "outcome", "arm", "status")
rep <- rep[, c(front, setdiff(names(rep), front))]
readr::write_csv(rep, file.path(TDIR, "calibration_gate_report.csv"))

cat("=== WS-D: calibration gate, applied retroactively ===\n")
print(as.data.frame(table(rep$status)), row.names = FALSE)
cat("\n--- cells that FAIL the gate ---\n")
f <- rep[rep$status == "calibration_failed", ]
print(as.data.frame(f[order(-f$abs_gap_pp),
  c("country","outcome","national_pred_pp","national_survey_pp","abs_gap_pp","threshold_pp")]),
  row.names = FALSE)
cat(sprintf("\n%d of %d cells fail. Sierra Leone child_vitA present: %s\n",
            nrow(f), nrow(rep),
            any(grepl("[Ss]ierra", f$country) & f$outcome == "child_vitA")))
