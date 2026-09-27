# =============================================================================
# dashboard/data-raw/00_read_targets.R
#
# One reader for the district outcome targets the dashboard fits on, sourced by
# 05_build_protocol_v2_bundles.R and 06_build_uncertainty_ensembles.R so the
# deployment ranking and its uncertainty ensembles see identical data.
#
# targets_v2.csv carries the protocol's 24 cells. The Malawi selenium and
# iodine outcomes (script 61, MW-SE/MW-IO) live in their own table with their
# own shape (p_low = share below the cut-off, n = respondents); they are
# adapted here to the targets_v2 columns the builder reads. They are Malawi-
# only, in-country evidence only (malawi_selenium_iodine_in_country.csv), and
# every consumer labels them so.
# =============================================================================

read_targets_with_extras <- function(p2 = "results/tables/protocol_v2") {
  TG <- read.csv(file.path(p2, "targets_v2.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  f <- file.path(p2, "malawi_selenium_iodine_targets.csv")
  if (file.exists(f)) {
    se <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    se <- se[is.finite(se$p_low) & is.finite(se$n_eff) & se$n_eff > 0, ]
    add <- data.frame(country = "Malawi", outcome = se$outcome,
                      Admin1 = se$Admin1, Admin2 = se$Admin2,
                      y_prev = se$p_low, n_eff = se$n_eff, n_eff_district = se$n_eff,
                      n_raw = se$n, n_psu = NA_real_,
                      y_level = se$mean_log, n_eff_cont = se$n_eff,
                      stringsAsFactors = FALSE)
    for (cc in setdiff(names(TG), names(add))) add[[cc]] <- NA
    TG <- rbind(TG, add[, names(TG)])
  }
  TG
}
