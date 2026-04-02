# =============================================================================
# scripts/export_admin2_predictions.R
#
# Export Admin-2 prediction datasets as clean CSVs for external sharing.
# Combines survey prevalence, individual-level SL predictions, and
# area-level SL predictions (MSE + NLL) into one file per country x outcome.
# =============================================================================

library(targets)
library(here)
library(dplyr)

TARGETS_STORE <- here("_targets_full")

countries <- c("gambia", "ghana", "sierraleone", "malawi")
outcomes  <- c("child_vitA", "women_vitA", "child_iron", "women_iron",
               "women_folate", "women_b12", "child_zinc", "women_zinc")
country_labels <- c(gambia = "Gambia", ghana = "Ghana",
                     sierraleone = "Sierra Leone", malawi = "Malawi")
outcome_labels <- c(
  child_vitA = "Vitamin A deficiency (children)",
  women_vitA = "Vitamin A deficiency (women)",
  child_iron = "Iron deficiency (children)",
  women_iron = "Iron deficiency (women)",
  women_folate = "Folate deficiency (women)",
  women_b12 = "Vitamin B12 deficiency (women)",
  child_zinc = "Zinc deficiency (children)",
  women_zinc = "Zinc deficiency (women)")

load_target <- function(name) {
  tryCatch(targets::tar_read_raw(name, store = TARGETS_STORE),
           error = function(e) NULL)
}

out_dir <- here("results", "shared")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_exported <- 0

for (ctry in countries) {
  for (oc in outcomes) {
    suffix <- paste0(ctry, "_", oc)

    # Survey prevalence
    svy <- load_target(paste0("svy_admin2_", suffix))
    if (is.null(svy) || nrow(svy) == 0) next

    # Individual-level SL aggregated to Admin-2
    sl <- load_target(paste0("admin2_sl_", suffix))

    # Area-level comparison (MSE + NLL)
    area <- load_target(paste0("area_comparison_", suffix))

    # GEE admin2 (for Admin1 names)
    gee <- load_target(paste0("gee_admin2_", ctry))

    # Build output
    out <- svy |>
      select(Admin2,
             svy_prev,
             svy_prev_se = any_of("svy_prev_se"),
             n_svy = any_of("n_svy"))

    # Add Admin1 from GEE data
    if (!is.null(gee) && "Admin1" %in% names(gee)) {
      out <- out |>
        left_join(gee |> select(Admin1, Admin2) |> distinct(),
                  by = "Admin2")
    }

    # Add individual SL predictions
    if (!is.null(sl) && "sl_prev" %in% names(sl)) {
      out <- out |> left_join(sl |> select(Admin2, sl_prev), by = "Admin2")
    }

    # Add area-level predictions
    if (!is.null(area)) {
      if (!is.null(area$area_mse) && !is.null(area$area_mse$area_preds)) {
        out <- out |>
          left_join(area$area_mse$area_preds |>
                      select(Admin2, area_mse_prev = area_pred),
                    by = "Admin2")
      }
      if (!is.null(area$area_nll) && !is.null(area$area_nll$area_preds)) {
        out <- out |>
          left_join(area$area_nll$area_preds |>
                      select(Admin2, area_nll_prev = area_pred),
                    by = "Admin2")
      }
    }

    # Add metadata
    out$country <- country_labels[ctry]
    out$outcome <- oc
    out$outcome_label <- outcome_labels[oc]
    out$has_survey <- TRUE

    # Reorder columns
    col_order <- intersect(
      c("country", "Admin1", "Admin2", "outcome", "outcome_label",
        "svy_prev", "svy_prev_se", "n_svy", "sl_prev",
        "area_mse_prev", "area_nll_prev", "has_survey"),
      names(out)
    )
    out <- out[, col_order]

    # Round numeric columns
    num_cols <- names(out)[sapply(out, is.numeric)]
    out <- out |> mutate(across(all_of(num_cols), ~round(., 4)))

    outfile <- file.path(out_dir, paste0("admin2_predictions_", suffix, ".csv"))
    write.csv(out, outfile, row.names = FALSE)
    n_exported <- n_exported + 1
    cat(sprintf("  Exported %s: %d Admin-2 areas\n", basename(outfile), nrow(out)))
  }
}

cat(sprintf("\nExported %d files to %s\n", n_exported, out_dir))
