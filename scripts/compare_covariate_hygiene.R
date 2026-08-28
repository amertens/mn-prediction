# =============================================================================
# scripts/compare_covariate_hygiene.R
#
# Does GEE covariate hygiene help, hurt, or make no difference?
#
# Runs the area-level LOCO transport recipe twice on the SAME cached Admin-2
# targets -- once with the current covariate set (v1) and once with the
# cross-band _annual_* summaries over non-commensurable bands and the static
# layers' per-year duplicates removed (v2, R/gee_band_semantics.R) -- and reports
# both side by side.
#
# This is the evidence for whether GEE_COVARIATE_HYGIENE should default to on.
# It reads the {targets} store and refits only the cheap elastic-net recipe, so
# it does NOT require a pipeline rebuild.
#
# Run:
#   Rscript scripts/compare_covariate_hygiene.R
#
# Writes results/sensitivity/covariate_hygiene_comparison_<scheme>.csv, tagged with a
# `scheme` column. The original 2026-08-18 covariate_hygiene_comparison.csv is a
# frozen baseline (results/tables/frozen_2026-08/) and is never overwritten;
# override the scheme with HYGIENE_SCHEME if you need a second run side by side.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})
# tar_source() rather than a hand-listed set: build_area_loco_dataset() now
# calls the Admin-2 key helpers (admin2_join_by, warn_if_join_multiplied,
# report_pair_join_losses), and a partial source list silently loses them.
targets::tar_source(here("R"))

OUTCOMES <- c("child_vitA", "women_vitA", "child_iron", "women_iron")

# Additive-results convention: every run is tagged and lands in its own file, so
# no previously published comparison is regenerated in place.
SCHEME   <- Sys.getenv("HYGIENE_SCHEME", "ws2b_2026-08")
OUT_PATH <- here("results", "sensitivity",
                 sprintf("covariate_hygiene_comparison_%s.csv", SCHEME))
if (basename(OUT_PATH) == "covariate_hygiene_comparison.csv")
  stop("Refusing to overwrite the frozen 2026-08-18 baseline. Pick a scheme id.")

read_target <- function(nm) tryCatch(targets::tar_read_raw(nm), error = function(e) NULL)

countries <- names(get_country_configs())
gee <- list()
for (cn in countries) {
  g <- read_target(paste0("gee_admin2_", tolower(cn)))
  if (!is.null(g) && ncol(g) > 2) gee[[cn]] <- g
}
message(sprintf("countries with cached GEE covariates: %s",
                paste(names(gee), collapse = ", ")))
if (length(gee) < 2) stop("Need >= 2 countries with gee_admin2_* in the store.")

#' Pooled area dataset in the shape run_area_transport_loco() expects.
build_pool <- function(otag) {
  svy <- list()
  for (cn in names(gee)) {
    s <- read_target(paste0("svy_admin2_", tolower(cn), "_", otag))
    if (!is.null(s) && nrow(s) > 0) svy[[cn]] <- s
  }
  if (length(svy) < 2) return(NULL)
  # build_area_loco_dataset() applies the hygiene pruning internally (it is a
  # no-op when GEE_COVARIATE_HYGIENE=false), so toggling the env var around
  # this call is what produces v1 vs v2.
  p <- build_area_loco_dataset(svy, gee[names(svy)])
  d <- p$pooled_data
  if (is.null(d) || !nrow(d)) return(NULL)
  if (!"n_svy" %in% names(d)) d$n_svy <- 1
  list(data = d, predictors = p$common_gee_vars,
       countries = p$country_names, outcome = otag)
}

run_variant <- function(otag, hygiene) {
  Sys.setenv(GEE_COVARIATE_HYGIENE = if (hygiene) "true" else "false")
  pool <- build_pool(otag)
  if (is.null(pool) || length(pool$predictors) < 5) return(NULL)
  res <- tryCatch(run_area_transport_loco(pool, AREA_TRANSPORT_RECIPE),
                  error = function(e) { message("  fit failed: ", conditionMessage(e)); NULL })
  if (is.null(res) || is.null(res$metrics) || !nrow(res$metrics)) return(NULL)
  m <- res$metrics
  m$variant     <- if (hygiene) "v2_hygiene" else "v1_current"
  m$outcome     <- otag
  m$n_predictors <- length(pool$predictors)
  m$scheme      <- SCHEME
  m
}

rows <- list()
for (otag in OUTCOMES) {
  message("\n=== ", otag, " ===")
  for (hy in c(FALSE, TRUE)) {
    message("-- ", if (hy) "v2 (hygiene)" else "v1 (current)")
    r <- run_variant(otag, hy)
    if (!is.null(r)) rows[[paste(otag, hy)]] <- r
  }
}
Sys.unsetenv("GEE_COVARIATE_HYGIENE")

if (!length(rows)) stop("No comparable results produced.")
all <- dplyr::bind_rows(rows)

out <- OUT_PATH
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(all, out)

num <- intersect(c("pearson_r", "spearman_r", "mae_pp", "rmse_pp", "n_selected"),
                 names(all))
summ <- all |>
  dplyr::group_by(outcome, variant) |>
  dplyr::summarise(n_predictors = dplyr::first(n_predictors),
                   dplyr::across(dplyr::all_of(num), ~ round(mean(.x, na.rm = TRUE), 3)),
                   .groups = "drop") |>
  dplyr::arrange(outcome, variant)

cat("\n================ mean across held-out countries ================\n")
print(as.data.frame(summ), row.names = FALSE)

if (all(c("pearson_r") %in% names(summ))) {
  w <- summ |>
    dplyr::select(outcome, variant, pearson_r) |>
    tidyr::pivot_wider(names_from = variant, values_from = pearson_r)
  if (all(c("v1_current", "v2_hygiene") %in% names(w))) {
    w$delta <- round(w$v2_hygiene - w$v1_current, 3)
    cat("\nchange in mean LOCO Pearson r (v2 - v1):\n")
    print(as.data.frame(w), row.names = FALSE)
    cat(sprintf("\noverall mean delta: %+.3f\n", mean(w$delta, na.rm = TRUE)))
  }
}
cat("\nwrote ", out, "\n", sep = "")
