# =============================================================================
# scripts/accuracy_impact/wsk1_model_anchor.R
#
# WS-K. Does a MODEL-PREDICTED national anchor rescue the no-survey case?
#
# THE QUESTION
# ------------
# WS-I established that under leave-one-country-out the level is everything and
# the covariate ranking is nearly free:
#
#   anchor = the country's TRUE national prevalence     district MAE  8.94 pp
#   anchor = the training countries' mean ("blind")     district MAE 17.55 pp
#
# The blind anchor is the weakest honest one, and it carries 14.18 pp of bias.
# Everything a country without a survey could ever get sits between those two
# numbers, and nothing in this project has measured where.
#
# The national VMNIS panel offers a third option. Its unit of observation is the
# country-survey, so a 69-country panel is a completely different estimation
# problem from 14 to 87 districts, and it performs accordingly: leave-one-
# country-out national MAE of 11.75 pp for vitamin A in preschool children
# (r 0.655, 108 surveys) against a null of 16.23.
#
# So: substitute that model's prediction for the blind anchor and re-score the
# districts. If 11.75 pp of national error buys most of the distance from 17.55
# down toward 8.94, a country with no survey has a usable product. If it does
# not, the honest answer is that the anchor has to come from data collection.
#
# WHAT THIS CAN AND CANNOT COVER
# ------------------------------
# VMNIS has no iron or ferritin panel, so child_iron and women_iron cannot be
# anchored this way at all. The overlap between NAT_PANEL_MAP and the outcomes
# the LOCO decomposition scores is vitamin A only, in two populations. Eight
# cells, and the script says so rather than implying wider coverage.
#
#   Rscript scripts/accuracy_impact/wsk1_model_anchor.R
# -> results/tables/model_anchor_loco.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results", "tables")
TILTS <- c(0, 0.2, 0.35)
SEED <- 20260928L
ISO <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
SURVEY_YEAR <- c(Gambia = 2018, Ghana = 2017, Malawi = 2016, SierraLeone = 2013)

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()

# The VMNIS panel and its national covariates, loaded once. Names taken from
# scripts/covariates/19_national_composition.R, which is the working caller;
# an earlier version of this script invented a loader that does not exist and
# the tryCatch turned every anchor into a silent NA.
NAT <- vmnis_national()
NATCOV <- load_national_covariates(Sys.getenv("NAT_COV", "wdi"))
NATV <- NATCOV$vars; NATLOGV <- NATCOV$log_vars
NATDAT <- dplyr::inner_join(
  NAT, NATCOV$df |> dplyr::select(iso3c, year, dplyr::all_of(NATV)),
  by = c("iso3c", "year"))
cat(sprintf("[wsk1] VMNIS panel %d rows / %d countries, %d covariates
",
            nrow(NATDAT), length(unique(NATDAT$iso3c)), length(NATV)))

rows <- list()
for (ocn in names(NAT_PANEL_MAP)) {
  # district data for this outcome, pooled
  dd <- list()
  for (cn in names(cfgs)) {
    cc <- cfgs[[cn]]
    if (!ocn %in% names(cc$outcomes)) next
    sv <- tryCatch(targets::tar_read_raw(
      paste0("svy_admin2_", tolower(cn), "_", ocn), store = STORE),
      error = function(e) NULL)
    if (is.null(sv) || !nrow(sv)) next
    hc <- H[H$country == cn, , drop = FALSE]; if (!nrow(hc)) next
    m <- dplyr::inner_join(sv, hc, by = admin2_join_by(sv, hc))
    if (nrow(m) < 8) next
    m$country <- cc$country; m$cn_key <- cn
    dd[[cn]] <- m
  }
  if (length(dd) < 3) {
    cat(sprintf("[wsk1] %-13s skipped: %d countries\n", ocn, length(dd))); next
  }
  d <- dplyr::bind_rows(dd)
  usable <- COVS[vapply(COVS, function(v) {
    x <- d[[v]]
    all(vapply(split(x, d$country), function(z)
      mean(is.finite(z)) > 0.8 && stats::sd(z, na.rm = TRUE) > 0, logical(1)))
  }, logical(1))]
  if (length(usable) < 5) next

  for (cn in unique(d$cn_key)) {
    te <- d[d$cn_key == cn, , drop = FALSE]
    tr <- d[d$cn_key != cn, , drop = FALSE]
    if (nrow(te) < 8 || length(unique(tr$country)) < 2) next
    Xtr <- as.matrix(tr[, usable, drop = FALSE]); Xte <- as.matrix(te[, usable, drop = FALSE])
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE)
      Xtr[!is.finite(Xtr[, j]), j] <- mu; Xte[!is.finite(Xte[, j]), j] <- mu
    }
    set.seed(SEED)
    f <- tryCatch(.ds_fit(Xtr, tr$svy_prev, k_screen = min(20L, ncol(Xtr))),
                  error = function(e) NULL)
    if (is.null(f)) next
    pred <- tryCatch(.ds_predict(f, Xte), error = function(e) NULL)
    if (is.null(pred)) next

    a_true  <- stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1), na.rm = TRUE)
    a_blind <- stats::weighted.mean(tr$svy_prev, pmax(tr$n_svy, 1), na.rm = TRUE)
    spread <- loco_training_spread(tr$svy_prev, tr$country)
    if (!is.finite(spread)) spread <- stats::sd(tr$svy_prev, na.rm = TRUE)

    # THE MODEL ANCHOR. Fitted on the VMNIS panel with this country excluded,
    # so it is a genuine out-of-country prediction of the national level.
    a_model <- tryCatch({
      pn <- national_panel_for(ocn)
      dpanel <- NATDAT[NATDAT$mn_group == pn[1] & NATDAT$pop == pn[2], , drop = FALSE]
      if (!nrow(dpanel)) stop("empty panel")
      predict_national_level(dpanel, NATV, NATLOGV,
                             iso3c = ISO[[cn]], year = SURVEY_YEAR[[cn]],
                             cov_df = NATCOV$df)
    }, error = function(e) { message("  [anchor] ", ocn, " ", cn, ": ",
                                     conditionMessage(e)); NA_real_ })

    for (nm in c("true", "blind", "model")) {
      a <- switch(nm, true = a_true, blind = a_blind, model = a_model)
      if (!is.finite(a)) next
      for (tl in TILTS) {
        est <- loco_anchor_tilt(a, pred, spread, tilt = tl)
        o <- is.finite(te$svy_prev) & is.finite(est)
        if (sum(o) < 4) next
        rows[[length(rows) + 1L]] <- data.frame(
          outcome = ocn, held_out = cn, anchor = nm, tilt = tl,
          n = sum(o),
          anchor_pp = round(100 * a, 2),
          anchor_err_pp = round(100 * (a - a_true), 2),
          mae_pp = round(100 * mean(abs(est[o] - te$svy_prev[o])), 2),
          bias_pp = round(100 * mean(est[o] - te$svy_prev[o]), 2),
          stringsAsFactors = FALSE)
      }
    }
    cat(sprintf("[wsk1] %-13s %-12s true=%5.1f blind=%5.1f model=%s\n",
                ocn, cn, 100 * a_true, 100 * a_blind,
                if (is.finite(a_model)) sprintf("%5.1f", 100 * a_model) else "   NA"))
  }
}

if (!length(rows)) stop("[wsk1] nothing scored")
out <- dplyr::bind_rows(rows)
readr::write_csv(out, file.path(TDIR, "model_anchor_loco.csv"))

cat("\n=== district error by anchor source, pooled ===\n")
s <- out |> group_by(anchor, tilt) |>
  summarise(cells = dplyr::n(),
            anchor_abs_err = round(mean(abs(anchor_err_pp)), 2),
            mae_pp = round(mean(mae_pp), 2),
            abs_bias = round(mean(abs(bias_pp)), 2), .groups = "drop") |>
  arrange(tilt, mae_pp)
print(as.data.frame(s), row.names = FALSE)

cat("\n=== how much of the gap from blind to true does the model anchor close? ===\n")
z <- s[s$tilt == 0, ]
if (all(c("true", "blind", "model") %in% z$anchor)) {
  tr <- z$mae_pp[z$anchor == "true"]; bl <- z$mae_pp[z$anchor == "blind"]
  md <- z$mae_pp[z$anchor == "model"]
  cat(sprintf("  blind %.2f -> model %.2f -> true %.2f pp\n", bl, md, tr))
  cat(sprintf("  gap closed: %.0f%%\n", 100 * (bl - md) / (bl - tr)))
}
cat(sprintf("\n-> %s\n", file.path("results", "tables", "model_anchor_loco.csv")))
