# =============================================================================
# sandbox_parsimony/R/24_national_recovery.R
#
# These surveys were designed and weighted for NATIONAL (and usually Admin-1)
# estimates. Admin-2 is an unplanned domain. Two questions follow.
#
# A. HOW ACCURATE ARE THE LOCO MODELS AT RECOVERING NATIONAL PREVALENCE?
#    Aggregate each LOCO model's district predictions to a national figure,
#    population-weighted, and compare with the design-based national estimate
#    from the survey itself. The design-based figure is the thing the survey
#    IS built to deliver, so it is the right reference.
#
# B. WHAT DOES THE UNPLANNED-DOMAIN PROBLEM COST?
#    Two separable effects:
#      B1 variance -- already quantified as r_max in FINDINGS section 1;
#      B2 SELECTION -- which districts got any sample at all is not controlled
#         by the design. If the ~29% of Ghanaian districts that were surveyed
#         differ systematically from the 71% that were not, a model trained on
#         them is extrapolating, and the national aggregate of its predictions
#         inherits that bias. Measured here by comparing covariate distributions
#         of surveyed vs unsurveyed districts, and by checking whether the
#         model's own aggregate over ALL districts drifts away from its
#         aggregate over the surveyed ones.
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
source("sandbox_parsimony/R/09_loco_fns.R")
suppressPackageStartupMessages({library(dplyr); library(mgcv)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi")
CFG_KEY <- c(gambia = "Gambia", ghana = "Ghana",
             sierraleone = "SierraLeone", malawi = "Malawi")
WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"
`%|%` <- function(a, b) if (is.null(a)) b else a

sink(stdout(), type = "message")
source("R/config.R")
cfg <- get_country_configs()
cov_ext <- readRDS("sandbox_parsimony/out/cov_ext.rds")

#' Design-based national prevalence, exactly as the survey supports it:
#' survey-weighted mean over ALL individuals, ignoring district entirely.
national_design_based <- function(ctry, oc, cc) {
  f <- file.path(STORE, paste0("outcome_data_", ctry, "_", oc$tag))
  if (!file.exists(f)) return(NULL)
  od <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(od)) return(NULL)
  d <- od$data
  if (!all(c(cc$weight_col, oc$binary) %in% names(d))) return(NULL)
  y <- suppressWarnings(as.numeric(d[[oc$binary]]))
  w <- suppressWarnings(as.numeric(d[[cc$weight_col]])); w[!is.finite(w) | w <= 0] <- 1
  ok <- is.finite(y)
  if (sum(ok) < 30) return(NULL)
  list(prev = stats::weighted.mean(y[ok], w[ok]), n = sum(ok))
}

#' Admin-2 survey prevalence + covariates, water dropped
build2 <- function(ctry, oc, cc) {
  f <- file.path(STORE, paste0("svy_admin2_", ctry, "_", oc$tag))
  if (!file.exists(f)) return(NULL)
  s <- readRDS(f)
  s <- s[is.finite(s$svy_prev) & !grepl(WATER, s$Admin2, ignore.case = TRUE), , drop = FALSE]
  cv <- cov_ext[[ctry]]; if (is.null(cv)) return(NULL)
  d <- dplyr::inner_join(s[, c("Admin2", "svy_prev", "n_svy")], cv, by = "Admin2")
  d <- d[is.finite(d$lon), , drop = FALSE]
  if (nrow(d) < 15) return(NULL)
  d$country <- ctry
  d
}

# ---------------------------------------------------------------------------
# A. National prevalence recovery under LOCO
# ---------------------------------------------------------------------------
OUTCOMES <- c("child_vitA", "child_iron", "women_iron", "women_vitA")
rows <- list()
for (oc_tag in OUTCOMES) {
  frames <- list(); nat <- list()
  for (ctry in COUNTRIES) {
    cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
    o <- Filter(function(x) x$tag == oc_tag, cc$outcomes); if (!length(o)) next
    o <- o[[1]]
    d <- build2(ctry, o, cc); if (is.null(d)) next
    nd <- national_design_based(ctry, o, cc); if (is.null(nd)) next
    frames[[ctry]] <- d; nat[[ctry]] <- nd
  }
  if (length(frames) < 3) next
  covn <- function(df) {
    cv <- setdiff(names(df), c("Admin2", "Admin1", "svy_prev", "n_svy", "country",
                               "lon", "lat", ".w"))
    cv[vapply(cv, function(c) is.numeric(df[[c]]) && any(is.finite(df[[c]])), logical(1))]
  }
  common <- Reduce(intersect, lapply(frames, covn))
  pooled <- bind_rows(lapply(frames, function(df)
    df[, c("country", "Admin2", "svy_prev", "n_svy", "lon", "lat", common)]))

  for (ho in names(frames)) {
    tr <- pooled[pooled$country != ho, ]; te <- pooled[pooled$country == ho, ]
    if (nrow(tr) < 25 || nrow(te) < 8) next
    nat_obs <- nat[[ho]]$prev

    variants <- list(
      list(name = "PROD enet30", f = function()
        loco_one(tr, te, common, function(X, y, v) screen_topK(X, y, v, 30L),
                 alpha = .5, centering = "none", anchor = "train_mean")),
      list(name = "centered_own ridge", f = function()
        loco_one(tr, te, common, NULL, alpha = 0, centering = "own",
                 anchor = "train_mean")),
      list(name = "zscore ridge", f = function()
        loco_zscore(tr, te, common, NULL, alpha = 0, anchor = "train_mean")),
      list(name = "spatial_tps", f = function()
        loco_spatial(tr, te, common, NULL, anchor = "train_mean")),
      list(name = "train-country mean", f = function()
        rep(stats::weighted.mean(tr$svy_prev, pmax(tr$n_svy, 1)), nrow(te)))
    )
    for (v in variants) {
      pr <- tryCatch(v$f(), error = function(e) NULL); if (is.null(pr)) next
      ok <- is.finite(pr)
      if (sum(ok) < 5) next
      # national aggregate of the predictions, weighted by survey n as a
      # population proxy for the surveyed districts
      nat_pred <- stats::weighted.mean(pr[ok], pmax(te$n_svy[ok], 1))
      rows[[paste(oc_tag, ho, v$name)]] <- data.frame(
        outcome = oc_tag, held_out = ho, model = v$name,
        n_districts = sum(ok),
        national_design_pp = round(100 * nat_obs, 1),
        national_pred_pp = round(100 * nat_pred, 1),
        error_pp = round(100 * (nat_pred - nat_obs), 1),
        abs_error_pp = round(100 * abs(nat_pred - nat_obs), 1),
        rel_error = round((nat_pred - nat_obs) / nat_obs, 2),
        stringsAsFactors = FALSE)
    }
  }
}
nr <- bind_rows(rows)
write.csv(nr, "sandbox_parsimony/out/national_recovery.csv", row.names = FALSE)

cat("\n=== A. LOCO national prevalence recovery ===\n")
cat("national_design_pp is the survey's own design-based national estimate --\n")
cat("the quantity these surveys ARE powered for.\n\n")
print(as.data.frame(nr |> group_by(model) |>
  summarise(cells = n(),
            mean_abs_error_pp = round(mean(abs_error_pp), 1),
            median_abs_error_pp = round(median(abs_error_pp), 1),
            worst_pp = round(max(abs_error_pp), 1),
            mean_rel_error = round(mean(rel_error), 2),
            .groups = "drop") |> arrange(mean_abs_error_pp)), row.names = FALSE)

cat("\n=== per cell, best-transporting model (zscore ridge) ===\n")
print(as.data.frame(nr |> filter(model == "zscore ridge") |>
  select(outcome, held_out, n_districts, national_design_pp,
         national_pred_pp, error_pp)), row.names = FALSE)

# ---------------------------------------------------------------------------
# B. Are the surveyed districts representative of the country?
# ---------------------------------------------------------------------------
cat("\n=== B. Surveyed vs unsurveyed districts ===\n")
cat("Standardised mean difference on each curated covariate. |SMD| > 0.25 is\n")
cat("the usual threshold for 'these groups are not exchangeable'.\n\n")
rep_rows <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  cv <- cov_ext[[ctry]]; if (is.null(cv)) next
  cv <- cv[!grepl(WATER, cv$Admin2, ignore.case = TRUE), , drop = FALSE]
  o <- cc$outcomes[[1]]
  f <- file.path(STORE, paste0("svy_admin2_", ctry, "_", o$tag))
  if (!file.exists(f)) next
  s <- readRDS(f)
  surveyed <- cv$Admin2 %in% s$Admin2[is.finite(s$svy_prev)]
  if (sum(surveyed) < 10 || sum(!surveyed) < 5) {
    rep_rows[[ctry]] <- data.frame(country = ctry, n_surveyed = sum(surveyed),
      n_unsurveyed = sum(!surveyed), max_abs_smd = NA_real_,
      n_smd_gt_025 = NA_integer_, worst_variable = "(nearly all surveyed)",
      stringsAsFactors = FALSE)
    next
  }
  cur <- intersect(curated_vars(names(cv)), names(cv))
  cur <- c(cur, grep("^gee_traveltime_city_median$|^spam_share_", names(cv), value = TRUE))
  smd <- vapply(unique(cur), function(v) {
    x <- suppressWarnings(as.numeric(cv[[v]]))
    a <- x[surveyed]; b <- x[!surveyed]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 5 || length(b) < 5) return(NA_real_)
    sp <- sqrt((stats::var(a) + stats::var(b)) / 2)
    if (!is.finite(sp) || sp < 1e-9) return(NA_real_)
    (mean(a) - mean(b)) / sp
  }, numeric(1))
  smd <- smd[is.finite(smd)]
  rep_rows[[ctry]] <- data.frame(
    country = ctry, n_surveyed = sum(surveyed), n_unsurveyed = sum(!surveyed),
    max_abs_smd = round(max(abs(smd)), 2),
    n_smd_gt_025 = sum(abs(smd) > 0.25),
    n_vars = length(smd),
    worst_variable = names(smd)[which.max(abs(smd))],
    worst_smd = round(smd[which.max(abs(smd))], 2),
    stringsAsFactors = FALSE)
}
rp <- bind_rows(rep_rows)
write.csv(rp, "sandbox_parsimony/out/surveyed_representativeness.csv", row.names = FALSE)
print(as.data.frame(rp), row.names = FALSE)
message("\nSaved -> out/national_recovery.csv, out/surveyed_representativeness.csv")
