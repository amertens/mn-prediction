# =============================================================================
# dashboard/data-raw/05_build_protocol_v2_bundles.R
#
# The dashboard's data, rebuilt from the corrected protocol (protocol v2,
# RR-10 tables of 2026-09-10). Until this script existed every bundle in
# dashboard/data/ came from the pre-audit pipeline (the person-level
# SuperLearner, the area-level recipe, the old leaderboard), so the app
# contradicted the decks and the manuscript. Everything the app now shows is
# built here, from results tables that are already committed, plus one fit:
# the deployment ranking, which no result table holds.
#
#   Rscript dashboard/data-raw/05_build_protocol_v2_bundles.R
#
# Writes to dashboard/data/ (gitignored; rebuilt before each deploy):
#   admin2_index.rds        the deployment ranking for every district of the
#                           four surveyed countries: the zero-tuning index fitted
#                           on all surveyed districts of a country and applied to
#                           every polygon, with the survey's own estimate where
#                           one exists, the out-of-fold worst-fifth probability
#                           (script scripts/policy_deck/07), a national anchor
#                           turning the ranking into a planning prevalence, and
#                           the exact per-column weights so a district's score
#                           can be decomposed into predictor contributions.
#   civ_index.rds           Cote d'Ivoire ranked from climate and soil alone
#                           (scripts/policy_deck/04-06), with rank uncertainty.
#   protocol_evidence.rds   the benchmark, targeting, ceiling, survey-design,
#                           importance and comparator tables the trust tabs read.
#   predictor_catalogue.rds every predictor with its metadata, definition,
#                           index weight per outcome, cross-country sign
#                           replication, and the within-country rank-normal
#                           values for the small maps.
#   metadata.rds            build stamp refreshed.
#
# Country spellings: the protocol tables write "SierraLeone", the population
# and boundary bundles "Sierra Leone". Everything written here uses the
# dashboard's labels ("Sierra Leone") and keys ("sierraleone").
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
setwd(here::here())
source("R/protocol_v2.R")
source("R/protocol_v2_weights.R")
source("R/protocol_v2_importance.R")
source("R/admin2_key_hygiene.R")

P2  <- "results/tables/protocol_v2"
PDK <- "results/tables/policy_deck"
CLU <- "results/tables/cluster_level"
HD  <- "data/covariates/harmonized"
OUT <- "dashboard/data"
BUILD_TIME <- Sys.time()

rd <- function(...) {
  f <- file.path(...)
  if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}
LABEL <- c(Gambia = "Gambia", Ghana = "Ghana", Malawi = "Malawi",
           SierraLeone = "Sierra Leone", `Sierra Leone` = "Sierra Leone")
KEY   <- c(Gambia = "gambia", Ghana = "ghana", Malawi = "malawi",
           SierraLeone = "sierraleone", `Sierra Leone` = "sierraleone")
relabel <- function(d, col = "country") {
  if (!is.null(d) && col %in% names(d)) d[[col]] <- unname(LABEL[d[[col]]])
  d
}
wilson <- function(p, n, z = qnorm(0.975)) {
  ok <- is.finite(p) & is.finite(n) & n > 0
  lo <- hi <- rep(NA_real_, length(p))
  den <- 1 + z^2 / n[ok]; cen <- (p[ok] + z^2 / (2 * n[ok])) / den
  half <- z * sqrt(p[ok] * (1 - p[ok]) / n[ok] + z^2 / (4 * n[ok]^2)) / den
  lo[ok] <- pmax(0, cen - half); hi[ok] <- pmin(1, cen + half)
  list(lo = lo, hi = hi)
}
meta <- readRDS(file.path(OUT, "metadata.rds"))
who_class_of <- function(prev, oc) {
  th <- meta$who_thresholds[[oc]]
  if (is.null(th)) return(rep(NA_character_, length(prev)))
  cut(prev, breaks = c(-Inf, sort(as.numeric(th)), Inf),
      labels = c("Low", "Mild", "Moderate", "Severe"), right = FALSE) |> as.character()
}

# ── A. The deployment ranking ───────────────────────────────────────────────
cat("A. deployment ranking\n")
TG  <- rd(P2, "targets_v2.csv")
S   <- rd(HD, "predictors_admin2_shared.csv")
MD  <- rd(HD, "predictors_admin2_shared_metadata.csv")
WF  <- rd(PDK, "worst_fifth_probability.csv")
NE  <- rd("results/tables/national_estimates_all.csv")
POP <- readRDS(file.path(OUT, "admin2_population.rds"))
PREDS <- intersect(MD$column, names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
stopifnot(length(PREDS) == nrow(MD))
S <- S[!is_water_admin2(S$Admin2), ]        # GADM ships Lake Malawi as districts
cat(sprintf("   %d predictors, %d districts after dropping water polygons\n", length(PREDS), nrow(S)))

# The rank-normal design matrix per country (what the model sees), kept for the
# catalogue's small maps and for the district decomposition.
XR <- list(); MU <- list()
for (ctry in unique(S$country)) {
  all_s <- S[S$country == ctry, ]
  Xr <- prep_predictors_v2(as.matrix(all_s[, PREDS]))
  rownames(Xr) <- paste(all_s$Admin1, all_s$Admin2, sep = "|")
  XR[[KEY[[ctry]]]] <- Xr
}

districts <- list(); fits <- list(); national <- list()
for (ctry in unique(TG$country)) {
  all_s <- S[S$country == ctry, ]
  Xr_all <- XR[[KEY[[ctry]]]]
  key_all <- paste(all_s$Admin1, all_s$Admin2, sep = "|")
  pop <- POP[POP$country == LABEL[[ctry]], ]
  pop_key <- paste(pop$Admin1, pop$Admin2, sep = "|")
  for (oc in unique(TG$outcome[TG$country == ctry])) {
    t <- TG[TG$country == ctry & TG$outcome == oc & is.finite(TG$y_prev) & is.finite(TG$n_eff), ]
    key_t <- paste(t$Admin1, t$Admin2, sep = "|")
    tr <- match(key_t, key_all); keep <- is.finite(tr); tr <- tr[keep]; t <- t[keep, ]
    if (length(tr) < 8) { cat(sprintf("   skip %s / %s: %d surveyed districts\n", ctry, oc, length(tr))); next }
    Y <- rep(NA_real_, nrow(all_s)); Y[tr] <- .v2_logit(t$y_prev)
    D <- domain_representation_v2(Xr_all, domain_of, sign_rows = tr)
    pred <- ARMS_V2[["domain_index"]](tr, seq_len(nrow(all_s)), Y, NULL, D,
                                      list(Admin1 = all_s$Admin1, y_nat = Y))
    # exact back-projection of the fitted weights onto the columns, on the
    # logit scale the prediction is reported on
    w <- .ws_z_pooled(tr, Y, D)
    beta <- index_backproject_v2(w, attr(D, "basis"), colnames(Xr_all))
    idx_tr <- as.numeric(D[tr, , drop = FALSE] %*% w)
    scale <- if (stats::sd(idx_tr) > 0) stats::sd(Y[tr]) / stats::sd(idx_tr) else 0
    mu <- colMeans(Xr_all[tr, , drop = FALSE])
    fits[[paste(KEY[[ctry]], oc)]] <- list(country_key = KEY[[ctry]], outcome = oc,
                                            beta = beta * scale, mu = mu,
                                            intercept = mean(Y[tr]) - sum(beta * scale * mu),
                                            n_train = length(tr))
    n <- nrow(all_s)
    rk <- rank(-pred, ties.method = "average")
    d <- data.frame(country = LABEL[[ctry]], country_key = KEY[[ctry]], outcome = oc,
                    Admin1 = all_s$Admin1, Admin2 = all_s$Admin2,
                    score_logit = pred, rank_worst = rk, n_districts = n,
                    priority = 100 * (n - rk + 0.5) / n,
                    prev_model = .v2_expit(pred), stringsAsFactors = FALSE)
    # the survey's own estimate where the district was surveyed
    i <- match(key_all, key_t)
    d$surveyed   <- is.finite(i)
    d$survey_prev <- t$y_prev[i]; d$n_resp <- t$n_raw[i]; d$n_clusters <- t$n_psu[i]; d$n_eff <- t$n_eff_district[i]
    ci <- wilson(d$survey_prev, d$n_eff); d$survey_lo <- ci$lo; d$survey_hi <- ci$hi
    # out-of-fold probability of the worst fifth (surveyed districts, 40 draws)
    d$p_worst_fifth <- NA_real_
    if (!is.null(WF)) {
      wf <- WF[WF$country == ctry & WF$outcome == oc, ]
      j <- match(key_all, paste(wf$Admin1, wf$Admin2, sep = "|"))
      d$p_worst_fifth <- wf$p_worst_fifth[j]
    }
    # population of the group the biomarker describes
    pj <- match(key_all, pop_key)
    d$population <- if (startsWith(oc, "child_")) pop$pop_child[pj] else pop$pop_women[pj]
    # national anchor: shift the logit scores so the population-weighted mean
    # equals the survey's national prevalence (the AR-01 design, at 100 percent
    # of the national sample). This is the planning prevalence; the ranking is
    # the product.
    p_nat <- NE$obs_prev[NE$country == LABEL[[ctry]] & NE$outcome == oc][1]
    if (!is.finite(p_nat)) p_nat <- stats::weighted.mean(t$y_prev, t$n_raw)
    ok <- is.finite(d$population) & d$population > 0
    f <- function(c) sum(d$population[ok] * .v2_expit(pred[ok] + c)) / sum(d$population[ok]) - p_nat
    shift <- tryCatch(stats::uniroot(f, c(-12, 12))$root, error = function(e) NA_real_)
    d$prev_anchored <- if (is.finite(shift)) .v2_expit(pred + shift) else NA_real_
    d$who_class <- who_class_of(d$prev_anchored, oc)
    d$people_affected <- d$prev_anchored * d$population
    districts[[length(districts) + 1]] <- d
    national[[length(national) + 1]] <- data.frame(
      country = LABEL[[ctry]], country_key = KEY[[ctry]], outcome = oc, national_prev = p_nat,
      n_surveyed = length(tr), n_districts = n, anchor_shift = shift, stringsAsFactors = FALSE)
    cat(sprintf("   %-12s %-13s surveyed %3d of %3d, national %.3f, anchor shift %+.2f\n",
                LABEL[[ctry]], oc, length(tr), n, p_nat, shift))
  }
}
districts <- bind_rows(districts); national <- bind_rows(national)
dup <- duplicated(paste(districts$country, districts$outcome, districts$Admin1, districts$Admin2))
if (any(dup)) stop(sprintf("admin2_index: %d duplicated country/outcome/Admin1/Admin2 keys", sum(dup)))
saveRDS(list(districts = districts, national = national, fits = fits, xr = XR,
             build_time = BUILD_TIME, protocol = "v2, RR-10 set (2026-09-10)"),
        file.path(OUT, "admin2_index.rds"))
cat(sprintf("   wrote admin2_index.rds: %d district rows, %d cells\n", nrow(districts), length(fits)))

# ── B. Cote d'Ivoire ────────────────────────────────────────────────────────
cat("B. Cote d'Ivoire\n")
civ <- rd(PDK, "civ_climate_soil_ranking.csv")
civ_u <- rd(PDK, "civ_rank_uncertainty.csv")
oos_old <- readRDS(file.path(OUT, "oos_cote_divoire.rds"))
if (!is.null(civ)) {
  n <- length(unique(paste(civ$Admin1, civ$Admin2)))
  civ$priority <- 100 - civ$pct
  civ$rank_worst <- civ$rank
  saveRDS(list(ranking = civ, uncertainty = civ_u, guards = rd(PDK, "civ_transport_guards.csv"),
               boundaries = oos_old$boundaries, n_districts = n, build_time = BUILD_TIME),
          file.path(OUT, "civ_index.rds"))
  cat(sprintf("   wrote civ_index.rds: %d districts, %d outcomes\n", n, length(unique(civ$outcome))))
}

# ── C. The evidence behind the trust tabs ───────────────────────────────────
cat("C. protocol evidence\n")
ev <- list(
  benchmarks_summary = relabel(rd(P2, "benchmarks_v2_summary.csv")),
  benchmarks_cells   = relabel(rd(P2, "benchmarks_v2_cells.csv")),
  transport_null     = rd(P2, "transport_null_calibration.csv"),
  admin1_transport   = relabel(rd(P2, "admin1_transport.csv")),
  climate_soil_admin1 = relabel(rd(P2, "climate_soil_admin1.csv")),
  nested_domains     = relabel(rd(P2, "nested_domain_selection.csv"), "heldout"),
  training_curve     = relabel(rd(P2, "training_country_curve.csv"), "heldout"),
  training_curve_cs  = relabel(rd(P2, "training_curve_climate_soil.csv"), "heldout"),
  targeting_summary  = rd(P2, "nce_targeting_summary.csv"),
  targeting_cells    = relabel(rd(P2, "nce_targeting_metrics.csv")),
  worst_fifth_calibration = rd(PDK, "worst_fifth_calibration.csv"),
  risk_summary       = rd(P2, "risk_category_accuracy_summary.csv"),
  risk_cells         = relabel(rd(P2, "risk_category_accuracy.csv")),
  ceiling            = relabel(rd(P2, "variance_components_ceiling.csv")),
  design_summary     = rd(P2, "anchor_and_rank_summary.csv"),
  importance_top     = rd(P2, "index_importance_top.csv"),
  importance_domains = rd(P2, "index_importance_domains.csv"),
  importance_patterns = rd(P2, "index_importance_patterns.csv"),
  domain_ablation    = rd(P2, "domain_ablation_loco_summary.csv"),
  source_ablation    = rd(P2, "source_ablation_loco_summary.csv"),
  weight_sources     = rd(P2, "weight_sources_summary.csv"),
  geostat_cells      = relabel(rd(CLU, "mbg_comparison_cells.csv")),
  individual_level   = relabel(rd(P2, "individual_level_models.csv")),
  build_time = BUILD_TIME
)
missing <- names(ev)[vapply(ev, is.null, logical(1))]
if (length(missing)) cat("   MISSING tables:", paste(missing, collapse = ", "), "\n")
saveRDS(ev, file.path(OUT, "protocol_evidence.rds"))
cat(sprintf("   wrote protocol_evidence.rds: %d tables\n", sum(!vapply(ev, is.null, logical(1))) - 1))

# ── D. The predictor catalogue ──────────────────────────────────────────────
cat("D. predictor catalogue\n")
VS <- rd(P2, "variable_sheet.csv")
IC <- rd(P2, "index_importance_columns.csv")
IT <- rd(P2, "index_importance_top.csv")
ID <- rd(P2, "index_importance_domains.csv")
P4 <- rd("results/tables/signal_probes", "p4_admin1_continuous_predictors.csv")
DA <- rd(P2, "domain_ablation_loco_summary.csv")

source_label <- function(s) dplyr::case_when(
  grepl("^DHS", s) ~ "DHS / MICS household surveys", grepl("AlphaEarth", s) ~ "AlphaEarth satellite embedding",
  grepl("SoilGrids", s) ~ "SoilGrids / iSDA soil", grepl("^GEE", s) ~ "Earth Engine (climate, land, built environment)",
  grepl("IHME", s) ~ "IHME modelled surfaces", grepl("Malaria Atlas", s) ~ "Malaria Atlas Project",
  grepl("MapSPAM", s) ~ "MapSPAM crops", grepl("Koppen", s) ~ "Koppen / agro-ecological zones",
  grepl("WFP", s) ~ "WFP market prices", grepl("FAOSTAT", s) ~ "FAOSTAT (national)",
  grepl("Livestock", s) ~ "Gridded Livestock of the World", grepl("ESPEN", s) ~ "WHO ESPEN helminths",
  grepl("Surface Water", s) ~ "Earth Engine (water and coast distance)", grepl("GFDx", s) ~ "GFDx fortification",
  grepl("HFID", s) ~ "HFID food prices", TRUE ~ "WorldPop / GHS / RWI")
temporal_label <- c(annual_series = "annual series, matched to the survey year", static = "static layer",
                    survey_year = "matched to the survey year", fieldwork_matched = "matched to the fieldwork months")

vars <- MD |>
  transmute(column, domain, source, source_label = source_label(source), n_countries,
            countries = gsub("[|;]", ", ", gsub("SierraLeone", "Sierra Leone", countries)),
            completeness, subnational = as.logical(subnational))
if (!is.null(VS)) {
  vs <- VS |> transmute(column = variable, definition, unit, temporal_kind, var_type,
                        median = suppressWarnings(as.numeric(median)), value_range,
                        pct_missing = suppressWarnings(as.numeric(pct_missing_overall)),
                        coverage_note, source_note, subdomain = proposed_subdomain,
                        ra_priority = suppressWarnings(as.integer(ra_priority)), flags)
  vars <- left_join(vars, vs, by = "column")
}
vars$definition[!nzchar(trimws(vars$definition %||% ""))] <- NA_character_
vars$climate_soil <- vars$domain %in% c("Climate and weather", "Soil characteristics")
# membership of the twenty-predictor composite, per outcome (pooled fit, biomarker level)
if (!is.null(IT)) {
  comp <- IT |> filter(target == "level", rank <= 20) |> group_by(column) |>
    summarise(composite_outcomes = paste(sort(unique(outcome)), collapse = ", "), .groups = "drop")
  vars <- left_join(vars, comp, by = "column")
}
# index weight per outcome (pooled four-country fit), both targets
weights <- if (!is.null(IC)) IC |> filter(scope == "pooled") |>
  transmute(column, outcome, target, beta_std, share, rank, fold_sign_agree, loco_sign_agree) else NULL
# cross-country sign replication of the district association (signal probe P4)
signal <- if (!is.null(P4)) P4 |> filter(group %in% unique(IT$outcome)) |>
  transmute(column = predictor, outcome = group, meta_z, mu, k_countries, sign_agree, q_bh_perm) else NULL
# a summary line per variable: how many outcomes it enters in the top 20, and its best rank
if (!is.null(weights)) {
  ws <- weights |> filter(target == "level") |> group_by(column) |>
    summarise(best_rank = min(rank), n_top20 = sum(rank <= 20), mean_abs_weight = mean(abs(beta_std)), .groups = "drop")
  vars <- left_join(vars, ws, by = "column")
}
if (!is.null(signal)) {
  sg <- signal |> group_by(column) |>
    summarise(n_replicated = sum(sign_agree >= k_countries & k_countries >= 3), .groups = "drop")
  vars <- left_join(vars, sg, by = "column")
}
domains <- vars |> group_by(domain) |>
  summarise(n_columns = n(), sources = paste(sort(unique(source_label)), collapse = "; "),
            n_defined = sum(!is.na(definition)), .groups = "drop")
if (!is.null(ID)) {
  dsh <- ID |> filter(scope == "pooled", target == "level") |> group_by(domain) |>
    summarise(share_mean = mean(share), share_min = min(share), share_max = max(share),
              n_axes = round(mean(n_axes)), .groups = "drop")
  domains <- left_join(domains, dsh, by = "domain")
}
if (!is.null(DA)) domains <- left_join(domains, DA |> filter(target == "level") |>
                                         transmute(domain, transport_cost = delta_drop, transport_alone = only), by = "domain")
sources <- vars |> group_by(source_label) |>
  summarise(n_columns = n(), domains = paste(sort(unique(domain)), collapse = "; "),
            n_defined = sum(!is.na(definition)), .groups = "drop")
saveRDS(list(variables = vars, weights = weights, signal = signal, domains = domains, sources = sources,
             domain_shares = if (!is.null(ID)) ID |> filter(scope == "pooled", target == "level") |> select(domain, outcome, share, n_axes) else NULL,
             xr = XR, build_time = BUILD_TIME),
        file.path(OUT, "predictor_catalogue.rds"))
cat(sprintf("   wrote predictor_catalogue.rds: %d variables (%d with a definition), %d domains, %d sources\n",
            nrow(vars), sum(!is.na(vars$definition)), nrow(domains), nrow(sources)))

# ── E. build stamp ──────────────────────────────────────────────────────────
meta$build_timestamp <- BUILD_TIME
meta$protocol <- "Protocol v2, RR-10 result set (2026-09-10)"
saveRDS(meta, file.path(OUT, "metadata.rds"))
cat("done\n")
