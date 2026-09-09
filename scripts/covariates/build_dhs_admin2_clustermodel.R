# =============================================================================
# scripts/covariates/build_dhs_admin2_clustermodel.R   [DS-01]
#
# DHS-DERIVED ADMIN-2 PREDICTORS THE DHS WAY: A CLUSTER-LEVEL SPATIAL MODEL,
# NO COVARIATES, EVERY DISTRICT.
#
# The shared predictor set's 145 DHS columns were direct estimates: a
# survey-weighted mean of the clusters that fall inside each district (77
# custom and 51 extra indicators), with an area-level Fay-Herriot smooth over
# them for the 18 StatCompiler-coded indicators. A district with no DHS
# cluster had no value, a district with one cluster had that cluster. This
# script re-estimates every BINARY indicator with surveyPrev's cluster-level
# model (clusterModel: binomial on the cluster counts, BYM2 area effects over
# the GADM Admin-2 adjacency, PC priors; the unit-level model the surveyPrev /
# SUMMER authors recommend for Admin-2 and a covariate-free cousin of the DHS
# Program's model-based geostatistics). Every district gets a posterior mean,
# including districts with no cluster. Continuous indicators (means, counts)
# are left to the direct estimate, since the model is binomial.
#
# Inputs, all on disk (no DHS API call):
#   rdhs cache recodes (IR / KR / PR / HR / BR) for the round nearest each
#     micronutrient survey, the same rounds the shared-set builder uses
#   data/DHS/clean/<Country>_<year>_cluster_admin_info.rds: surveyPrev
#     cluster.info (cluster GPS in GADM polygons) and admin.info2 (adjacency)
#   the indicator definitions: surveyPrev built-ins (getDHSindicator), the 77
#     custom derivations (functions parsed out of
#     src/DHS/DHS_custom_admin2_indicators.R, whose main loop is not run) and
#     the extra derivers (src/DHS/DHS_extra_indicators_2026-09.R)
# Only indicators that exist as dhs_<name> columns in the current shared set
# are fitted, so the leakage exclusions carry over unchanged.
#
#   DHS_COUNTRY=Gambia Rscript -e "source('scripts/covariates/build_dhs_admin2_clustermodel.R')"
#   (one country per process; run the four in parallel)
# -> data/covariates/harmonized/dhs_admin2_clustermodel_<Country>.csv
#    data/covariates/harmonized/dhs_admin2_clustermodel_<Country>_log.csv
# The shared-set builder (scripts/covariates/build_shared_predictor_set.R)
# overrides its direct DHS estimates with these wherever they exist
# (DHS_ADMIN2_MODEL=direct restores the direct estimates).
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tibble); library(surveyPrev); library(INLA)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/survey_years.R")

CN <- Sys.getenv("DHS_COUNTRY", "Gambia")
CACHE <- "C:/Users/andre/AppData/Local/andre/rdhs/Cache/datasets"
HDIR  <- "data/covariates/harmonized"
DHS_SETS <- list(   # same rounds and GPS files as build_shared_predictor_set.R
  Gambia      = list(IR = "GMIR81DT", KR = "GMKR81DT", PR = "GMPR81DT", HR = "GMHR81DT", BR = "GMBR81DT", year = 2019, clean = "Gambia_2019"),
  Ghana       = list(IR = "GHIR72DT", KR = "GHKR72DT", PR = "GHPR72DT", HR = "GHHR72DT", BR = "GHBR72DT", year = 2014, clean = "Ghana_2014"),
  Malawi      = list(IR = "MWIR7ADT", KR = "MWKR7ADT", PR = "MWPR7ADT", HR = "MWHR7ADT", BR = "MWBR7ADT", year = 2015, clean = "Malawi_2015"),
  SierraLeone = list(IR = "SLIR61DT", KR = "SLKR61DT", PR = "SLPR61DT", HR = "SLHR61DT", BR = "SLBR61DT", year = 2013, clean = "Sierra Leone_2013"))
sp <- DHS_SETS[[CN]]; if (is.null(sp)) stop("unknown DHS_COUNTRY: ", CN)
BUILTINS <- c("AN_NUTS_W_THN", "CH_DIAT_C_ORT", "CH_VACC_C_BAS", "CH_VACC_C_DP1", "CH_VACC_C_DP3", "CH_VACC_C_MSL", "CH_VACC_C_NON",
              "CM_ECMR_C_NNR", "CN_BRFS_C_EXB", "CN_NUTS_C_HA2", "CN_NUTS_C_WH2", "FP_CUSA_W_MOD", "FP_NADA_W_UNT", "ML_NETP_H_IT2",
              "RH_ANCN_W_N4P", "RH_DELA_C_SKP", "WS_TLET_H_IMP", "WS_TLET_P_BAS")

# ── the derivation functions, without running either script's main loop ──────
ex <- parse("src/DHS/DHS_custom_admin2_indicators.R"); n_fun <- 0L
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) { eval(e, globalenv()); n_fun <- n_fun + 1L }
source("src/DHS/DHS_extra_indicators_2026-09.R")     # defines EXTRA_DERIVERS (list by recode) and derive_* functions
cat(sprintf("[DS-01] %s: %d custom derivation functions, extra derivers for %s\n", CN, n_fun, paste(names(EXTRA_DERIVERS), collapse = "/")))

# ── inputs ───────────────────────────────────────────────────────────────────
rd <- function(x) { p <- file.path(CACHE, paste0(x, ".rds")); if (file.exists(p)) readRDS(p) else NULL }
t0 <- Sys.time()
dl <- list(IRdata = rd(sp$IR), KRdata = rd(sp$KR), PRdata = rd(sp$PR), HRdata = rd(sp$HR), BRdata = rd(sp$BR))
cat(sprintf("[DS-01] recodes: %s (%.0f s)\n", paste(sprintf("%s %d", names(dl), sapply(dl, function(d) if (is.null(d)) 0L else nrow(d))), collapse = ", "), as.numeric(Sys.time() - t0, units = "secs")))
info <- readRDS(file.path("data/DHS/clean", paste0(sp$clean, "_cluster_admin_info.rds")))
cluster.info <- info$cluster.info; admin.info2 <- info$admin.info2
cat(sprintf("[DS-01] spatial: %d clusters with GPS, %d Admin-2 areas\n", nrow(cluster.info$data), nrow(admin.info2$data)))
# guard: recode clusters must be the GPS file's clusters (numbering repeats across rounds)
irc <- unique(as_num(dl$IRdata$v001)); ov <- mean(irc %in% cluster.info$data$cluster)
cat(sprintf("[DS-01] recode/GPS cluster overlap %.0f%%\n", 100 * ov)); if (ov < 0.95) stop("cluster overlap below 95%: wrong round")
SHARED <- names(read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), nrows = 2, check.names = FALSE))
want <- sub("^dhs_", "", grep("^dhs_", SHARED, value = TRUE))

# ── the indicator table: name -> surveyPrev-format data (cluster, householdID, weight, value) ──
IND <- list(); origin <- character()
for (id in intersect(BUILTINS, want)) {
  got <- NULL
  for (rc in c("IRdata", "KRdata", "HRdata", "PRdata", "BRdata")) {
    if (is.null(dl[[rc]])) next
    r <- tryCatch(suppressWarnings(suppressMessages(getDHSindicator(dl[[rc]], indicator = id))), error = function(e) NULL)
    if (!is.null(r) && "value" %in% names(r) && sum(!is.na(r$value)) > 0) { got <- r; break }
  }
  if (!is.null(got)) { IND[[id]] <- got; origin[id] <- "builtin" }
}
cr <- derive_all_indicators(dl)
for (nm in intersect(names(cr$standardized), want)) { IND[[nm]] <- cr$standardized[[nm]]; origin[nm] <- "custom" }
for (rc in c("IR", "KR", "HR")) {
  d <- dl[[paste0(rc, "data")]]; if (is.null(d)) next
  add <- character()
  for (fn in EXTRA_DERIVERS[[rc]]) { r <- fn(d); d <- r$df; add <- c(add, r$added) }
  for (v in intersect(add, want)) { s <- standardize_indicator(d, rc, v); if (!is.null(s) && nrow(s)) { IND[[v]] <- s; origin[v] <- "extra" } }
}
isbin <- vapply(IND, function(d) all(stats::na.omit(d$value) %in% c(0, 1)), logical(1))
cat(sprintf("[DS-01] indicators matched to the shared set: %d (builtin %d, custom %d, extra %d); binary %d, continuous %d (left to the direct estimate)\n",
            length(IND), sum(origin == "builtin"), sum(origin == "custom"), sum(origin == "extra"), sum(isbin), sum(!isbin)))

# ── fits ─────────────────────────────────────────────────────────────────────
areas <- admin.info2$data[, c("admin2.name.full", "admin1.name", "admin2.name")]
out <- data.frame(country = CN, Admin1 = areas$admin1.name, Admin2 = areas$admin2.name, stringsAsFactors = FALSE)
LOG <- list(); i <- 0L
for (nm in names(IND)) {
  i <- i + 1L
  d <- IND[[nm]]
  rec <- list(country = CN, column = paste0("dhs_", nm), origin = origin[[nm]], binary = isbin[[nm]], n_obs = sum(is.finite(d$value)),
              n_clusters = length(unique(d$cluster[is.finite(d$value)])), method = NA_character_, areas = NA_integer_, secs = NA_real_, message = "")
  if (!isbin[[nm]]) { rec$method <- "direct (continuous)"; LOG[[nm]] <- rec; next }
  dd <- as.data.frame(d[, intersect(c("cluster", "householdID", "weight", "value"), names(d))])
  dd <- dd[stats::complete.cases(dd), ]
  t1 <- Sys.time()
  fit <- tryCatch(suppressWarnings(suppressMessages(clusterModel(data = dd, cluster.info = cluster.info, admin.info = admin.info2, admin = 2, model = "bym2", aggregation = FALSE))),
                  error = function(e) { rec$message <<- conditionMessage(e); NULL })
  rec$secs <- round(as.numeric(Sys.time() - t1, units = "secs"), 1)
  if (is.null(fit)) { rec$method <- "FAILED"; LOG[[nm]] <- rec
    cat(sprintf("  [%3d/%d] %-32s FAILED after %5.1f s: %s\n", i, length(IND), nm, rec$secs, substr(rec$message, 1, 80))); flush.console(); next }
  est <- fit$res.admin2$mean[match(areas$admin2.name.full, fit$res.admin2$admin2.name.full)]
  out[[paste0("dhs_", nm)]] <- as.numeric(est)
  rec$method <- "clusterModel bym2"; rec$areas <- sum(is.finite(est)); LOG[[nm]] <- rec
  cat(sprintf("  [%3d/%d] %-32s %3d areas  %5.1f s\n", i, length(IND), nm, rec$areas, rec$secs)); flush.console()
  if (i %% 10 == 0) write.csv(out, file.path(HDIR, paste0("dhs_admin2_clustermodel_", CN, "_partial.csv")), row.names = FALSE)
}
write.csv(out, file.path(HDIR, paste0("dhs_admin2_clustermodel_", CN, ".csv")), row.names = FALSE)
write.csv(bind_rows(lapply(LOG, as.data.frame)), file.path(HDIR, paste0("dhs_admin2_clustermodel_", CN, "_log.csv")), row.names = FALSE)
unlink(file.path(HDIR, paste0("dhs_admin2_clustermodel_", CN, "_partial.csv")))
cat(sprintf("[DS-01] %s done: %d columns modelled over %d areas in %.1f min\n", CN, ncol(out) - 3L, nrow(out), as.numeric(Sys.time() - t0, units = "mins")))
