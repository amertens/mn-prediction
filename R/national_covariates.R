# =============================================================================
# R/national_covariates.R
#
# National (country x year) covariates for the VMNIS prevalence model, and the
# fold-safe machinery it needs: nearest-value carry, leakage-safe prescreening,
# missingness indicators, and kNN imputation with training-only donors.
#
# Promoted 2026-08-31 from sandbox_parsimony/R/37_national_covariates.R, whose
# only changes are the cache paths, which now live under data/national/ so the
# targets pipeline can depend on them.
#
# THREE SOURCES. "wdi" is the 17 a-priori World Bank indicators chosen on the
# same causal story as the curated Admin-2 set. "panel" is the project's own
# country x year panel (~1,400 usable columns, essentially the DHS set).
# "both" is the union, restricted to the countries both cover.
#
# ANAEMIA COLUMNS ARE DROPPED. Predicting iron deficiency from measured anaemia
# is close to predicting a deficiency from a deficiency; supplementation and
# dietary columns are kept, because those are exposures and interventions.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

PANEL_DTA <- Sys.getenv(
  "NAT_PANEL_DTA",
  "C:/Users/andre/OneDrive/Documents/mn-proxies/data/national/combined_dataset.dta")
NAT_DIR     <- here::here("data", "national")
PANEL_CACHE <- file.path(NAT_DIR, "panel_national.rds")
WDI_CACHE   <- file.path(NAT_DIR, "wdi_national.rds")

WDI_VARS <- c(
  gdp_pc        = "NY.GDP.PCAP.PP.KD", undernourish  = "SN.ITK.DEFC.ZS",
  food_prod_idx = "AG.PRD.FOOD.XD",    cereal_yield  = "AG.YLD.CREL.KG",
  ag_land_pct   = "AG.LND.AGRI.ZS",    urban_pct     = "SP.URB.TOTL.IN.ZS",
  electricity   = "EG.ELC.ACCS.ZS",    water_basic   = "SH.H2O.BASW.ZS",
  sanitation    = "SH.STA.BASS.ZS",    u5_mortality  = "SH.DYN.MORT",
  life_exp      = "SP.DYN.LE00.IN",    dtp3          = "SH.IMM.IDPT",
  health_exp_pc = "SH.XPD.CHEX.PC.CD", fertility     = "SP.DYN.TFRT.IN",
  female_school = "SE.SEC.ENRR.FE",    pop_total     = "SP.POP.TOTL",
  rural_pop_pct = "SP.RUR.TOTL.ZS")

WDI_LOG_VARS <- c("gdp_pc", "pop_total", "health_exp_pc", "cereal_yield")

# ---------------------------------------------------------------------------
# Nearest-value carry within a country, matrix-wise.
#
# The per-column vapply version in 31/33 is fine for 17 indicators and far too
# slow for 1,400. This walks rows instead of columns, so every column is filled
# in one vectorised pass, and tracks the distance carried so `max_gap` is
# enforced exactly as before: last observation carried forward, then next
# observation carried backward for whatever is still missing.
# ---------------------------------------------------------------------------
.carry <- function(M, max_gap, from_last = FALSE) {
  n <- nrow(M)
  if (n < 2L) return(M)
  idx <- if (from_last) seq.int(n - 1L, 1L) else seq.int(2L, n)
  prev_step <- if (from_last) 1L else -1L
  gap <- matrix(0L, n, ncol(M))
  for (i in idx) {
    j  <- i + prev_step
    na <- is.na(M[i, ])
    if (!any(na)) next
    M[i, na]   <- M[j, na]
    gap[i, na] <- gap[j, na] + 1L
  }
  M[gap > max_gap] <- NA_real_
  M
}

fill_near_matrix <- function(M, group, max_gap = 5L) {
  for (g in unique(group)) {
    r <- which(group == g)
    if (length(r) < 2L) next
    sub <- M[r, , drop = FALSE]
    f   <- .carry(sub, max_gap, from_last = FALSE)
    b   <- .carry(sub, max_gap, from_last = TRUE)
    f[is.na(f)] <- b[is.na(f)]
    M[r, ] <- f
  }
  M
}

# ---------------------------------------------------------------------------
# The 17 a-priori WDI indicators (the original source, kept as the control)
# ---------------------------------------------------------------------------
load_wdi_covariates <- function(max_gap = 5L) {
  if (!file.exists(WDI_CACHE)) {
    stopifnot(requireNamespace("WDI", quietly = TRUE))
    message("downloading ", length(WDI_VARS), " WDI indicators ...")
    w <- WDI::WDI(indicator = WDI_VARS, start = 1985, end = 2023, extra = TRUE)
    saveRDS(w, WDI_CACHE)
  }
  wdi <- readRDS(WDI_CACHE) |>
    filter(!is.na(iso3c), region != "Aggregates") |>
    mutate(iso3c = as.character(iso3c)) |>
    select(iso3c, year, all_of(names(WDI_VARS))) |>
    arrange(iso3c, year)

  M <- fill_near_matrix(as.matrix(wdi[, names(WDI_VARS)]), wdi$iso3c, max_gap)
  wdi[, names(WDI_VARS)] <- as.data.frame(M)

  list(df = wdi, vars = names(WDI_VARS), log_vars = WDI_LOG_VARS,
       source = "wdi")
}

# ---------------------------------------------------------------------------
# The project's national panel
# ---------------------------------------------------------------------------
build_panel_covariates <- function(min_coverage = 0.5,
                                   yr_range     = c(1990L, 2022L),
                                   max_gap      = 5L,
                                   drop_anaemia = TRUE) {
  stopifnot(requireNamespace("haven", quietly = TRUE))
  if (!file.exists(PANEL_DTA))
    stop("national panel not found at: ", PANEL_DTA)

  message("reading national panel (this takes a minute) ...")
  raw <- haven::read_dta(PANEL_DTA)
  message(sprintf("  raw: %d rows x %d columns", nrow(raw), ncol(raw)))

  # Column names are STATcompiler codes (AN_ANEM_C_ANY, ML_HEMO_C_HL8); the
  # semantics are in the Stata variable labels, so anything that filters on
  # meaning has to read the labels.
  labs <- vapply(raw, function(x) {
    l <- attr(x, "label"); if (is.null(l)) "" else as.character(l)[1]
  }, character(1))

  raw <- raw |>
    mutate(iso3c = suppressWarnings(countrycode::countrycode(
             as.character(country_name), "country.name", "iso3c", warn = FALSE)),
           year = as.integer(year)) |>
    filter(!is.na(iso3c), !is.na(year)) |>
    arrange(iso3c, year)

  num  <- vapply(raw, is.numeric, logical(1))
  vars <- setdiff(names(raw)[num], c("year", "population"))

  n_anaemia <- 0L
  if (drop_anaemia) {
    # Measured anaemia/haemoglobin STATUS only. Supplementation and dietary
    # columns ("received iron tablets", "vitamin A-rich foods", "micronutrient
    # powder") mention the nutrients but are interventions and exposures, which
    # are legitimate -- and mechanistically among the most relevant -- proxies.
    hit <- vars[grepl("an(a)?emi|h(a)?emoglobin", labs[vars], ignore.case = TRUE)]
    n_anaemia <- length(hit)
    vars <- setdiff(vars, hit)
  }

  M <- as.matrix(raw[, vars])
  storage.mode(M) <- "double"
  M <- fill_near_matrix(M, raw$iso3c, max_gap)

  in_win <- raw$year >= yr_range[1] & raw$year <= yr_range[2]
  cover  <- colMeans(!is.na(M[in_win, , drop = FALSE]))
  nuniq  <- apply(M[in_win, , drop = FALSE], 2,
                  function(x) length(unique(x[is.finite(x)])))
  keep   <- cover >= min_coverage & nuniq > 1L

  message(sprintf(
    "  numeric %d | anaemia/Hb status dropped %d | >=%.0f%% covered after +/-%dy carry: %d",
    sum(num), n_anaemia, 100 * min_coverage, max_gap, sum(keep)))

  kept <- vars[keep]
  out <- data.frame(iso3c = raw$iso3c, year = raw$year,
                    as.data.frame(M[, keep, drop = FALSE]),
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[-(1:2)] <- make.names(names(out)[-(1:2)], unique = TRUE)

  # Keep the labels: the codes are unreadable, and the importance table is the
  # point of running this at all.
  attr(out, "labels")            <- setNames(unname(labs[kept]), names(out)[-(1:2)])
  attr(out, "n_anaemia_dropped") <- n_anaemia
  attr(out, "n_raw_columns")     <- sum(num)
  out
}

load_panel_covariates <- function(min_coverage = 0.5, refresh = FALSE, ...) {
  key <- sprintf("%s_cov%02d.rds", sub("\\.rds$", "", PANEL_CACHE),
                 round(100 * min_coverage))
  if (!refresh && file.exists(key)) {
    p <- readRDS(key)
    message(sprintf("using cached panel (%d covariates)", ncol(p) - 2L))
  } else {
    p <- build_panel_covariates(min_coverage = min_coverage, ...)
    saveRDS(p, key)
  }
  v <- setdiff(names(p), c("iso3c", "year"))
  # right-skewed money/population/count columns; log1p only where non-negative
  lg <- v[vapply(p[v], function(x) {
    x <- x[is.finite(x)]
    length(x) > 0 && min(x) >= 0 && stats::sd(x) > 0 &&
      (max(x) / (stats::median(x[x > 0]) + 1)) > 50
  }, logical(1))]
  list(df = p, vars = v, log_vars = lg, source = "panel",
       labels = attr(p, "labels"))
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
load_national_covariates <- function(source = Sys.getenv("NAT_COV", "wdi"),
                                     min_coverage = as.numeric(
                                       Sys.getenv("NAT_COV_MIN", "0.5"))) {
  source <- match.arg(source, c("wdi", "panel", "both"))
  if (source == "wdi")   return(load_wdi_covariates())
  if (source == "panel") return(load_panel_covariates(min_coverage))

  w <- load_wdi_covariates(); p <- load_panel_covariates(min_coverage)
  clash <- intersect(w$vars, p$vars)
  if (length(clash)) names(p$df)[match(clash, names(p$df))] <-
    paste0("pnl_", clash)
  p$vars     <- ifelse(p$vars %in% clash, paste0("pnl_", p$vars), p$vars)
  p$log_vars <- ifelse(p$log_vars %in% clash, paste0("pnl_", p$log_vars),
                       p$log_vars)
  df <- dplyr::inner_join(w$df, p$df, by = c("iso3c", "year"))
  list(df = df, vars = c(w$vars, p$vars),
       log_vars = c(w$log_vars, p$log_vars), source = "both")
}

# ---------------------------------------------------------------------------
# Leakage-safe prescreen.
#
# With ~50 held-out folds and up to ~1,400 candidate columns, the screen has to
# happen INSIDE the fold or it sees the held-out country. Ranked by absolute
# Spearman correlation with the training outcome -- rank-based so it is not
# driven by the heavy right skew several of these columns carry.
# ---------------------------------------------------------------------------
prescreen_cols <- function(Xtr, ytr, k) {
  if (ncol(Xtr) <= k) return(seq_len(ncol(Xtr)))
  s <- suppressWarnings(apply(Xtr, 2, function(x) {
    if (!any(is.finite(x)) || stats::sd(x) < 1e-12) return(0)
    abs(stats::cor(x, ytr, method = "spearman"))
  }))
  s[!is.finite(s)] <- 0
  sort(order(s, decreasing = TRUE)[seq_len(k)])
}

# ---------------------------------------------------------------------------
# Common support.
#
# The panel holds 57 countries (essentially the DHS set); WDI holds ~200. Left
# alone, the two sources are scored on different countries -- the panel's subset
# is smaller, more African and higher-burden -- so a raw wdi-vs-panel comparison
# confounds the covariates with the country sample. Restricting both to the
# intersection is what makes the contrast interpretable.
# ---------------------------------------------------------------------------
common_countries <- function(min_coverage = as.numeric(
                               Sys.getenv("NAT_COV_MIN", "0.5"))) {
  w <- load_wdi_covariates()$df
  p <- load_panel_covariates(min_coverage)$df
  sort(intersect(unique(w$iso3c), unique(p$iso3c)))
}

# ---------------------------------------------------------------------------
# Missingness indicators.
#
# For a DHS-derived panel, "this country never reported this indicator" is
# itself informative -- reporting is patterned by survey programme, era and
# capacity, all of which correlate with nutritional status. Median-imputing
# throws that signal away and quietly pretends the value was average. An
# explicit 0/1 column lets the model use it instead, and lets the importance
# table show whether it did.
#
# Rates are computed on TRAINING rows only, so which indicators get one is not
# influenced by the held-out country. The indicator VALUES for the held-out
# rows are observable metadata (was this reported or not), never its outcome.
# ---------------------------------------------------------------------------
missing_indicators <- function(X, tr, lo = 0.02, hi = 0.98) {
  rate <- colMeans(is.na(X[tr, , drop = FALSE]))
  use  <- which(rate > lo & rate < hi)
  if (!length(use)) return(NULL)
  M <- matrix(as.numeric(is.na(X[, use, drop = FALSE])), nrow = nrow(X))
  colnames(M) <- paste0("mi_", colnames(X)[use])
  # an indicator that is constant across the training rows carries nothing
  keep <- apply(M[tr, , drop = FALSE], 2, function(x) length(unique(x)) > 1L)
  if (!any(keep)) return(NULL)
  M[, keep, drop = FALSE]
}

# ---------------------------------------------------------------------------
# k-nearest-neighbour imputation, fold-safe.
#
# Donors are TRAINING rows only, and the standardisation that defines the
# distance is fitted on training rows only, so nothing about the held-out
# country informs how it or anyone else is filled.
#
# Distance is mean squared difference over the dimensions two rows both
# observe, which is the usable definition when the missingness is this patchy
# (the median column is a third populated). Rows that overlap a donor on fewer
# than `min_overlap` dimensions cannot be compared meaningfully, so they fall
# back to the training mean rather than borrowing from an arbitrary neighbour.
# ---------------------------------------------------------------------------
knn_impute_fold <- function(X, tr, k = 5L, min_overlap = 5L) {
  mu <- colMeans(X[tr, , drop = FALSE], na.rm = TRUE)
  sdv <- apply(X[tr, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0; sdv[!is.finite(sdv) | sdv < 1e-9] <- 1

  Z  <- sweep(sweep(X, 2, mu, "-"), 2, sdv, "/")
  D  <- Z[tr, , drop = FALSE]                     # donor pool
  obsD <- !is.na(D)
  k <- min(k, nrow(D))

  for (i in seq_len(nrow(Z))) {
    na_i <- which(is.na(Z[i, ]))
    if (!length(na_i)) next
    obs_i <- !is.na(Z[i, ])
    # exclude self when the row is itself a donor
    donors <- if (i %in% tr) which(tr != i) else seq_len(nrow(D))
    if (!length(donors)) { Z[i, na_i] <- 0; next }

    Dd <- D[donors, , drop = FALSE]; od <- obsD[donors, , drop = FALSE]
    both <- sweep(od, 2, obs_i, "&")
    n_ov <- rowSums(both)
    diff <- sweep(Dd, 2, Z[i, ], "-")
    diff[!both] <- 0
    dist <- ifelse(n_ov >= min_overlap, rowSums(diff^2) / pmax(n_ov, 1), Inf)

    if (all(!is.finite(dist))) { Z[i, na_i] <- 0; next }
    nn <- donors[order(dist)[seq_len(min(k, sum(is.finite(dist))))]]
    fill <- colMeans(D[nn, na_i, drop = FALSE], na.rm = TRUE)
    fill[!is.finite(fill)] <- 0        # no neighbour observed it either
    Z[i, na_i] <- fill
  }
  sweep(sweep(Z, 2, sdv, "*"), 2, mu, "+")
}
