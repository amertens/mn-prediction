# =============================================================================
# R/brinda_adjustment.R
#
# Uniform BRINDA-style inflammation adjustment for RBP (vitamin A), so every
# country's vitamin-A-deficiency outcome uses ONE consistent definition for the
# cross-country (LOCO) transport analysis — instead of the previous mix of
# Thurnham-adjusted (Gambia/Ghana), unspecified-"Adj" (Sierra Leone), and RAW
# (Malawi) RBP (DC-H2).
#
# Method (BRINDA regression correction; Larson et al. 2017):
#   1. Regress log(RBP) on log(CRP) + log(AGP) within the population group.
#   2. Reference = max of the lowest CRP/AGP decile (10th percentile).
#   3. adjusted log(RBP) = log(RBP)
#        - b_CRP * max(log(CRP) - log(CRP_ref), 0)
#        - b_AGP * max(log(AGP) - log(AGP_ref), 0)
#      Coefficients are CLAMPED to <= 0: inflammation can only depress RBP, so a
#      positive (collinearity-driven) coefficient is set to 0 rather than
#      inflating apparent deficiency.
#   4. VAD = adjusted RBP < 0.70 umol/L.
#
# Validation against published survey VAD prevalences is in
# docs/dc_h2_brinda_validation.md.
# =============================================================================

#' BRINDA-adjust RBP for inflammation (vitamin A).
#'
#' @param rbp,crp,agp numeric vectors (RBP in umol/L; CRP in mg/L; AGP in g/L),
#'   same length, one population group (e.g. preschool children OR women).
#' @param min_n minimum complete cases to fit the regression (else return raw).
#' @return numeric vector of inflammation-adjusted RBP (raw where markers NA).
brinda_adjust_rbp <- function(rbp, crp, agp = NULL, min_n = 30L) {
  rbp <- suppressWarnings(as.numeric(rbp))
  crp <- suppressWarnings(as.numeric(crp))
  # CRP-only mode: some surveys (e.g. TDHS 2010) assay CRP but not AGP. The
  # regression then has one inflammation term instead of two. This is a
  # DIFFERENT, weaker correction than the two-marker BRINDA — never silently mix
  # it with the two-marker version in a cross-country comparison without saying
  # so; see brinda_country_method().
  crp_only <- is.null(agp)
  agp <- if (crp_only) rep(NA_real_, length(rbp)) else suppressWarnings(as.numeric(agp))
  adj <- rbp
  ok <- is.finite(rbp) & is.finite(crp) & rbp > 0 & crp > 0
  if (!crp_only) ok <- ok & is.finite(agp) & agp > 0
  if (sum(ok) < min_n) return(adj)
  lrbp <- log(rbp[ok]); lcrp <- log(crp[ok])
  crp_ref <- as.numeric(stats::quantile(crp[ok], 0.10, na.rm = TRUE))
  if (crp_only) {
    b <- tryCatch(stats::coef(stats::lm(lrbp ~ lcrp)),
                  error = function(e) c(`(Intercept)` = NA, lcrp = 0))
    bC <- min(unname(b["lcrp"]), 0, na.rm = TRUE)
    corr <- bC * pmax(lcrp - log(crp_ref), 0)
  } else {
    lagp <- log(agp[ok])
    agp_ref <- as.numeric(stats::quantile(agp[ok], 0.10, na.rm = TRUE))
    b <- tryCatch(stats::coef(stats::lm(lrbp ~ lcrp + lagp)),
                  error = function(e) c(`(Intercept)` = NA, lcrp = 0, lagp = 0))
    bC <- min(unname(b["lcrp"]), 0, na.rm = TRUE)
    bA <- min(unname(b["lagp"]), 0, na.rm = TRUE)
    corr <- bC * pmax(lcrp - log(crp_ref), 0) + bA * pmax(lagp - log(agp_ref), 0)
  }
  adj[ok] <- exp(lrbp - corr)
  adj
}

#' Overwrite a VitA outcome's binary in-place with the uniform BRINDA-adjusted
#' VAD (RBP < 0.70). No-op for non-VitA outcomes or when the raw RBP/CRP/AGP
#' columns are absent. `d` must already be filtered to ONE population (children
#' or women). Called from build_outcome_dataset() (national + within-country +
#' area transport) AND from the LOCO pooling in transportability.R (individual
#' transport) so EVERY VitA path shares one definition (DC-H2).
#' @return d with d[[oc$binary]] replaced by the BRINDA VAD where applicable.
apply_brinda_vita_binary <- function(d, cc, oc, label = "") {
  if (is.null(oc$tag) || !grepl("vitA", oc$tag, ignore.case = TRUE)) return(d)
  if (is.null(oc$binary)) return(d)
  newbin <- brinda_vad_binary(d, cc, oc, label = label)
  if (is.null(newbin)) return(d)
  d[[oc$binary]] <- newbin
  d
}

#' Derive the harmonized VAD binary (adjusted RBP < 0.70) for one population.
#'
#' Single source of truth shared by apply_brinda_vita_binary() (national /
#' within-country / area paths) and compute_svy_admin2() (survey-weighted
#' Admin-2 prevalence), so every VitA path uses one definition. `d` must already
#' be filtered to ONE population.
#'
#' @return integer 0/1 vector the length of nrow(d), or NULL if the country's
#'   biomarker columns are unavailable (caller keeps the configured binary).
brinda_vad_binary <- function(d, cc, oc, cutoff = 0.70, label = "") {
  pop  <- if (grepl("^child", oc$tag)) "child" else "women"
  spec <- brinda_rbp_cols(cc$country)[[pop]]
  if (is.null(spec)) {
    warning(sprintf("[brinda]%s %s: no RBP/CRP(/AGP) column map for this country; ",
                    label, cc$country),
            "keeping the configured binary")
    return(NULL)
  }
  need <- unname(unlist(spec[intersect(names(spec), c("rbp", "crp", "agp"))]))
  if (!all(need %in% colnames(d))) {
    warning(sprintf("[brinda]%s %s %s: missing column(s) %s; keeping configured binary",
                    label, cc$country, oc$tag,
                    paste(setdiff(need, colnames(d)), collapse = ", ")))
    return(NULL)
  }

  adj <- brinda_adjust_rbp(d[[spec$rbp]], d[[spec$crp]],
                           if (is.null(spec$agp)) NULL else d[[spec$agp]])

  # Surveys that assayed the inflammation marker on a SUBSAMPLE leave the rest
  # of the rows uncorrected. Rather than mixing corrected and raw values in one
  # outcome, fall back to the survey agency's own adjusted biomarker on those
  # rows (declared per country as `fallback`). Rows with neither stay NA.
  n_fallback <- 0L
  if (!is.null(spec$fallback) && spec$fallback %in% colnames(d)) {
    no_corr <- is.na(d[[spec$crp]]) | !is.finite(suppressWarnings(as.numeric(d[[spec$crp]])))
    fb <- suppressWarnings(as.numeric(d[[spec$fallback]]))
    use <- no_corr & is.finite(fb)
    adj[use]  <- fb[use]
    adj[no_corr & !is.finite(fb)] <- NA_real_
    n_fallback <- sum(use)
  }

  newbin <- as.integer(adj < cutoff)
  cat(sprintf("  [brinda]%s %s %s: VAD = %s RBP<%.2f (%d/%d = %.1f%%)%s\n",
              label, cc$country, oc$tag, brinda_country_method(cc$country), cutoff,
              sum(newbin == 1, na.rm = TRUE), sum(!is.na(newbin)),
              100 * mean(newbin, na.rm = TRUE),
              if (n_fallback > 0)
                sprintf(" [%d rows on survey-provided adjustment]", n_fallback) else ""))
  newbin
}

#' Short label for the inflammation-adjustment method a country actually gets.
#' Report this in the manuscript: it is NOT the same method everywhere.
brinda_country_method <- function(country) {
  spec <- brinda_rbp_cols(country)
  if (is.null(spec)) return("configured (no BRINDA)")
  if (is.null(spec$child$agp)) "BRINDA CRP-only" else "BRINDA CRP+AGP"
}

# =============================================================================
# WS3 (2026-08): uniform inflammation adjustment for FERRITIN.
#
# Vitamin A was harmonized in DC-H2: apply_brinda_vita_binary() overwrites every
# country's VAD binary with one BRINDA CRP+AGP definition, and
# brinda_country_method() now returns "BRINDA CRP+AGP" for all four active
# countries. Iron was not. UNIFORM_TRANSPORT_TAGS (R/admin2_analysis.R:20)
# applies a uniform CUTOFF to iron, but it applies it to each country's own
# configured continuous column, and those are adjusted four different ways:
#
#   Gambia        gw_LogFerAdj       log-scale, survey-agency adjustment
#   Ghana         gw_*FerrAdjThurn   Thurnham
#   Sierra Leone  gw_cFerrAdj (children) / gw_wFerAdjBR1 (women)
#   Malawi        sf_reg             regression-adjusted serum ferritin
#
# A uniform cutoff on non-uniform adjustments is not a uniform outcome. Since
# the dominant LOCO failure mode is a national LEVEL offset, and raw child
# ferritin medians run from 11.3 (Gambia) to 71.4 (Sierra Leone) while Sierra
# Leone's children also carry far the highest CRP, adjustment heterogeneity is a
# candidate explanation that has to be measured rather than assumed.
#
# These functions are ADDITIVE. Nothing below is called from the production
# path; they are driven from scripts/run_uniform_brinda.R under the
# `uniform_brinda` scheme so the configured outcomes are unchanged.
# =============================================================================

#' BRINDA-adjust ferritin for inflammation (iron).
#'
#' Same regression correction as brinda_adjust_rbp(), with the sign reversed.
#' Inflammation RAISES ferritin (it is an acute-phase reactant), so the
#' coefficients are clamped to >= 0 and the correction lowers the adjusted
#' value. In brinda_adjust_rbp() inflammation depresses RBP, the coefficients
#' are clamped to <= 0, and the correction raises it. Clamping in the wrong
#' direction would let collinearity manufacture or erase deficiency, so the two
#' clamps are not interchangeable.
#'
#' @param fer,crp,agp numeric vectors (ferritin ug/L; CRP mg/L; AGP g/L), same
#'   length, one population group.
#' @param min_n minimum complete cases to fit the regression (else return raw).
#' @return numeric vector of inflammation-adjusted ferritin (raw where markers NA)
brinda_adjust_ferritin <- function(fer, crp, agp = NULL, min_n = 30L) {
  fer <- suppressWarnings(as.numeric(fer))
  crp <- suppressWarnings(as.numeric(crp))
  crp_only <- is.null(agp)
  agp <- if (crp_only) rep(NA_real_, length(fer)) else suppressWarnings(as.numeric(agp))
  adj <- fer
  ok <- is.finite(fer) & is.finite(crp) & fer > 0 & crp > 0
  if (!crp_only) ok <- ok & is.finite(agp) & agp > 0
  if (sum(ok) < min_n) return(adj)
  lfer <- log(fer[ok]); lcrp <- log(crp[ok])
  crp_ref <- as.numeric(stats::quantile(crp[ok], 0.10, na.rm = TRUE))
  if (crp_only) {
    b <- tryCatch(stats::coef(stats::lm(lfer ~ lcrp)),
                  error = function(e) c(`(Intercept)` = NA, lcrp = 0))
    bC <- max(unname(b["lcrp"]), 0, na.rm = TRUE)
    corr <- bC * pmax(lcrp - log(crp_ref), 0)
  } else {
    lagp <- log(agp[ok])
    agp_ref <- as.numeric(stats::quantile(agp[ok], 0.10, na.rm = TRUE))
    b <- tryCatch(stats::coef(stats::lm(lfer ~ lcrp + lagp)),
                  error = function(e) c(`(Intercept)` = NA, lcrp = 0, lagp = 0))
    bC <- max(unname(b["lcrp"]), 0, na.rm = TRUE)
    bA <- max(unname(b["lagp"]), 0, na.rm = TRUE)
    corr <- bC * pmax(lcrp - log(crp_ref), 0) + bA * pmax(lagp - log(agp_ref), 0)
  }
  adj[ok] <- exp(lfer - corr)
  adj
}

#' Per-country RAW ferritin / CRP / AGP column map.
#'
#' Raw, not the survey-adjusted columns the configs point at: the whole point of
#' the uniform scheme is to re-derive the adjustment from the same starting
#' point everywhere. Units verified as ferritin ug/L, CRP mg/L, AGP g/L across
#' all four countries.
#'
#' Malawi stores one biomarker column per marker for both populations and splits
#' by its `population` text column, so child and women entries are identical.
brinda_ferritin_cols <- function(country) {
  switch(country,
    "Gambia"       = list(child = list(fer = "gw_cFER",      crp = "gw_cCRP", agp = "gw_cAGP"),
                          women = list(fer = "gw_wFER",      crp = "gw_wCRP", agp = "gw_wAGP")),
    "Ghana"        = list(child = list(fer = "gw_cFerr",     crp = "gw_cCRP", agp = "gw_cAGP"),
                          women = list(fer = "gw_wFerr",     crp = "gw_wCRP", agp = "gw_wAGP")),
    "Sierra Leone" = list(child = list(fer = "gw_cFerritin", crp = "gw_cCRP", agp = "gw_cAGP"),
                          women = list(fer = "gw_wFerritin", crp = "gw_wCRP", agp = "gw_wAGP")),
    "Malawi"       = list(child = list(fer = "fer", crp = "crp", agp = "agp"),
                          women = list(fer = "fer", crp = "crp", agp = "agp")),
    NULL)
}

#' Which inflammation adjustment does a country's IRON outcome actually get?
#'
#' The analogue of brinda_country_method() for iron. Unlike vitamin A, the
#' configured iron columns are not harmonized, so this reports the configured
#' column rather than a single method name.
brinda_iron_method <- function(country, population = "child") {
  spec <- brinda_ferritin_cols(country)
  if (is.null(spec)) return("configured (no uniform ferritin map)")
  if (is.null(spec[[population]]$agp)) "BRINDA CRP-only" else "BRINDA CRP+AGP"
}

#' Derive the uniform iron-deficiency binary from BRINDA-adjusted ferritin.
#'
#' `d` must already be filtered to ONE population. Cutoffs are the WHO values
#' the configs already use: 12 ug/L in children, 15 ug/L in women.
#'
#' @return integer 0/1 vector, or NULL when the country's raw columns are absent
#'   (caller keeps the configured binary rather than silently substituting).
brinda_id_binary <- function(d, cc, oc, label = "") {
  pop  <- if (grepl("^child", oc$tag)) "child" else "women"
  spec <- brinda_ferritin_cols(cc$country)[[pop]]
  if (is.null(spec)) {
    warning(sprintf("[brinda-fe]%s %s: no raw ferritin/CRP/AGP map; keeping configured binary",
                    label, cc$country))
    return(NULL)
  }
  need <- unname(unlist(spec[intersect(names(spec), c("fer", "crp", "agp"))]))
  if (!all(need %in% colnames(d))) {
    warning(sprintf("[brinda-fe]%s %s %s: missing column(s) %s; keeping configured binary",
                    label, cc$country, oc$tag, paste(setdiff(need, colnames(d)), collapse = ", ")))
    return(NULL)
  }
  cutoff <- if (pop == "child") 12 else 15
  adj <- brinda_adjust_ferritin(d[[spec$fer]], d[[spec$crp]],
                                if (is.null(spec$agp)) NULL else d[[spec$agp]])
  newbin <- as.integer(adj < cutoff)
  cat(sprintf("  [brinda-fe]%s %s %s: ID = %s ferritin<%g (%d/%d = %.1f%%)\n",
              label, cc$country, oc$tag, brinda_iron_method(cc$country, pop), cutoff,
              sum(newbin == 1, na.rm = TRUE), sum(!is.na(newbin)),
              100 * mean(newbin, na.rm = TRUE)))
  list(binary = newbin, adjusted = adj, cutoff = cutoff)
}

#' Enumerate, per country and outcome, exactly which adjustment is applied.
#'
#' Reads the configs rather than asserting from memory, so the table cannot
#' drift from the code. `configured_*` describes what the pipeline uses today;
#' `uniform_*` describes what the uniform_brinda scheme would use instead.
#'
#' @param configs get_country_configs() output
#' @return data.frame, one row per country x outcome
adjustment_inventory <- function(configs = get_country_configs()) {
  rows <- list()
  for (cn in names(configs)) {
    cc <- configs[[cn]]
    for (on in names(cc$outcomes)) {
      oc  <- cc$outcomes[[on]]
      pop <- if (grepl("^child", oc$tag)) "child" else "women"
      is_vita <- grepl("vitA", oc$tag, ignore.case = TRUE)
      is_iron <- grepl("iron", oc$tag, ignore.case = TRUE)

      rbp_spec <- brinda_rbp_cols(cc$country)[[pop]]
      fer_spec <- brinda_ferritin_cols(cc$country)[[pop]]

      rows[[length(rows) + 1L]] <- data.frame(
        country              = cc$country,
        outcome              = oc$tag,
        population           = oc$population,
        configured_continuous = oc$continuous %||% NA_character_,
        configured_binary     = oc$binary %||% NA_character_,
        cutoff                = oc$cutoff %||% NA_real_,
        cutoff_scale          = oc$cutoff_scale %||% NA_character_,
        uniform_cutoff_applied = oc$tag %in% UNIFORM_TRANSPORT_TAGS,
        # Vitamin A is already re-derived at runtime; iron is not.
        adjustment_harmonized_today = is_vita,
        harmonized_method_today = if (is_vita) brinda_country_method(cc$country) else NA_character_,
        uniform_raw_marker   = if (is_vita) (rbp_spec$rbp %||% NA_character_)
                               else if (is_iron) (fer_spec$fer %||% NA_character_)
                               else NA_character_,
        uniform_crp          = if (is_vita) (rbp_spec$crp %||% NA_character_)
                               else if (is_iron) (fer_spec$crp %||% NA_character_)
                               else NA_character_,
        uniform_agp          = if (is_vita) (rbp_spec$agp %||% NA_character_)
                               else if (is_iron) (fer_spec$agp %||% NA_character_)
                               else NA_character_,
        uniform_method       = if (is_vita) brinda_country_method(cc$country)
                               else if (is_iron) brinda_iron_method(cc$country, pop)
                               else "not applicable (no inflammation adjustment defined)",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Per-country raw RBP / CRP / AGP column map for the uniform VitA transport
#' outcome (gw_ countries split by sex via gw_child_flag; Malawi by the
#' `population` text column).
#'
#' Each population entry is a list with:
#'   rbp      raw (unadjusted) RBP column, umol/L
#'   crp      CRP column, mg/L
#'   agp      AGP column, g/L — NULL when the survey did not assay AGP, which
#'            switches brinda_adjust_rbp() to its weaker CRP-only correction
#'   fallback optional: the survey agency's own adjusted RBP, used for rows
#'            where the inflammation marker was not assayed (subsample designs)
#'
#' Tanzania (TDHS 2010) is the CRP-only case: no AGP was measured, and CRP was
#' assayed on ~27% of the RBP sample, so most rows fall back to the DHS-supplied
#' `rbpadcrp`. Its VAD is therefore NOT method-identical to the other four
#' countries — see docs/transportability_loco_methods.md.
brinda_rbp_cols <- function(country) {
  gw_sexed <- list(child = list(rbp = "gw_cRBP", crp = "gw_cCRP", agp = "gw_cAGP"),
                   women = list(rbp = "gw_wRBP", crp = "gw_wCRP", agp = "gw_wAGP"))
  switch(country,
    "Gambia"       = gw_sexed,
    "Ghana"        = gw_sexed,
    "Sierra Leone" = gw_sexed,
    "Malawi"       = list(child = list(rbp = "rbp", crp = "crp", agp = "agp"),
                          women = list(rbp = "rbp", crp = "crp", agp = "agp")),
    "Tanzania"     = list(
      child = list(rbp = "gw_rbp_raw_umol", crp = "gw_crp", agp = NULL,
                   fallback = "gw_cRBPAdjCRP"),
      women = list(rbp = "gw_rbp_raw_umol", crp = "gw_crp", agp = NULL,
                   fallback = "gw_wRBPAdjCRP")),
    NULL)
}
