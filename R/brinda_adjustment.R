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
