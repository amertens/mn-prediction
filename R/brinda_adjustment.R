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
brinda_adjust_rbp <- function(rbp, crp, agp, min_n = 30L) {
  rbp <- suppressWarnings(as.numeric(rbp))
  crp <- suppressWarnings(as.numeric(crp))
  agp <- suppressWarnings(as.numeric(agp))
  adj <- rbp
  ok <- is.finite(rbp) & is.finite(crp) & is.finite(agp) &
        rbp > 0 & crp > 0 & agp > 0
  if (sum(ok) < min_n) return(adj)
  lrbp <- log(rbp[ok]); lcrp <- log(crp[ok]); lagp <- log(agp[ok])
  crp_ref <- as.numeric(stats::quantile(crp[ok], 0.10, na.rm = TRUE))
  agp_ref <- as.numeric(stats::quantile(agp[ok], 0.10, na.rm = TRUE))
  b <- tryCatch(stats::coef(stats::lm(lrbp ~ lcrp + lagp)),
                error = function(e) c(`(Intercept)` = NA, lcrp = 0, lagp = 0))
  bC <- min(unname(b["lcrp"]), 0, na.rm = TRUE)
  bA <- min(unname(b["lagp"]), 0, na.rm = TRUE)
  corr <- bC * pmax(lcrp - log(crp_ref), 0) + bA * pmax(lagp - log(agp_ref), 0)
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
  pop  <- if (grepl("^child", oc$tag)) "child" else "women"
  cols <- brinda_rbp_cols(cc$country)[[pop]]
  if (is.null(cols) || !all(cols %in% colnames(d))) {
    warning(sprintf("[brinda]%s %s %s: RBP/CRP/AGP cols missing (%s); keeping configured binary",
                    label, cc$country, oc$tag, paste(cols, collapse = ",")))
    return(d)
  }
  adj    <- brinda_adjust_rbp(d[[cols[1]]], d[[cols[2]]], d[[cols[3]]])
  newbin <- as.integer(adj < 0.70)
  d[[oc$binary]] <- newbin
  cat(sprintf("  [brinda]%s %s %s: VAD = uniform BRINDA RBP<0.70 (%d/%d = %.1f%%)\n",
              label, cc$country, oc$tag, sum(newbin == 1, na.rm = TRUE),
              sum(!is.na(newbin)), 100 * mean(newbin, na.rm = TRUE)))
  d
}

#' Per-country raw RBP/CRP/AGP column map for the uniform VitA transport outcome.
#' (gw_ countries split by sex via gw_child_flag; Malawi by the `population`
#' text column.) Used by compute_svy_admin2() when building the harmonized
#' women_vitA / child_vitA transport outcome.
brinda_rbp_cols <- function(country) {
  switch(country,
    "Gambia"       = list(child = c("gw_cRBP","gw_cCRP","gw_cAGP"),
                          women = c("gw_wRBP","gw_wCRP","gw_wAGP")),
    "Ghana"        = list(child = c("gw_cRBP","gw_cCRP","gw_cAGP"),
                          women = c("gw_wRBP","gw_wCRP","gw_wAGP")),
    "Sierra Leone" = list(child = c("gw_cRBP","gw_cCRP","gw_cAGP"),
                          women = c("gw_wRBP","gw_wCRP","gw_wAGP")),
    "Malawi"       = list(child = c("rbp","crp","agp"),
                          women = c("rbp","crp","agp")),
    NULL)
}
