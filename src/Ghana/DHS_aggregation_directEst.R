# ------------------------------------------------------------------------------
# Hard-coded, MIS-friendly indicator creation for Ghana 2016 MIS
# - No surveyPrev built-ins
# - No sourcing external repo
# - Robust to missing variables: creates only what is possible and logs the rest
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rlang)
  library(tibble)
  library(stringr)
})

# --- Small helpers ------------------------------------------------------------

has_all <- function(df, vars) all(vars %in% names(df))
nz_ <- function(x) length(x) && !is.null(x)

as_num <- function(x) {
  # robust numeric coercion that tolerates labelled/factors/characters
  if (inherits(x, "haven_labelled")) {
    # avoid requiring 'labelled' / 'haven' at runtime
    suppressWarnings(as.numeric(as.character(x)))
  } else if (is.factor(x)) {
    suppressWarnings(as.numeric(as.character(x)))
  } else {
    suppressWarnings(as.numeric(x))
  }
}

add_log <- function(log, dataset, var, ok, msg = "") {
  bind_rows(log, tibble(dataset = dataset, var = var, ok = ok, message = msg))
}

# Weights & survey ids are needed later for prevalence mapping
add_design_cols <- function(df, which) {
  if (is.null(df) || !nrow(df)) return(df)
  out <- df
  if (which %in% c("PRdata", "HRdata")) {
    if ("hv005" %in% names(out) && !"wt" %in% names(out)) out <- mutate(out, wt = hv005 / 1e6)
    # add common ids if absent
    for (nm in c("hv001","hv002","hv024","hv025")) if (!nm %in% names(out)) out[[nm]] <- NA
  } else { # IR/KR/MR
    if ("v005" %in% names(out) && !"wt" %in% names(out)) out <- mutate(out, wt = v005 / 1e6)
    for (nm in c("v001","v002","v024","v025")) if (!nm %in% names(out)) out[[nm]] <- NA
  }
  out
}

# --- Classification helpers (used only if the recodes exist) ------------------

# Basic drinking-water service (very conservative; skips if hv201 absent)
# NOTE: Codes follow standard DHS recode; if codes differ in your MIS, this will just be skipped.
water_basic <- function(hv201) {
  # Improved sources (core DHS; may vary slightly by country)
  # 11 piped to dwelling, 12 yard/plot, 13 public tap, 31 tubewell/borehole,
  # 32 protected well, 41 protected spring, 51 tanker (often 'other'),
  # 61 surface water (UNIMPROVED) -- intentionally excluded
  improved <- c(11,12,13,14,31,32,41)  # 14 sometimes "neighbor/other piped"
  x <- as_num(hv201)
  out <- ifelse(is.na(x), NA_integer_,
                ifelse(x %in% improved, 1L,
                       ifelse(x %in% c(33,42,43,51,61,71,96,99), 0L, NA_integer_)))
  out
}

# Basic sanitation (improved AND not shared); needs hv205 (+ hv225 if present)
sanitation_basic <- function(hv205, hv225 = NULL) {
  # Improved tech if not shared: flush to sewer/septic/pit (11,12,13),
  # ventilated improved pit (21), pit latrine with slab (22), composting (41)
  improved <- c(11,12,13,21,22,41)
  tlet <- as_num(hv205)
  shared <- if (!is.null(hv225)) as_num(hv225) else NA_real_
  # hv225==1 means shared (not basic), 0 not shared (basic)
  imp <- ifelse(is.na(tlet), NA_integer_,
                ifelse(tlet %in% improved, 1L,
                       ifelse(tlet %in% c(14,15,23,31,51,61,71,96,99), 0L, NA_integer_)))
  if (!is.null(hv225)) {
    out <- ifelse(is.na(imp) | is.na(shared), NA_integer_,
                  ifelse(imp == 1L & shared == 0, 1L,
                         ifelse(imp == 1L & shared == 1, 0L, imp)))
  } else {
    # If we don't have 'shared', return "improved" as a proxy and log that it's not exact
    out <- imp
  }
  out
}

# Any bednet in household (needs hv227)
any_net_hh <- function(hv227) {
  x <- as_num(hv227)
  ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, ifelse(x %in% c(0, 8, 9), 0L, NA_integer_)))
}

# Child slept under *a* net last night (very common MIS item in PR: hml12)
u5_slept_any_net <- function(hv105, hml12) {
  age <- as_num(hv105)
  used <- as_num(hml12) # 1 yes, 0 no, 8/9 DK/missing (set to 0 per convention or NA; we use 0)
  out <- ifelse(is.na(age), NA_integer_,
                ifelse(age < 5, ifelse(used == 1, 1L,
                                       ifelse(used %in% c(0,8,9), 0L, NA_integer_)),
                       NA_integer_))
  out
}

# --- Main: add MIS indicators safely -----------------------------------------

augment_mis_indicators <- function(PRdata = NULL, HRdata = NULL, IRdata = NULL, KRdata = NULL) {

  log <- tibble(dataset = character(), var = character(), ok = logical(), message = character())
  added <- list(PRdata = character(), HRdata = character(), IRdata = character(), KRdata = character())

  # 0) Ensure weights / common IDs exist (harmless no-ops if already present)
  PRdata <- add_design_cols(PRdata, "PRdata")
  HRdata <- add_design_cols(HRdata, "HRdata")
  IRdata <- add_design_cols(IRdata, "IRdata")
  KRdata <- add_design_cols(KRdata, "KRdata")

  # 1) HR: Any net in household (hv227)
  if (nz_(HRdata)) {
    if ("hv227" %in% names(HRdata)) {
      HRdata <- HRdata %>% mutate(ml_hh_anynet = any_net_hh(hv227))
      log <- add_log(log, "HRdata", "ml_hh_anynet", TRUE, "")
      added$HRdata <- c(added$HRdata, "ml_hh_anynet")
    } else {
      log <- add_log(log, "HRdata", "ml_hh_anynet", FALSE, "missing hv227")
    }
  }

  # 2) HR: Basic water source (hv201) -> ws_source_basic
  if (nz_(HRdata)) {
    if ("hv201" %in% names(HRdata)) {
      HRdata <- HRdata %>% mutate(ws_source_basic = water_basic(hv201))
      log <- add_log(log, "HRdata", "ws_source_basic", TRUE, "classification conservative; country codes may vary")
      added$HRdata <- c(added$HRdata, "ws_source_basic")
    } else {
      log <- add_log(log, "HRdata", "ws_source_basic", FALSE, "missing hv201")
    }
  }

  # 3) HR: Basic sanitation (hv205 [+ hv225 if available]) -> ws_toi_basic
  if (nz_(HRdata)) {
    if ("hv205" %in% names(HRdata)) {
      HRdata <- HRdata %>%
        mutate(ws_toi_basic = sanitation_basic(hv205, if ("hv225" %in% names(HRdata)) hv225 else NULL))
      msg <- if ("hv225" %in% names(HRdata)) "" else "no hv225; using 'improved' as proxy"
      log <- add_log(log, "HRdata", "ws_toi_basic", TRUE, msg)
      added$HRdata <- c(added$HRdata, "ws_toi_basic")
    } else {
      log <- add_log(log, "HRdata", "ws_toi_basic", FALSE, "missing hv205")
    }
  }

  # 4) PR: U5 slept under any net last night (hml12) -> ml_u5_anynet_use
  if (nz_(PRdata)) {
    needs <- c("hv105", "hml12")
    if (has_all(PRdata, needs)) {
      PRdata <- PRdata %>% mutate(ml_u5_anynet_use = u5_slept_any_net(hv105, hml12))
      log <- add_log(log, "PRdata", "ml_u5_anynet_use", TRUE, "")
      added$PRdata <- c(added$PRdata, "ml_u5_anynet_use")
    } else {
      log <- add_log(log, "PRdata", "ml_u5_anynet_use", FALSE,
                     paste("missing", paste(setdiff(needs, names(PRdata)), collapse = ",")))
    }
  }

  # 5) IR: (Optional) currently‑pregnant women slept under a net (if present)
  # Often MISs have v213 (currently pregnant) but may not have a direct net-use item.
  # If PR has member-level net use but you want a women-only version, you can later
  # build it by linking PR <-> IR via v001/v002/v003; we keep it out for simplicity here.

  # Return updated objects and a compact “what we added” summary
  list(
    PRdata = PRdata,
    HRdata = HRdata,
    IRdata = IRdata,
    KRdata = KRdata,
    added_cols = added,
    log = log
  )
}

# --- Convenience: standardize ONE 0/1 variable for surveyPrev later -----------
# This matches the format surveyPrev::getDHSindicator() produces (cluster/HH/region/weight/strata/value).
standardize_binary <- function(df, which, value_var) {
  stopifnot(value_var %in% names(df))
  if (which %in% c("PRdata","HRdata")) {
    tibble(
      cluster     = df[["hv001"]],
      householdID = df[["hv002"]],
      v024        = df[["hv024"]],
      weight      = df[["hv005"]] / 1e6,
      strata      = if ("hv023" %in% names(df)) df[["hv023"]] else df[["hv025"]],
      value       = df[[value_var]]
    ) %>% filter(!is.na(value))
  } else {
    tibble(
      cluster     = df[["v001"]],
      householdID = df[["v002"]],
      v024        = df[["v024"]],
      weight      = df[["v005"]] / 1e6,
      strata      = if ("v023" %in% names(df)) df[["v023"]] else df[["v025"]],
      value       = df[[value_var]]
    ) %>% filter(!is.na(value))
  }
}


d <- readRDS(here::here("data/DHS/dhs_Ghana_2016.RDS"))

PRdata <- d$Ghana_2016$PRdata
HRdata <- d$Ghana_2016$HRdata
IRdata <- d$Ghana_2016$IRdata
KRdata <- d$Ghana_2016$KRdata


colnames(PRdata)
colnames(HRdata)
colnames(IRdata)
colnames(KRdata)


ncol(PRdata)
ncol(HRdata)
ncol(IRdata)
ncol(KRdata)

# 1) Add MIS indicators directly to your data.frames
res <- augment_mis_indicators(PRdata = PRdata, HRdata = HRdata, IRdata = IRdata, KRdata = KRdata)

# Reassign to keep working with the updated objects
PRdata <- res$PRdata
HRdata <- res$HRdata
IRdata <- res$IRdata
KRdata <- res$KRdata

ncol(PRdata)
ncol(HRdata)
ncol(IRdata)
ncol(KRdata)

# 2) Quick diagnostics: what got created?
cat("\nAdded columns:\n")
print(vapply(res$added_cols, length, integer(1)))
cat("\nAdded to HRdata:\n"); print(res$added_cols$HRdata)
cat("\nAdded to PRdata:\n"); print(res$added_cols$PRdata)

# 3) Inspect skips/failures so we know what to add next
cat("\nCreation log (failures at bottom):\n")
print(res$log %>% arrange(ok, dataset, var))

# 4) (Optional) Standardize an indicator for admin-2 direct estimates with surveyPrev later:
#    Example: HR household has any net (ml_hh_anynet), if it was created
if ("ml_hh_anynet" %in% names(HRdata)) {
  ml_anynet_std <- standardize_binary(HRdata, which = "HRdata", value_var = "ml_hh_anynet")
}

#    Example: PR children <5 slept under any net (ml_u5_anynet_use), if it was created
if ("ml_u5_anynet_use" %in% names(PRdata)) {
  ml_u5_use_std <- standardize_binary(PRdata, which = "PRdata", value_var = "ml_u5_anynet_use")
}

# 5) When you’re ready to compute admin-2 prevalence, plug one of the standardized
#    data frames into surveyPrev::directEST(data=..., cluster.info=..., admin=2).
#    (See surveyPrev vignette for details.) :contentReference[oaicite:0]{index=0}
