
rm(list=ls())
# load  packages
library(rdhs)
library(ggplot2)
library(here)

metadata<-read.csv(here("data/DHS/clean/dhs_indicators_metadata.csv"))

metadata$ShortName[grepl("deficiency", metadata$ShortName)]

rdhs::dhs_countries()

#Test
d2015 <- dhs_data(countryIds = "MW",
                  indicatorIds = "FE_FRTR_W_A15",
                  surveyYearStart = 2015,
                  breakdown = "subnational")
d2015

#NOTE! RDHS only includes core indicators, not micronutrients





# ---- Setup ----
library(haven)
library(dplyr)
library(labelled)
library(survey)
library(stringr)
library(purrr)
library(broom)

# ---------- Helper utils ----------
pick_var <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit)) hit[1] else NA_character_
}

has_all <- function(df, vars) all(vars %in% names(df))

# BRINDA internal regression adjustment:
# outcome_var is "log ferritin" or "log RBP" depending on biomarker supplied.
# biom: one of "ferritin", "rbp"
# ref values from BRINDA papers (can tweak): CRP_ref=0.1 mg/L, AGP_ref=0.63 g/L
brinda_adjust <- function(dat, biom_value, crp, agp,
                          biom = c("ferritin","rbp"),
                          CRP_ref = 0.1, AGP_ref = 0.63) {
  biom <- match.arg(biom)
  d <- dat %>%
    mutate(
      biom_val = !!sym(biom_value),
      crp_val  = !!sym(crp),
      agp_val  = !!sym(agp)
    ) %>%
    filter(is.finite(biom_val), is.finite(crp_val), is.finite(agp_val)) %>%
    mutate(
      lnb   = log(biom_val),
      lncrp = log(pmax(crp_val, 1e-6)),
      lnagp = log(pmax(agp_val, 1e-6))
    )
  if (nrow(d) < 50) return(rep(NA_real_, nrow(dat))) # not enough data to fit

  fit <- lm(lnb ~ lncrp + lnagp, data = d)
  # Adjust to reference inflammation:
  # ln(b)_adj = ln(b) - b1*(lnCRP - lnCRPref) - b2*(lnAGP - lnAGPref)
  # (Same formula for ferritin and RBP; the signs will be learned by the model)
  lnCRP_ref <- log(CRP_ref)
  lnAGP_ref <- log(AGP_ref)

  adjs <- rep(NA_real_, nrow(dat))
  # predict per row using estimated coefs
  co <- coef(fit)
  b1 <- unname(co["lncrp"])
  b2 <- unname(co["lnagp"])

  # build lnb for all rows (even if some missing -> keep NA)
  lnb_all   <- log(dat[[biom_value]])
  lncrp_all <- log(pmax(dat[[crp]], 1e-6))
  lnagp_all <- log(pmax(dat[[agp]], 1e-6))

  lnb_adj <- lnb_all - b1 * (lncrp_all - lnCRP_ref) - b2 * (lnagp_all - lnAGP_ref)
  exp(lnb_adj) # returns adjusted concentration
}

# altitude adjustment for hemoglobin (simple CDC/DHS table implementation).
# If altitude is missing, returns unadjusted Hb.
hb_altitude_adjust <- function(hb, altitude_m) {
  if (all(is.na(altitude_m))) return(hb)
  # Formula (approx): adjustment (g/dL) = 0 if alt < 1000 m
  #  = -0.032*alt + 0.022*(alt^2) for alt in km above 1.  (rounded DHS/CDC)
  # Use DHS standard table approximation:
  alt_km_over1 <- pmax(altitude_m/1000 - 1, 0)
  adj <- 0.032 * alt_km_over1*10 + 0.0 # simple linear approx per 1000m (tame)
  # Use subtraction because Hb decreases at altitude-adjusted sea-level equivalent
  hb - adj
}

# generic survey-weighted prevalence helper
svy_prev <- function(design, var) {
  if (!var %in% names(design$variables)) return(tibble(var = var, prev = NA_real_, se = NA_real_))
  est <- suppressWarnings(svymean(~get(var), design, na.rm = TRUE))
  tibble(var = var,
         prev = as.numeric(coef(est)),
         se   = as.numeric(SE(est)))
}

# ---------- Load PR file (Malawi 2015–16 example) ----------
# Replace with your actual .DTA path (PR file)
# Example: pr_path <- "MWPR7HFL.DTA"


#NOTE! I need to get access to the Malawi biomarker dataset

malawi_dhs <- readRDS(here("data/DHS/dhs_Malawi_2015.RDS"))

pr <- malawi_dhs$Malawi_2015$PRdata

library(haven)
library(labelled)
library(dplyr)

# Function to extract names + labels
get_var_labels <- function(df) {
  data.frame(
    name  = names(df),
    label = sapply(df, function(x) {
      lab <- attr(x, "label", exact = TRUE)
      if (is.null(lab)) "" else lab
    }),
    stringsAsFactors = FALSE
  )
}

# Example usage:
# df <- read_sav("yourfile.sav")
labels_df <- get_var_labels(PR)
labels_df2 <- get_var_labels(malawi_dhs$Malawi_2015$BRdata)
head(labels_df)

labels_df[grepl("crp",labels_df$label),]
labels_df[grepl("vit",labels_df$label),]
labels_df[grepl("rbp",labels_df$label),]



# ---------- Identify variables (flexible name matching) ----------
# weights / strata / psu
wgt  <- pick_var(pr, c("hv005","hv005a"))
strt <- pick_var(pr, c("hv022","v022"))
psu  <- pick_var(pr, c("hv021","v021"))

# region (admin-1) for subnational estimates
region <- pick_var(pr, c("hv024","v024","region"))

# age (in months) & sex to subset children 6–59 months
age_m <- pick_var(pr, c("hc1","b19","age", "hvidx_age_months"))
sex   <- pick_var(pr, c("hc27","sex","hv104"))

# altitude (meters) if available at cluster level
alt_m <- pick_var(pr, c("hv040","altitude"))

# Biomarkers (common DHS/HW candidates)
hb_child   <- pick_var(pr, c("hc1","hw53","hgb","hgbc"))       # g/dL
hb_women   <- pick_var(pr, c("ha1","hw54","hgbw"))             # g/dL


#NOTE! Need to debug and find these variables
labels(pr)

rbp_child  <- pick_var(pr, c("hw55","rbpc"))
rbp_women  <- pick_var(pr, c("hw56","rbpw"))

crp_child  <- pick_var(pr, c("hw57","crpc"))
crp_women  <- pick_var(pr, c("hw58","crpw"))

agp_child  <- pick_var(pr, c("hw59","agpc"))
agp_women  <- pick_var(pr, c("hw60","agpw"))

fer_child  <- pick_var(pr, c("hw70","ferritinc","ferc"))
fer_women  <- pick_var(pr, c("hw71","ferritinw","ferw"))

stfr_child <- pick_var(pr, c("hw72","stfrc"))
stfr_women <- pick_var(pr, c("hw73","stfrw"))

zinc_child <- pick_var(pr, c("hw74","zincc"))
zinc_women <- pick_var(pr, c("hw75","zincw"))

iodine_child <- pick_var(pr, c("hu1","uicc"))
iodine_women <- pick_var(pr, c("hu2","uicw"))

# ---------- Construct child 6–59 mo analytic set ----------



children <- pr %>%
  mutate(
    wt        = .data[[wgt]] / 1e6,
    strata    = .data[[strt]],
    psu       = .data[[psu]],
    region    = .data[[region]],
    age_m     = .data[[age_m]],
    sex       = .data[[sex]],
    altitude_m= if (!is.na(alt_m)) .data[[alt_m]] else NA_real_,

    hb        = if (!is.na(hb_child))   as.numeric(.data[[hb_child]])   else NA_real_,
    rbp       = if (!is.na(rbp_child))  as.numeric(.data[[rbp_child]])  else NA_real_,
    crp       = if (!is.na(crp_child))  as.numeric(.data[[crp_child]])  else NA_real_,
    agp       = if (!is.na(agp_child))  as.numeric(.data[[agp_child]])  else NA_real_,
    ferritin  = if (!is.na(fer_child))  as.numeric(.data[[fer_child]])  else NA_real_,
    stfr      = if (!is.na(stfr_child)) as.numeric(.data[[stfr_child]]) else NA_real_,
    zinc      = if (!is.na(zinc_child)) as.numeric(.data[[zinc_child]]) else NA_real_,
    uic       = if (!is.na(iodine_child)) as.numeric(.data[[iodine_child]]) else NA_real_
  ) %>%
  filter(between(age_m, 6, 59))


summary(children$hb)
summary(children$rbp)
summary(children$crp)

# adjust Hb for altitude if present
children <- children %>%
  mutate(hb_adj = if (!all(is.na(altitude_m))) hb_altitude_adjust(hb, altitude_m) else hb)

# BRINDA adjustments where possible
if (!any(is.na(children$crp)) && !any(is.na(children$agp)) && !all(is.na(children$ferritin))) {
  children$ferritin_adj <- brinda_adjust(children, "ferritin", "crp", "agp", biom = "ferritin")
} else {
  children$ferritin_adj <- NA_real_
}

if (!any(is.na(children$crp)) && !any(is.na(children$agp)) && !all(is.na(children$rbp))) {
  children$rbp_adj <- brinda_adjust(children, "rbp", "crp", "agp", biom = "rbp")
} else {
  children$rbp_adj <- NA_real_
}

# ---------- Child indicators (DHS-style) ----------
# Anemia (WHO child): <11.0 g/dL
children <- children %>%
  mutate(
    anemia = ifelse(!is.na(hb_adj), hb_adj < 11.0, NA),
    # Iron deficiency (after adjustment if available): ferritin < 12 µg/L
    fer_use = ifelse(!is.na(ferritin_adj), ferritin_adj, ferritin),
    iron_def = ifelse(!is.na(fer_use), fer_use < 12, NA),
    ida = ifelse(!is.na(iron_def) & !is.na(anemia), iron_def & anemia, NA),
    # Vitamin A deficiency via RBP (after adjustment if available): RBP < 0.7 µmol/L
    rbp_use = ifelse(!is.na(rbp_adj), rbp_adj, rbp),
    vad = ifelse(!is.na(rbp_use), rbp_use < 0.7, NA),
    # Zinc deficiency (VERY context-dependent; placeholder cutoff 65 µg/dL AM non-fasting etc.) — optional:
    zinc_def = ifelse(!is.na(zinc), NA, NA),  # left as NA unless you set context-specific cutoffs
    # Iodine deficiency at individual level is not standard; usually report population median UIC.
    uic_lt100 = ifelse(!is.na(uic), uic < 100, NA)
  )

summary(children$rbp_use)
summary(children$ferritin)
summary(children$ferritin_adj)


# survey design
des_child <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt,
                       data = children, nest = TRUE)

# National prevalences
child_inds <- c("anemia","iron_def","ida","vad","uic_lt100")
child_nat <- bind_rows(lapply(child_inds, \(v) svy_prev(des_child, v))) %>%
  mutate(pop = "Children 6-59m", level = "National")

# Regional prevalences (admin-1)
# --- Children: regional prevalences (admin-1)
child_reg <- purrr::map_dfr(child_inds, function(v) {
  regs <- unique(stats::na.omit(children$region))
  purrr::map_dfr(regs, function(r) {
    des_sub <- subset(des_child, region == r)
    est <- tryCatch(
      svymean(as.formula(paste0("~", v)), des_sub, na.rm = TRUE),
      error = function(e) NULL
    )
    tibble::tibble(
      var    = v,
      region = r,
      prev   = if (is.null(est)) NA_real_ else as.numeric(coef(est)),
      se     = if (is.null(est)) NA_real_ else as.numeric(SE(est)),
      pop    = "Children 6-59m",
      level  = "Region"
    )
  })
})


# ---------- Women 15–49 analytic set ----------
# Women commonly live in IR file, but when biomarkers are merged to PR you can subset by sex/age
# sex codes often: 1 = male, 2 = female (check labels)
# Decide column names up front
age_var <- pick_var(pr, c("hv105","v012","age"))   # years
# sex code typically: 1 = male, 2 = female in PR
# all other *var objects (wgt, strt, psu, region, alt_m, hb_women, etc.) already defined earlier

women <- pr %>%
  mutate(
    wt        = .data[[wgt]] / 1e6,
    strata    = .data[[strt]],
    psu       = .data[[psu]],
    region    = .data[[region]],
    age_y     = if (!is.na(age_var)) as.numeric(.data[[age_var]]) else NA_real_,
    sex       = .data[[sex]],
    altitude_m= if (!is.na(alt_m)) .data[[alt_m]] else NA_real_,
    hb        = if (!is.na(hb_women))   as.numeric(.data[[hb_women]])   else NA_real_,
    rbp       = if (!is.na(rbp_women))  as.numeric(.data[[rbp_women]])  else NA_real_,
    crp       = if (!is.na(crp_women))  as.numeric(.data[[crp_women]])  else NA_real_,
    agp       = if (!is.na(agp_women))  as.numeric(.data[[agp_women]])  else NA_real_,
    ferritin  = if (!is.na(fer_women))  as.numeric(.data[[fer_women]])  else NA_real_,
    stfr      = if (!is.na(stfr_women)) as.numeric(.data[[stfr_women]]) else NA_real_,
    zinc      = if (!is.na(zinc_women)) as.numeric(.data[[zinc_women]]) else NA_real_,
    uic       = if (!is.na(iodine_women)) as.numeric(.data[[iodine_women]]) else NA_real_
  ) %>%
  filter(sex == 2, dplyr::between(age_y, 15, 49)) %>%
  # always initialize adjustment columns to correct length (even if 0 rows)
  mutate(
    ferritin_adj = NA_real_,
    rbp_adj      = NA_real_
  )

# If there are rows, compute BRINDA adjustments only where complete cases exist
if (nrow(women) > 0) {
  # ferritin
  idx_f <- stats::complete.cases(women[, c("ferritin","crp","agp")])
  if (sum(idx_f) >= 30) {  # need some data to fit the model
    adj_f <- brinda_adjust(women[idx_f, ], "ferritin", "crp", "agp", biom = "ferritin")
    women$ferritin_adj[idx_f] <- adj_f
  }
  # RBP
  idx_r <- stats::complete.cases(women[, c("rbp","crp","agp")])
  if (sum(idx_r) >= 30) {
    adj_r <- brinda_adjust(women[idx_r, ], "rbp", "crp", "agp", biom = "rbp")
    women$rbp_adj[idx_r] <- adj_r
  }
}

# proceed as before...
women <- women %>%
  mutate(hb_adj = if (!all(is.na(altitude_m))) hb_altitude_adjust(hb, altitude_m) else hb,
         fer_use = dplyr::if_else(!is.na(ferritin_adj), ferritin_adj, ferritin),
         rbp_use = dplyr::if_else(!is.na(rbp_adj), rbp_adj, rbp),
         anemia  = dplyr::if_else(!is.na(hb_adj), hb_adj < 12.0, NA),
         iron_def= dplyr::if_else(!is.na(fer_use), fer_use < 15, NA),
         ida     = dplyr::if_else(!is.na(iron_def) & !is.na(anemia), iron_def & anemia, NA),
         vad     = dplyr::if_else(!is.na(rbp_use), rbp_use < 0.7, NA),
         zinc_def = ifelse(!is.na(zinc), NA, NA),
         uic_lt100 = ifelse(!is.na(uic), uic < 100, NA))





des_women <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt,
                       data = women, nest = TRUE)

women_inds <- c("anemia","iron_def","ida","vad","uic_lt100")
women_nat <- bind_rows(lapply(women_inds, \(v) svy_prev(des_women, v))) %>%
  mutate(pop = "Women 15-49", level = "National")


# --- Women: regional prevalences (admin-1)
women_reg <- purrr::map_dfr(women_inds, function(v) {
  regs <- unique(stats::na.omit(women$region))
  purrr::map_dfr(regs, function(r) {
    des_sub <- subset(des_women, region == r)
    est <- tryCatch(
      svymean(as.formula(paste0("~", v)), des_sub, na.rm = TRUE),
      error = function(e) NULL
    )
    tibble::tibble(
      var    = v,
      region = r,
      prev   = if (is.null(est)) NA_real_ else as.numeric(coef(est)),
      se     = if (is.null(est)) NA_real_ else as.numeric(SE(est)),
      pop    = "Women 15-49",
      level  = "Region"
    )
  })
})

# ---------- Results ----------
nat <- bind_rows(child_nat, women_nat)
reg <- bind_rows(child_reg, women_reg)

print(nat)
print(reg)
