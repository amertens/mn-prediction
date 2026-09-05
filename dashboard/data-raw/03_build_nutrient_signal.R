# =============================================================================
# dashboard/data-raw/03_build_nutrient_signal.R
#
# Builds the bundle behind the "What tracks which nutrient" tab.
#
# WHY THIS TAB EXISTS
# -------------------
# The dashboard could show where deficiency is predicted, how models compare and
# how the pipeline works, but nowhere could a user see WHICH INDICATORS TRACK
# WHICH NUTRIENT - which is the question a nutrition programme actually brings.
# The answer is also the project's most defensible output: the associations are
# nutrient-specific, mechanistically coherent, and agree in direction across all
# four countries.
#
# EVIDENCE STANDARD: SIGN AGREEMENT, NOT P-VALUES
# -----------------------------------------------
# Under correction for having examined 451 predictors at 14-87 districts per
# country, few individual indicators clear conventional significance; the audit
# measured the per-cell screen's median power at 0.005. The same indicator with
# the same sign in four independent national surveys - countries with different
# diets, agriculture and survey teams - is both stronger evidence and a more
# honest claim. The tab therefore leads with the agreement count and carries the
# family-wise p only as a secondary column.
#
#   Rscript dashboard/data-raw/03_build_nutrient_signal.R
# -> dashboard/data/nutrient_signal.rds
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
setwd(here::here())
PROBE <- here::here("results", "tables", "signal_probes")
DASH  <- here::here("dashboard", "data")
MD <- read.csv(here::here("data","covariates","harmonized",
                          "predictors_admin2_shared_metadata.csv"),
               stringsAsFactors = FALSE)
dm <- stats::setNames(MD$domain, MD$column)

rd <- function(f) {
  p <- file.path(PROBE, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}
OUT_LABEL <- c(child_iron = "Iron — children", child_vitA = "Vitamin A — children",
               women_iron = "Iron — women", women_vitA = "Vitamin A — women",
               women_folate = "Folate — women", women_b12 = "Vitamin B12 — women",
               child_zinc = "Zinc — children", women_zinc = "Zinc — women")

tidy <- function(d, scale_lab) {
  if (is.null(d)) return(NULL)
  d |> filter(group %in% names(OUT_LABEL)) |>
    transmute(outcome = unname(OUT_LABEL[group]), scale = scale_lab,
              indicator = predictor,
              domain = ifelse(is.na(dm[predictor]), NA_character_, unname(dm[predictor])),
              meta_z, countries = k_countries, agree = sign_agree,
              direction = ifelse(meta_z > 0, "more deficiency", "less deficiency"),
              p_fwer)
}
pred <- bind_rows(tidy(rd("p1_admin1_scan_predictors.csv"), "Deficiency prevalence"),
                  tidy(rd("p4_admin1_continuous_predictors.csv"), "Biomarker concentration"))
dom  <- bind_rows(tidy(rd("p1_admin1_scan_domains.csv"), "Deficiency prevalence"),
                  tidy(rd("p4_admin1_continuous_domains.csv"), "Biomarker concentration")) |>
  rename(domain_name = indicator) |> select(-domain)

# Headline, read off the scans rather than curated from expectation.
# CORRECTED 2026-09-02. The previous version of this table asserted that legume
# and fruit-and-vegetable consumption tracked LESS child vitamin A deficiency,
# and cattle ownership LESS child iron deficiency. Both are backwards. In
# p4_admin1_continuous_predictors.csv the outcome is negated (see the script
# header) so a POSITIVE meta_z means MORE deficiency, and legumes (+4.87, 4/4)
# and cattle (+4.40, 3/3) are both positive. The signs replicate in the
# independent binary scan (p1). The honest reading is that these are markers of
# rural subsistence agro-ecology, which is what an ecological fallacy looks like
# when you check it - so the table now says that instead of implying that
# district associations are dietary effects.
headline <- data.frame(
  outcome = c("Vitamin A — children", "Iron — children", "Iron — women",
              "Vitamin A — women", "Vitamin B12 — women"),
  tracks = c("Legume consumption; warm night-time temperatures",
             "Cattle ownership; cereal-dominant cropping (root crops track less)",
             "Warm night-time temperatures; soil aluminium heterogeneity",
             "Child wasting in the same district; rainfall",
             "Vegetable production; goat ownership (soil organic carbon tracks less)"),
  direction = c("more legumes → MORE deficiency",
                "more cattle → MORE deficiency",
                "→ more deficiency",
                "more wasting → more; more rain → less",
                "→ more deficiency"),
  agreeing = c("4 of 4", "4 of 4", "4 of 4", "4 of 4", "3 of 3"),
  reading = c(
    "Counterintuitive, and replicated in both scans. At district level legume-eating and cattle-keeping mark rural subsistence areas, which are poorer and more deficient. Do NOT read this as 'legumes cause deficiency' - it is the opposite sign to the individual-level diet relationship.",
    "Same pattern as vitamin A: cattle ownership marks the rural districts, not protected children. Cereal-dominant cropping tracking more deficiency is consistent with staple-heavy diets; root crops track less.",
    "Agro-ecology markers. Soil is a marker of farming system here, not a demonstrated causal pathway.",
    "The one clearly interpretable pair, because both are genuine district properties: child wasting and women's vitamin A deficiency co-locate and may be targetable together.",
    "Also agro-ecological rather than dietary. Fewer countries measure B12, so this row rests on three."),
  stringsAsFactors = FALSE)

saveRDS(list(headline = headline, predictors = pred, domains = dom,
             note = paste("Evidence is cross-country sign agreement, not",
                          "significance: under correction for 451 predictors",
                          "at 14-87 districts, the per-cell screen had median",
                          "power 0.005. Several replicated associations run",
                          "OPPOSITE to individual-level nutrition biology.",
                          "These are district-level markers of deficient areas,",
                          "useful for targeting, and must not be read as",
                          "individual dietary effects."),
             build_time = Sys.time()),
        file.path(DASH, "nutrient_signal.rds"))
cat(sprintf("-- Nutrient signal --\n  headline %d | predictor rows %d | domain rows %d\n",
            nrow(headline), nrow(pred %||% data.frame()), nrow(dom %||% data.frame())))
