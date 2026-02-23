
# ── 0.  Packages ──────────────────────────────────────────────────────────
rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(survey)
library(srvyr)       # wrapper around survey
library(fastDummies)
library(purrr)
library(stringr)
library(here)
library(haven)
library(caret)
source(paste0(here::here(),"/src/0-functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_functions.R"))
source(paste0(here::here(),"/src/DHS/DHS_variable_recode.R"))

options(survey.lonely.psu = "adjust")



d <-read_dta(here("data/MICS/Ghana 2017/Ghana_cleaned.dta")) %>% clean_names()

#drop near zero variance columns
d <- d [,-caret::nzv(d)]

meta <- makeVlist(d)
write.csv(meta,
          file = here("metadata/MICS/misc_metadata.csv"),
          row.names = FALSE,
          na = "")



# ── 2 · Define thematic keywords ──────────────────────────────────────────
themes <- list(
  nutrition        = c("breast", "diet", "meal", "weight",
                       "height", "stunt", "wast", "underweight",
                       "vitamin", "iron", "folic"),
  wash             = c("water", "sanitation", "toilet", "latrine",
                       "hygiene", "handwash", "soap","coli","cook","cooking","vessel","storage",
                       "handwashing","detergent","energy"),
  food_security    = c("food security", "hunger", "fies", "fcs"),
  socio_economic   = c("education", "literacy", "wealth", "quintile",
                       "income", "electricity", "asset", "floor",
                       "roof", "fuel", "household", "number","percentile group",
                       "material","member","internet","goats","sheep","chickens","wscore"),
  health_services  = c("age","malaria", "fever", "diarrhoea", "diarrhea",
                       "ari", "anc", "antenatal", "delivery",
                       "postnatal", "vaccin", "immuni","iodization", "mosquitos")
)

# Build ONE regex:  \b(keyword1|keyword2|…)\b
key_regex <- paste0("\\b(", paste(unlist(themes), collapse = "|"), ")\\b")

# ── 3 · Pick ≤ 100 candidate variables from metadata ─────────────────────
vars_kept <- meta %>%
  mutate(label = tolower(label)) %>%         # normalise once
  filter(str_detect(label, key_regex), !(name %in% c("hhname","hhaddr")))

vars_dropped <- meta %>%
  mutate(label = tolower(label)) %>%         # normalise once
  filter(!str_detect(label, key_regex))

vars_kept <- vars_kept %>%   # keep only relevant labels
  #slice_head(n = 100) %>%                    # hard cap of 100
  pull(name)

# OPTIONAL: preview the short‑list
print(vars_kept)

# ── 4 · Keep only those variables (plus region & weight) in the micro‑data ─
needed <- c("region", grep("weight$", names(d), value = TRUE), vars_kept, "treat_any", "san_shared","animals","floor")
d_trim <- d %>% select(any_of(needed)) %>% clean_names()

#clean region
# value         label
# 1       WESTERN
# 2       CENTRAL
# 3 GREATER ACCRA
# 4         VOLTA
# 5       EASTERN
# 6       ASHANTI
# 7   BRONG AHAFO
# 8      NORTHERN
# 9    UPPER EAST
# 10    UPPER WEST
d_trim$region <- factor(d_trim$region, levels=1:10, labels=c("Western",'Central',"Greater Accra","Volta","Eastern","Ashanti","Brong Ahafo","Northern","Upper East","Upper West"))

# identify the weight column we just kept
wt_var <- grep("weight$", names(d_trim), value = TRUE)[1]

# ── 5 · Labelled integers → factors so they become categoricals ───────────
d_trim <- d_trim %>%
  mutate(
    across(
      where(~ is.labelled(.x) && !is.numeric(.x)),  # *categorical* labelled
      as_factor
    ),
    across(                                         # strip labels from numeric vars
      where(~ is.labelled(.x) &&  is.numeric(.x)),
      ~ haven::zap_labels(.x, preserve_attributes = FALSE)
    )
  )

svy <- d_trim %>%
  as_survey_design(weights = !!sym(wt_var))   # srvyr magic

# ----------------------- 3 · NUMERIC FEATURES ------------------------------

num_vars <- d_trim %>%
  select(-all_of(wt_var), -region) %>%  # drop weight + group var
  select(where(is.numeric)) %>%
  names()

num_summary <- svy %>%
  group_by(region) %>%
  summarise(
    across(
      all_of(num_vars),
      ~ survey_mean(., na.rm = TRUE, vartype = "se"),  # mean + SE
      .names = "{.col}"
    ),
    .groups = "drop"
  )


# ── 8 · Categorical summaries (weighted proportions) ─────────────────────
cat_vars <- d_trim %>% select(-region) %>%
  select(-where(is.numeric)) %>%
  names()


cat_dummy <- d_trim %>%
  select( all_of(wt_var), all_of(cat_vars)) %>%
  fastDummies::dummy_cols(remove_first_dummy = FALSE,
                          remove_selected_columns = TRUE) %>%
  clean_names()
colnames(cat_dummy)

cat_dummy <- data.frame(region=d_trim$region, cat_dummy)
head(cat_dummy)

des_cat <- svydesign(ids = ~1,
                     weights = as.formula(paste0("~", wt_var)),
                     data    = cat_dummy)

rhs <- setdiff(colnames(cat_dummy), c("region", wt_var))
cat_tbl <- svyby(
  as.formula(paste("~", paste(rhs, collapse = "+"))),
  ~region, des_cat, svymean, na.rm = TRUE, vartype = NULL
)

region_summary <- left_join(num_summary, cat_tbl, by = "region")

colnames(region_summary) <- paste0("mics_",colnames(region_summary))

write_csv(region_summary, here("data/MICS/mics_ghana_2017_region_summary.csv"))
