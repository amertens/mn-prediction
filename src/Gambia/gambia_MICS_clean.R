
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



d <-read_dta(here("data/MICS/Gambia 2018/Gambia_cleaned.dta")) %>% clean_names()

#drop near zero variance columns
d <- d [,-caret::nzv(d)]

meta <- makeVlist(d)
write.csv(meta,
          file = here("metadata/MICS/gambia_misc_metadata.csv"),
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

# #clean region
# # value       label
# # 1      Banjul
# # 2    Kanifing
# # 3     Brikama
# # 4  Mansakonko
# # 5     Kerewan
# # 6     Kuntaur
# # 7 Janjanbureh
# # 8       Basse
# d_trim$region <- factor(d_trim$region, levels=1:8, labels=c("Banjul","Kanifing","Brikama","Mansakonko","Kerewan","Kuntaur","Janjanbureh","Basse"))

# Harmonize MICS region codes to Admin1 labels

# Original MICS codes:
# 1 Banjul
# 2 Kanifing
# 3 Brikama        -> Western
# 4 Mansakonko    -> Lower River
# 5 Kerewan       -> North Bank
# 6 Kuntaur       -> Maccarthy Island
# 7 Janjanbureh   -> Maccarthy Island
# 8 Basse         -> Upper River

d_trim$region <- factor(
  d_trim$region,
  levels = 1:8,
  labels = c(
    "Banjul",
    "Kanifing",
    "Western",
    "Lower River",
    "North Bank",
    "Maccarthy Island",
    "Maccarthy Island",
    "Upper River"
  )
)


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

write_csv(region_summary, here("data/MICS/mics_gambia_2018_region_summary.csv"))
