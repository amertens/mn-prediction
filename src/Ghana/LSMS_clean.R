# ── 0.  Setup ───────────────────────────────────────────────────────────────
rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(survey)
library(fastDummies)   # for categorical dummies
library(here)

hh  <- read_dta(here("data/LSMS/g7aggregates_hhlevel.dta"))  %>% clean_names()
ind <- read_dta(here("data/LSMS/g7aggregates_indlevel.dta")) %>% clean_names()

# NB: a few manuals you uploaded earlier have expired from the session cache.
# If you still need me to parse them for extra variable definitions, please re‑upload.

# ── 1.  Households: identify variable types ─────────────────────────────────
id_hh   <- c("hid","clust","country","survemo","surveyr",
             "region","loc2","loc5","loc5g","loc7","loc7g",
             "rururb","ez","wta_s")
hh_num  <- hh %>% select(where(is.numeric)) %>% select(-any_of(id_hh)) %>% names()
#hh_cat  <- hh %>% select(-where(is.numeric)) %>% select(-any_of(id_hh)) %>% names()

# ── 2.  Build survey design & summarise HOUSEHOLD numerics ──────────────────
des_hh      <- svydesign(ids = ~1, weights = ~wta_s, data = hh)

hh_num_mean <- svyby(
  as.formula(paste("~", paste(hh_num, collapse = "+"))),
  ~region, des_hh, svymean, na.rm = TRUE, vartype = NULL) %>%
  rename_with(~ paste0("hh_", .x), -region)

# check categorical - seems to be just IDs
 hh %>% select(-where(is.numeric))



# ── 5.  INDIVIDUAL‑level summaries (numeric + categorical) ─────────────────
id_ind  <- c("hid","pid","country","surveyr","district",
             "region","loc2","wta_s")
ind_num <- ind %>% select(where(is.numeric)) %>% select(-any_of(id_ind)) %>% names()

des_ind       <- svydesign(ids = ~1, weights = ~wta_s, data = ind)

ind_num_mean  <- svyby(
  as.formula(paste("~", paste(ind_num, collapse = "+"))),
  ~region, des_ind, svymean, na.rm = TRUE, vartype = NULL) %>%
  rename_with(~ paste0("ind_", .x), -region)



# ── 6.  Merge household & individual feature sets ──────────────────────────
admin1_ml_features <- hh_num_mean %>%
  left_join(ind_num_mean, by = "region")

#remove near zero variance columns
admin1_ml_features <- admin1_ml_features[,-caret::nzv(admin1_ml_features)]

#code regions
table(admin1_ml_features$region)
admin1_ml_features <- admin1_ml_features %>% mutate(admin1 = factor(region, levels=c(1:10),
                                    labels=c("Western", "Central", "Greater Accra", "Volta", "Eastern", "Ashanti", "Brong Ahafo", "Northern", "Upper East", "Upper West")))
table(admin1_ml_features$admin1)

#append column names for source
colnames(admin1_ml_features) <- paste0("lsms_", colnames(admin1_ml_features))


# save
saveRDS(admin1_ml_features, here("data/LSMS/Ghana_LSMS_clean.RDS"))


