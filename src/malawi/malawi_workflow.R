# Workflow to clean, merge, and prepare DHS data for Malawi
suppressPackageStartupMessages({
  library(here)
})

source(here::here("src/DHS/country_workflow_utils.R"))

run_country_workflow(
  country = "Malawi",
  iso2 = "MW",
  iso3 = "MWI",
  year = 2019,
  dhs_path = here("data/DHS/dhs_Malawi_2019.RDS"),
  gee_path = here("data/DHS/clean/DHS_GEE_merge_MW_2019_mkts.csv"),
  price_path = here("data/food_price/wfp_food_prices_mwi.csv"),
  malaria_prefix = "MWI"
)
