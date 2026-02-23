# Workflow to clean, merge, and prepare DHS data for The Gambia
suppressPackageStartupMessages({
  library(here)
})

source(here::here("src/country_workflow_utils.R"))

# run_country_workflow(
#   country = "The Gambia",
#   iso2 = "GM",
#   iso3 = "GMB",
#   year = 2019,
#   dhs_path = here("data/DHS/dhs_TheGambia_2019.RDS"),
#   gee_path = here("data/DHS/clean/DHS_GEE_merge_GM_2019_mkts.csv"),
#   price_path = here("data/food_price/wfp_food_prices_gmb.csv"),
#   malaria_prefix = "GMB"
# )


country = "The Gambia"
iso2 = "GM"
iso3 = "GMB"
year = 2019
dhs_path = here("data/DHS/dhs_Gambia_2019.RDS")
#gee_path = here("data/DHS/clean/DHS_GEE_merge_GH_2019_mkts.csv") $this is ghana
price_path = here("data/food_price/wfp_food_prices_gmb.csv")
malaria_prefix = "GMB"

df <- clean_country_dhs(dhs_path)
df <- attach_compiled_gps(df, gps_path, country, year)


gee_df <- read.csv(gee_path) %>%
  rename(cluster = DHSCLUST)

dim(df)
dim(gee_df)
df$
df2 <- left_join(df, gee_df, by = c("cluster","hh"))
dim(df2)

if (!is.null(gee_path)) {
  df <- attach_gee_features(df, gee_path)
}

if (!is.null(price_path)) {
  df <- attach_food_prices(df, price_path)
}

if (!is.null(malaria_prefix)) {
  df <- attach_malaria_atlas(df, raster_dir, malaria_prefix)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

outfile <- file.path(output_dir, sprintf("dhs_%s_%s_clean.RDS", gsub(" ", "", country), year))
saveRDS(df, outfile)
message(sprintf("Saved cleaned data to %s", outfile))

