#!/usr/bin/env Rscript
# Test WFP food price extraction and inject into cached external predictors
library(here)
library(janitor)
library(tidyr)

.empty_result <- function(admin2_names = character(0)) {
  data.frame(Admin2 = admin2_names, stringsAsFactors = FALSE)
}

# ── Inline the function to avoid sourcing the full external_data.R ──
extract_wfp_foodprices_standalone <- function(gadm_admin2_names, gadm_admin1_names,
                                               country_name, survey_year,
                                               price_dir = here::here("data", "food_price")) {
  country_codes <- c(
    "Gambia" = "gmb", "Ghana" = "gha",
    "Sierra Leone" = "sle", "Malawi" = "mwi"
  )
  code <- country_codes[[country_name]]
  if (is.null(code)) {
    cat(sprintf("[WFP prices] No mapping for '%s'\n", country_name))
    return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))
  }

  csv_path <- file.path(price_dir, paste0("wfp_food_prices_", code, ".csv"))
  if (!file.exists(csv_path)) {
    cat(sprintf("[WFP prices] File not found: %s\n", csv_path))
    return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))
  }

  cat(sprintf("[WFP prices] Loading %s\n", basename(csv_path)))
  raw <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  raw <- raw[!grepl("^#", raw[[1]]), ]
  cat(sprintf("[WFP prices] %d rows after removing metadata\n", nrow(raw)))

  raw$year <- as.integer(substr(raw$date, 1, 4))
  raw <- raw[!is.na(raw$year) & raw$year >= (survey_year - 1L) & raw$year <= survey_year, ]
  cat(sprintf("[WFP prices] %d rows in %d-%d\n", nrow(raw), survey_year - 1L, survey_year))

  if (nrow(raw) == 0) return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))

  raw$usdprice <- suppressWarnings(as.numeric(raw$usdprice))
  raw <- raw[!is.na(raw$usdprice) & raw$usdprice > 0, ]

  if ("category" %in% colnames(raw)) {
    raw <- raw[!grepl("non-food|miscellaneous", raw$category, ignore.case = TRUE), ]
  }
  cat(sprintf("[WFP prices] %d food price observations\n", nrow(raw)))
  if (nrow(raw) == 0) return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))

  # Determine admin level
  admin_col <- if ("admin2" %in% colnames(raw) &&
                   sum(nchar(trimws(raw$admin2)) > 0) > nrow(raw) / 2) {
    "admin2"
  } else if ("admin1" %in% colnames(raw)) {
    cat("[WFP prices] admin2 sparse — using admin1\n")
    "admin1"
  } else {
    return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))
  }

  raw <- raw[nchar(trimws(raw[[admin_col]])) > 0, ]
  if (nrow(raw) == 0) return(data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE))

  raw$commodity_clean <- tolower(gsub("[^a-zA-Z0-9]+", "_", raw$commodity))
  raw$commodity_clean <- gsub("_+", "_", raw$commodity_clean)
  raw$commodity_clean <- gsub("^_|_$", "", raw$commodity_clean)
  raw$wfp_admin <- raw[[admin_col]]

  agg <- stats::aggregate(
    usdprice ~ wfp_admin + commodity_clean,
    data = raw, FUN = function(x) mean(x, na.rm = TRUE)
  )

  wide <- tidyr::pivot_wider(agg, names_from = commodity_clean, values_from = usdprice)

  price_cols <- setdiff(colnames(wide), "wfp_admin")
  for (pc in price_cols) {
    colnames(wide)[colnames(wide) == pc] <- paste0("wfp_price_", pc)
  }
  cat(sprintf("[WFP prices] %d admin areas x %d commodities\n", nrow(wide), length(price_cols)))

  # Match names
  wfp_names <- wide$wfp_admin
  match_target <- if (admin_col == "admin1") gadm_admin1_names else gadm_admin2_names

  name_map <- setNames(rep(NA_character_, length(gadm_admin2_names)), gadm_admin2_names)

  for (i in seq_along(match_target)) {
    nm <- match_target[i]
    if (nm %in% wfp_names) name_map[gadm_admin2_names[i]] <- nm
  }

  unmatched_targets <- unique(match_target[is.na(name_map)])
  for (nm in unmatched_targets) {
    fuzzy <- agrep(nm, wfp_names, max.distance = 0.2, value = TRUE)
    if (length(fuzzy) >= 1) {
      best <- fuzzy[which.min(abs(nchar(fuzzy) - nchar(nm)))]
      idx <- which(match_target == nm & is.na(name_map))
      name_map[gadm_admin2_names[idx]] <- best
    }
  }

  n_matched <- sum(!is.na(name_map))
  cat(sprintf("[WFP prices] Matched %d/%d GADM units\n", n_matched, length(gadm_admin2_names)))

  result <- data.frame(Admin2 = gadm_admin2_names, stringsAsFactors = FALSE)
  result$.wfp_admin <- name_map[gadm_admin2_names]

  wfp_price_cols <- grep("^wfp_price_", colnames(wide), value = TRUE)
  result <- merge(result, wide[, c("wfp_admin", wfp_price_cols), drop = FALSE],
                  by.x = ".wfp_admin", by.y = "wfp_admin", all.x = TRUE, sort = FALSE)
  result$.wfp_admin <- NULL

  result
}


# ── Update each country's external cache ──
countries <- list(
  list(name = "Gambia", year = 2018L, cache = "gambia_external_predictors.rds"),
  list(name = "Ghana", year = 2017L, cache = "ghana_external_predictors.rds"),
  list(name = "Sierra Leone", year = 2013L, cache = "sierraleone_external_predictors.rds"),
  list(name = "Malawi", year = 2016L, cache = "malawi_external_predictors.rds")
)

cache_dir <- here::here("data", "external_cache")

for (cc in countries) {
  cache_file <- file.path(cache_dir, cc$cache)
  if (!file.exists(cache_file)) {
    cat(sprintf("\n[SKIP] Cache not found: %s\n", cc$cache))
    next
  }

  cat(sprintf("\n=== %s ===\n", cc$name))
  ext <- readRDS(cache_file)

  # Remove any existing wfp_price_ columns
  old_price_cols <- grep("^wfp_price_", colnames(ext), value = TRUE)
  if (length(old_price_cols) > 0) {
    ext <- ext[, !colnames(ext) %in% old_price_cols, drop = FALSE]
    cat(sprintf("  Removed %d existing wfp_price_ columns\n", length(old_price_cols)))
  }

  fp <- extract_wfp_foodprices_standalone(
    ext$Admin2, ext$Admin1, cc$name, cc$year
  )

  # Merge price columns into existing cache
  price_cols <- grep("^wfp_price_", colnames(fp), value = TRUE)
  if (length(price_cols) > 0) {
    ext <- merge(ext, fp[, c("Admin2", price_cols), drop = FALSE],
                 by = "Admin2", all.x = TRUE, sort = FALSE)
    cat(sprintf("  Added %d price columns to cache (%d total cols now)\n",
                length(price_cols), ncol(ext)))

    saveRDS(ext, cache_file)
    cat(sprintf("  Saved: %s\n", basename(cache_file)))
  } else {
    cat("  No price columns extracted\n")
  }
}

cat("\nDone.\n")
