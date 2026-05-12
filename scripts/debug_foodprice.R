library(janitor)
raw <- read.csv("data/food_price/wfp_food_prices_gmb.csv", stringsAsFactors = FALSE)
raw <- raw[!grepl("^#", raw[[1]]), ]
raw$year <- as.integer(substr(raw$date, 1, 4))
raw <- raw[!is.na(raw$year) & raw$year >= 2017 & raw$year <= 2018, ]
raw$usdprice <- suppressWarnings(as.numeric(raw$usdprice))
raw <- raw[!is.na(raw$usdprice) & raw$usdprice > 0, ]
raw <- raw[!grepl("non-food|miscellaneous", raw$category, ignore.case = TRUE), ]
raw$commodity_clean <- janitor::make_clean_names(raw$commodity)
raw$wfp_admin <- raw$admin2
cat("nrow:", nrow(raw), "\n")
cat("unique admin:", length(unique(raw$wfp_admin)), "\n")
cat("unique commodity:", length(unique(raw$commodity_clean)), "\n")
agg <- aggregate(usdprice ~ wfp_admin + commodity_clean, data = raw, FUN = mean, na.rm = TRUE)
cat("agg nrow:", nrow(agg), "\n")
wide <- tidyr::pivot_wider(agg, names_from = commodity_clean, values_from = usdprice)
cat("wide:", nrow(wide), "rows x", ncol(wide), "cols\n")
cat("commodities:", paste(head(setdiff(colnames(wide), "wfp_admin"), 10), collapse=", "), "\n")
