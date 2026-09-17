# =============================================================================
# scripts/covariates/build_hces_diet_block.R   [HC-01, 2026-09-15]
#
# HOUSEHOLD DIET FROM THE HOUSEHOLD CONSUMPTION AND EXPENDITURE SURVEYS
#
# The one domain every ablation said was missing: what households eat. Four
# HCES, one per country, within a year of three of the four micronutrient
# surveys (Sierra Leone's is five years after):
#
#   Malawi        IHS4 2016-17 (LSMS-ISA)   module G 7-day food recall, 141 items
#                 with purchase value / own production / gifts; consumption
#                 aggregate; EA coordinates (displaced)      -> TA (Admin-2)
#   The Gambia    IHS 2015/16                7-day food recall, 165 items with
#                 purchase value / own production; non-food purchases at
#                 7-day / 1-month / 3-month / 12-month recall; no aggregate
#                 -> district (Admin-2, name-matched to GADM)
#   Sierra Leone  SLIHS 2011 (HC-02)         diary purchases at five visits, 166
#                 labelled items, own consumption per visit; non-food diary and
#                 12-month purchases; no aggregate -> district (Admin-2, by code).
#                 Two years before the 2013 survey; replaces SLIHS 2018 (five after)
#   Ghana         GLSS7 2016/17              household expenditure aggregates
#                 by COICOP class (food, non-purchased food, ...); region only
#                 -> 10 old regions -> 16 (crosswalk) -> districts
#
# Household indicators (weighted means per area; the household table per
# country is written for audit):
#   hces_food_share          food (purchased + own + gifts) / total consumption
#   hces_own_prod_share      own production share of food: by VALUE where the
#                            aggregate splits it (Sierra Leone, Ghana), by
#                            ITEM (share of consumed items with any own-
#                            production quantity) where it does not (Malawi,
#                            The Gambia) - recorded per country in `assumption`
#   hces_log_cons_pae_rel    log real consumption per adult equivalent (per
#                            capita for Malawi, the aggregate's own basis),
#                            centred on the national household-weighted mean
#                            (the aggregates are in four currencies)
#   hces_hdds                mean count of the 12 HDDS food groups consumed
#                            in the 7-day recall (Malawi, The Gambia)
#   hces_any_<group>         share of households consuming any: asf (meat,
#                            fish, eggs, dairy), fish, meat, eggs, dairy,
#                            pulses_nuts, fruit, veg, dgl (dark green leafy),
#                            vita_fv (vitamin-A-rich fruit and vegetables)
#   hces_asf_purchase_share  animal-source share of food PURCHASE value
# Item -> HDDS group maps are written to metadata/hces_food_groups_<country>.csv
# so the classification can be reviewed line by line.
#
# Nothing is filtered for coverage: Ghana has NA for the recall-based
# indicators and that is recorded in the metadata.
#
#   Rscript -e "source('scripts/covariates/build_hces_diet_block.R')"
# -> data/covariates/harmonized/predictors_admin2_hces.csv (+ _metadata.csv)
# -> data/covariates/harmonized/hces_household_<Country>.csv
# -> metadata/hces_food_groups_{Malawi,Gambia,SierraLeone}.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(haven); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
targets::tar_source("R")
RA   <- "data/RA_2026-09/extracted"
HDIR <- "data/covariates/harmonized"
MIN_HH <- 15L
kk  <- function(x) tolower(gsub("[^a-z]", "", tolower(as.character(x))))
num <- function(x) suppressWarnings(as.numeric(unclass(x)))
wm  <- function(x, w) { ok <- is.finite(x) & is.finite(w) & w > 0; if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
lab <- function(x) { l <- attr(x, "labels"); if (is.null(l)) return(NULL); stats::setNames(names(l), as.character(l)) }

SH <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
spine <- SH[, c("country", "Admin1", "Admin2")]
GROUPS <- c("cereals", "roots", "veg", "fruit", "meat", "eggs", "fish", "pulses_nuts", "dairy", "oils", "sugar", "misc")

# ── food-group classification from item labels ──────────────────────────────
classify_item <- function(label, code = NA_integer_, block = NA_character_) {
  s <- tolower(label)
  has <- function(p) grepl(p, s, perl = TRUE)
  # order matters: composite baked goods before fish ("meat pie/fish pie"), garden egg
  # before eggs, oils before pulses ("groundnut oil"), roots before sugar ("sweet
  # potato") and before fruit ("plantain, cooking banana")
  # condiments, stimulants and drink concentrates are HDDS "spices, condiments,
  # beverages" whatever plant they come from
  if (has("meat pie|fish pie|sausage roll")) return("cereals")
  # SLIHS 2011 labels (HC-02): composite meals, bread listed with cassava bread, guinea corn, avocado as "butter pear"
  if (has("^bread|guinea corn|\\btuo\\b")) return("cereals")
  if (has("foofoo|acheke|starchy")) return("roots")
  if (has("avocado")) return("fruit")
  if (has("cashew nut|^other pulses$")) return("pulses_nuts")
  if (has("guinea fowl|wild bird|wild game|\\bgame\\b|frog|\\bmeet\\b")) return("meat")
  if (has("pumkin")) return("veg")
  if (has("kola|cola nut|locust bean|netetu|neteetu|dawadawa|tomato puree|tomato paste|tomato sauce|hot sauce|chilli|black pepper|dry pepper|powder pepper|pepper powder|garlic|squash \\(|sobo|concentrate")) return("misc")
  if (has("ice cream")) return("sugar")
  if (has("garden egg|aubergine|egg ?plant")) return("veg")
  if (has("leaf\\b|leaves")) return("veg")
  if (has("fish|bonga|tilapia|shrimp|oyster|crab|couta|tenny|sardine|seafood|snail|barracuda|grouper|ladyfish")) return("fish")
  if (has("\\begg")) return("eggs")
  if (has("milk|yogh|yogurt|cheese|chambiko|vitalait|cream\\b|dairy|infant feeding formula|baby milk")) return("dairy")
  if (has("beef|goat|mutton|sheep|pork|chicken|poultry|duck|\\bmeat|rabbit|mice|insect|termite|caterpillar|sausage|offal|liver")) return("meat")
  if (has("oil\\b|oils\\b|margarine|butter\\b|\\bfats?\\b|ghee|lard|mayonnaise|palm nut|palm kernel") && !has("peanut butter|groundnut paste")) return("oils")
  if (has("other fruits|cashew|cabaa|\\bkaba\\b")) return("fruit")
  if (has("\\bbeans?\\b|\\bpeas?\\b|pigeon ?pea|cowpea|chickpea|groundnut|peanut|soya|soy\\b|\\bnuts?\\b|lentil|nzama|nandolo|\\bseeds?\\b|sesame|benniseed")) return("pulses_nuts")
  if (has("cassava|\\byam|potato|cocoyam|plantain|masimbi|tuber|\\broot|taro|fufu|attieke|gari|chips")) return("roots")
  if (has("mango|banana|citrus|orange|lemon|lime\\b|pineapple|ananas|papaya|paw ?-? ?paw|guava|avocado|apple|fruit|watermelon|melon|grape|baobab|malambe|masau|dates|coconut|tamarind|daharr|dakhar|plum|berry|pear\\b")) return("fruit")
  if (has("onion|cabbage|\\brape\\b|tanaposi|nkhwani|leafy|tomato|cucumber|pumpkin|okra|therere|mushroom|vegetable|pepper|carrot|sorrel|bitter|spinach|lettuce|jakatu|bissap|bisap|kren|salad|jute|moringa|amaranth|kale")) return("veg")
  if (has("sugar|honey|jams?\\b|jelly|sweets|candy|chocolate|chewing gum|mint stick")) return("sugar")
  if (has("maize|rice|millet|sorghum|findi|fonio|wheat|flour|bread|bun|scone|pasta|spaghetti|macaroni|cereal|porridge|couscous|noodle|grain|biscuit|cake|doughnut|mandazi|samosa|popcorn|zikondamoyo|nkate|meal eaten|restaurant")) return("cereals")
  if (has("salt|spice|seasoning|maggi|cube|yeast|baking|sauce|vinegar|curry|mustard|tea|coffee|cocoa|milo|drink|juice|water|beer|wine|liquor|spirits|stout|kachasu|thobwa|maheu|chibuku|soda|soft|beverage|tobacco|cigarette|ice|mint")) return("misc")
  if (!is.na(block)) return(block)
  "misc"
}
dgl_item  <- function(label) grepl("\\brape\\b|tanaposi|nkhwani|green leaf|leafy|leaves|spinach|sorrel|bitter leaf|cassava leaf|potato leaf|jute|moringa|bissap|bisap|kren|pumpkin leaf|amaranth|kale|chinese cabbage", tolower(label), perl = TRUE)
# vitamin-A-rich fruit and vegetables (red palm oil is vitamin-A rich but is not a fruit or vegetable)
vita_item <- function(label) grepl("orange sweet potato|pumpkin|pumkin|mango|papaya|pawpaw|paw-paw|paw paw|paw - paw|carrot|apricot|butternut", tolower(label)) | dgl_item(label)
asf_groups <- c("meat", "fish", "eggs", "dairy")

hh_from_items <- function(items, id, group, consumed, value_purch, own_qty, dgl, vita) {
  # items: long data frame; returns one row per household with the diet indicators
  it <- data.frame(id = id, group = group, consumed = consumed, value = value_purch, own = own_qty, dgl = dgl, vita = vita, stringsAsFactors = FALSE)
  it <- it[!is.na(it$id), ]
  it$consumed <- is.finite(it$consumed) & it$consumed == 1
  it$value[!is.finite(it$value) | it$value < 0] <- 0
  it$own <- is.finite(it$own) & it$own > 0
  g <- it |> group_by(id) |>
    summarise(hdds = n_distinct(group[consumed]),
              any_asf = any(consumed & group %in% asf_groups), any_fish = any(consumed & group == "fish"),
              any_meat = any(consumed & group == "meat"), any_eggs = any(consumed & group == "eggs"),
              any_dairy = any(consumed & group == "dairy"), any_pulses_nuts = any(consumed & group == "pulses_nuts"),
              any_fruit = any(consumed & group == "fruit"), any_veg = any(consumed & group == "veg"),
              any_dgl = any(consumed & dgl), any_vita_fv = any(consumed & vita),
              n_items_consumed = sum(consumed), own_prod_item_share = if (sum(consumed)) mean(own[consumed]) else NA_real_,
              food_purch_value = sum(value), asf_purchase_share = if (sum(value) > 0) sum(value[group %in% asf_groups]) / sum(value) else NA_real_,
              .groups = "drop")
  as.data.frame(g)
}

# ── MALAWI IHS4 ─────────────────────────────────────────────────────────────
cat("\n[Malawi IHS4 2016-17]\n")
mw_path <- file.path(RA, "LSMS/MWI_2016/malawi_aggregated_hhlevel.dta")
mw_labs <- {
  r <- haven::read_dta(mw_path, n_max = 1); l <- sapply(r, function(x) { a <- attr(x, "label"); if (is.null(a)) "" else a })
  data.frame(var = names(l), label = unname(l), stringsAsFactors = FALSE)
}
g01 <- mw_labs[grepl("^hh_g01_", mw_labs$var), ]
g01$code <- as.integer(sub("^hh_g01_", "", g01$var)); g01$item <- sub(" - Did you.*$", "", g01$label)
blk <- function(code) c(`1` = "cereals", `2` = "roots", `3` = "pulses_nuts", `4` = "veg", `5` = "meat", `6` = "fruit", `7` = "dairy", `8` = "misc", `9` = "misc")[as.character(code %/% 100)]
g01$block <- ifelse(g01$code >= 5000, "fish", blk(pmin(g01$code, 999)))
g01$group <- mapply(function(l, c, b) classify_item(l, c, b), g01$item, g01$code, g01$block)
g01$dgl <- dgl_item(g01$item); g01$vita <- vita_item(g01$item)
write.csv(g01[, c("code", "item", "block", "group", "dgl", "vita")], "metadata/hces_food_groups_Malawi.csv", row.names = FALSE)
cat(sprintf("  %d items classified: %s\n", nrow(g01), paste(names(table(g01$group)), table(g01$group), collapse = " ")))
cols <- c("case_id", "ea_id", "district", "hh_wgt", "lat_modified", "lon_modified", "rexpagg", "hhsize", "adulteq",
          paste0("rexp_cat", sprintf("%02d", 1:12)), g01$var, sub("g01", "g05", g01$var), sub("g01", "g06a", g01$var))
cols <- intersect(cols, mw_labs$var)
mw <- haven::read_dta(mw_path, col_select = all_of(cols)) |> as.data.frame()
codes <- g01$code
long <- data.frame(id = rep(mw$case_id, times = length(codes)), code = rep(codes, each = nrow(mw)),
                   consumed = as.vector(sapply(g01$var, function(v) num(mw[[v]]))),
                   value = as.vector(sapply(sub("g01", "g05", g01$var), function(v) if (v %in% names(mw)) num(mw[[v]]) else NA_real_)),
                   own = as.vector(sapply(sub("g01", "g06a", g01$var), function(v) if (v %in% names(mw)) num(mw[[v]]) else NA_real_)))
long$consumed <- ifelse(long$consumed == 1, 1, 0)
long$group <- g01$group[match(long$code, g01$code)]; long$dgl <- g01$dgl[match(long$code, g01$code)]; long$vita <- g01$vita[match(long$code, g01$code)]
H_mw <- hh_from_items(long, long$id, long$group, long$consumed, long$value, long$own, long$dgl, long$vita)
tot <- rowSums(mw[, paste0("rexp_cat", sprintf("%02d", 1:12))], na.rm = TRUE)
base <- data.frame(id = mw$case_id, w = num(mw$hh_wgt), lon = num(mw$lon_modified), lat = num(mw$lat_modified), district_code = num(mw$district),
                   food_share = ifelse(tot > 0, num(mw$rexp_cat01) / tot, NA_real_),
                   log_cons_pae = log(pmax(num(mw$rexpagg), 1)), stringsAsFactors = FALSE)
H_mw <- left_join(base, H_mw, by = "id"); H_mw$own_prod_share <- H_mw$own_prod_item_share
# EA points -> GADM TA (Admin-2) polygons
poly <- sf::st_transform(load_gadm_cached("MWI", level = 2), 4326)
pts <- sf::st_as_sf(H_mw[is.finite(H_mw$lon) & is.finite(H_mw$lat), c("id", "lon", "lat")], coords = c("lon", "lat"), crs = 4326)
ix <- vapply(sf::st_within(pts, poly, sparse = TRUE), function(z) if (length(z)) z[1] else NA_integer_, integer(1))
pts$Admin1 <- as.character(poly$NAME_1)[ix]; pts$Admin2 <- as.character(poly$NAME_2)[ix]
H_mw <- left_join(H_mw, sf::st_drop_geometry(pts)[, c("id", "Admin1", "Admin2")], by = "id")
cat(sprintf("  %d households, %d with a TA (%.1f%%), %d TAs, weighted food share %.2f, HDDS %.1f, any ASF %.2f\n",
            nrow(H_mw), sum(!is.na(H_mw$Admin2)), 100 * mean(!is.na(H_mw$Admin2)), n_distinct(H_mw$Admin2), wm(H_mw$food_share, H_mw$w), wm(H_mw$hdds, H_mw$w), wm(H_mw$any_asf, H_mw$w)))
H_mw$country <- "Malawi"

# ── THE GAMBIA IHS 2015/16 ───────────────────────────────────────────────────
cat("\n[The Gambia IHS 2015/16]\n")
gm_dir <- file.path(RA, "LSMS/GMB_2015")
fd <- haven::read_dta(file.path(gm_dir, "Part B Section 1A-Food_consumption expenditure.dta"))
il <- lab(fd$s1aq1); items <- data.frame(code = as.integer(names(il)), item = unname(il), stringsAsFactors = FALSE)
items$group <- mapply(function(l, c) classify_item(l, c, NA_character_), items$item, items$code)
items$dgl <- dgl_item(items$item); items$vita <- vita_item(items$item)
write.csv(items, "metadata/hces_food_groups_Gambia.csv", row.names = FALSE)
cat(sprintf("  %d items classified: %s\n", nrow(items), paste(names(table(items$group)), table(items$group), collapse = " ")))
fd$code <- num(fd$s1aq1); fd$group <- items$group[match(fd$code, items$code)]
fd$dgl <- items$dgl[match(fd$code, items$code)]; fd$vita <- items$vita[match(fd$code, items$code)]
fd$consumed <- ifelse(num(fd$s1aq2) == 1, 1, 0)
H_gm <- hh_from_items(fd, as.character(fd$hid), fd$group, fd$consumed, num(fd$s1aq5), num(fd$s1aq7a), fd$dgl, fd$vita)
# Engel share from purchases: 7-day food x 52 over food + annualised non-food purchases (7d x 52, 1m x 12, 3m x 4, 12m x 1)
nf <- list(c("2A-Nonfood last 7 days", "s2aq3", 52), c("2B-Nonfood last 1 month", "s2bq3", 12), c("2C-Nonfood last 3 months", "s2cq3", 4), c("2D-Nonfood last 12 months", "s2dq3", 1))
nonfood <- bind_rows(lapply(nf, function(z) { d <- haven::read_dta(file.path(gm_dir, sprintf("Part B Section %s.dta", z[1]))); v <- num(d[[z[2]]]); v[!is.finite(v) | v < 0] <- 0
  data.frame(id = as.character(d$hid), nonfood = v * as.numeric(z[3])) })) |> group_by(id) |> summarise(nonfood_annual = sum(nonfood), .groups = "drop")
geo <- fd |> group_by(id = as.character(hid)) |> summarise(eanum = first(eanum), lga = first(as.character(haven::as_factor(lga))), district = first(as.character(haven::as_factor(district))), .groups = "drop")
wts <- haven::read_dta(file.path(gm_dir, "Household adjusted_weight.dta")); geo$w <- num(wts$hhweight)[match(num(geo$eanum), num(wts$eanum))]
H_gm <- left_join(geo, H_gm, by = "id") |> left_join(nonfood, by = "id")
H_gm$food_annual <- H_gm$food_purch_value * 52
H_gm$food_share <- ifelse(is.finite(H_gm$food_annual) & (H_gm$food_annual + coalesce(H_gm$nonfood_annual, 0)) > 0,
                          H_gm$food_annual / (H_gm$food_annual + coalesce(H_gm$nonfood_annual, 0)), NA_real_)
H_gm$log_cons_pae <- NA_real_; H_gm$own_prod_share <- H_gm$own_prod_item_share
# district names -> GADM Admin-2 (spine names), within the division
sp_gm <- spine[spine$country == "Gambia", ]
# GBoS 2015/16 districts that post-date or subdivide the GADM 4.1 units: the
# Kanifing municipality wards, and the districts carved from Upper / Lower
# Baddibu (North Bank) and from Fulladu East (Upper River), the same aliases
# the ESPEN build used
gm_alias <- c(bakau = "Kanifing", oldjeshwang = "Kanifing", newjeshwang = "Kanifing", serekundacentral = "Kanifing", serekundaeast = "Kanifing", serekundawest = "Kanifing",
              illiasa = "Upper Baddibu", sabachsanjar = "Upper Baddibu", sabachsanjal = "Upper Baddibu",   # Upper Baddibu was split into Illiasa and Sabach Sanjal
              basse = "Fulladu East", jimara = "Fulladu East", tumana = "Fulladu East",                 # Fulladu East -> Basse, Jimara, Tumana
              wulieast = "Wuli", wuliwest = "Wuli", lowerfulladuwest = "Fulladu West", lowerfuladuwest = "Fulladu West", upperfuladuwest = "Fulladu West", kombonorth = "Kombo Saint Mary", kombonorthstmary = "Kombo Saint Mary")
map_gm <- data.frame(district = unique(H_gm$district), stringsAsFactors = FALSE)
map_gm$target <- ifelse(kk(map_gm$district) %in% names(gm_alias), unname(gm_alias[kk(map_gm$district)]), map_gm$district)
dm <- stringdist::stringdistmatrix(kk(map_gm$target), kk(sp_gm$Admin2), method = "jw", p = 0.1)
best <- apply(dm, 1, which.min); dist <- apply(dm, 1, min)
map_gm$Admin2 <- ifelse(dist <= 0.15, sp_gm$Admin2[best], NA_character_); map_gm$jw <- round(dist, 3)
map_gm$Admin1 <- sp_gm$Admin1[match(map_gm$Admin2, sp_gm$Admin2)]
map_gm$n_hh <- as.integer(table(H_gm$district)[map_gm$district])
write.csv(map_gm[order(map_gm$Admin1, map_gm$Admin2), ], "metadata/hces_gambia_district_map.csv", row.names = FALSE)
cat("  district -> GADM Admin-2 matches:", sum(!is.na(map_gm$Admin2)), "of", nrow(map_gm), "| unmatched:", paste(map_gm$district[is.na(map_gm$Admin2)], collapse = ", "), "\n")
print(map_gm[order(map_gm$Admin1, map_gm$Admin2), c("district", "Admin2", "jw", "n_hh")], row.names = FALSE)
H_gm <- left_join(H_gm, map_gm[, c("district", "Admin1", "Admin2")], by = "district")
H_gm$country <- "Gambia"
cat(sprintf("  %d households, weighted food (purchase) share %.2f, HDDS %.1f, any ASF %.2f\n", nrow(H_gm), wm(H_gm$food_share, H_gm$w), wm(H_gm$hdds, H_gm$w), wm(H_gm$any_asf, H_gm$w)))

# ── SIERRA LEONE SLIHS 2011 ──────────────────────────────────────────────────
# HC-02 (2026-09-16): SLIHS 2011 replaces SLIHS 2018. It sits two years before
# the 2013 micronutrient survey instead of five after, and its diary item codes
# are labelled (166 purchase items, 106 own-consumption items), so Sierra Leone
# gets the food-group indicators the 2018 release could not give. The diary
# runs over five visits about five days apart (visits 2-6 record purchases
# since the previous visit); own consumption is recorded per visit in the same
# module. Definitions follow The Gambia: food share from purchases (diary food
# over diary food + diary non-food + the 12-month non-food scaled to the diary
# window), own production by ITEM (share of consumed items with any own-
# produced quantity in the diary), consumption level from annualised purchases
# per capita (no aggregate is released for 2011), centred nationally.
cat("\n[Sierra Leone SLIHS 2011]\n")
sl_dir <- file.path(RA, "LSMS/SLE_2011/SLE_2011_SLIHS_v01_M_Stata")
fp <- haven::read_dta(file.path(sl_dir, "Section 13 Part A Food purchases.dta"),
                      col_select = c("district", "sector", "slihseacode", "hid", "s13aitem", "s13aq2", "s13aq3", "s13aq4", "s13aq5", "s13aq6", "s13aq7")) |> as.data.frame()
il <- lab(fp$s13aitem); items <- data.frame(code = as.integer(names(il)), item = unname(il), stringsAsFactors = FALSE)
items$group <- mapply(function(l, c) classify_item(l, c, NA_character_), items$item, items$code)
items$dgl <- dgl_item(items$item); items$vita <- vita_item(items$item)
write.csv(items, "metadata/hces_food_groups_SierraLeone.csv", row.names = FALSE)
cat(sprintf("  %d purchase items classified: %s\n", nrow(items), paste(names(table(items$group)), table(items$group), collapse = " ")))
vis <- sapply(c("s13aq3", "s13aq4", "s13aq5", "s13aq6", "s13aq7"), function(v) { z <- num(fp[[v]]); z[!is.finite(z) | z < 0] <- 0; z })
fp$value <- rowSums(vis); fp$code <- num(fp$s13aitem); fp$id <- as.character(fp$hid)
fp$group <- items$group[match(fp$code, items$code)]; fp$dgl <- items$dgl[match(fp$code, items$code)]; fp$vita <- items$vita[match(fp$code, items$code)]
# own consumption during the diary: any home-produced quantity at visits 2-6 (module 12H, same item codes)
oc <- haven::read_dta(file.path(sl_dir, "Section 12 Part H Own food consumption.dta"),
                      col_select = c("hid", "s12hitem", "s12hq5", "s12hq6", "s12hq7", "s12hq8", "s12hq9")) |> as.data.frame()
oq <- sapply(c("s12hq5", "s12hq6", "s12hq7", "s12hq8", "s12hq9"), function(v) { z <- num(oc[[v]]); z[!is.finite(z) | z < 0] <- 0; z })
own <- data.frame(id = as.character(oc$hid), code = num(oc$s12hitem), own_qty = rowSums(oq)) |> filter(own_qty > 0) |> distinct(id, code, .keep_all = TRUE)
fp <- left_join(fp, own, by = c("id", "code")); fp$own_qty[is.na(fp$own_qty)] <- 0
fp$consumed <- ifelse(num(fp$s13aq2) == 1 | fp$own_qty > 0, 1, 0)
H_sl <- hh_from_items(fp, fp$id, fp$group, fp$consumed, fp$value, fp$own_qty, fp$dgl, fp$vita)
# diary window per household from the visit dates (visits 2-6, about five days apart; the first interval is
# the five days before visit 2); 25 days where the dates are missing or inconsistent
dt <- haven::read_dta(file.path(sl_dir, "Section 13 Part A Food purchases-interview dates.dta")) |> as.data.frame()
mk <- function(dd, mm, yy) as.Date(sprintf("%04d-%02d-%02d", 2000L + as.integer(num(yy)), as.integer(num(mm)), as.integer(num(dd))), optional = TRUE)
d2 <- mk(dt$day2_dd, dt$day2_mm, dt$day2_yr); d6 <- mk(dt$day6_dd, dt$day6_mm, dt$day6_yr)
span <- as.numeric(d6 - d2); span[!is.finite(span) | span < 10 | span > 60] <- 20
diary <- data.frame(id = as.character(dt$hid), diary_days = span * 5 / 4, stringsAsFactors = FALSE) |> distinct(id, .keep_all = TRUE)
# non-food: frequent purchases over the diary (module 13B1) and infrequent purchases over 12 months (13B2)
b1 <- haven::read_dta(file.path(sl_dir, "Section 13 Part B-B1 NonFood freq purchases.dta"), col_select = c("hid", "s13b1q3", "s13b1q4", "s13b1q5", "s13b1q6", "s13b1q7")) |> as.data.frame()
nf1 <- data.frame(id = as.character(b1$hid), v = rowSums(sapply(c("s13b1q3", "s13b1q4", "s13b1q5", "s13b1q6", "s13b1q7"), function(v) { z <- num(b1[[v]]); z[!is.finite(z) | z < 0] <- 0; z }))) |>
  group_by(id) |> summarise(nonfood_diary = sum(v), .groups = "drop")
b2 <- haven::read_dta(file.path(sl_dir, "Section 13 Part B-B2 NonFood infreq purchases.dta"), col_select = c("hid", "s13b2q4")) |> as.data.frame()
nf2 <- data.frame(id = as.character(b2$hid), v = { z <- num(b2$s13b2q4); z[!is.finite(z) | z < 0] <- 0; z }) |> group_by(id) |> summarise(nonfood_12m = sum(v), .groups = "drop")
# household size, EA weight, district
ro <- haven::read_dta(file.path(sl_dir, "Section 1-5 Individual-level.dta"), col_select = c("hid")) |> as.data.frame()
hhsize <- data.frame(id = as.character(ro$hid)) |> count(id, name = "hhsize")
wts <- haven::read_dta(file.path(sl_dir, "Household Cluster Weights.dta")) |> as.data.frame()
geo <- fp |> group_by(id) |> summarise(ea = first(as.character(slihseacode)), district_code = first(num(district)), urban = first(num(sector)), .groups = "drop")
geo$w <- num(wts$weight)[match(geo$ea, as.character(wts$slihseacode))]
dl <- lab(fp$district); geo$district <- unname(dl[as.character(geo$district_code)])
H_sl <- geo |> left_join(H_sl, by = "id") |> left_join(diary, by = "id") |> left_join(nf1, by = "id") |> left_join(nf2, by = "id") |> left_join(hhsize, by = "id")
H_sl$diary_days[is.na(H_sl$diary_days)] <- 25
nf <- coalesce(H_sl$nonfood_diary, 0) + coalesce(H_sl$nonfood_12m, 0) * H_sl$diary_days / 365
H_sl$food_share <- ifelse(is.finite(H_sl$food_purch_value) & (H_sl$food_purch_value + nf) > 0, H_sl$food_purch_value / (H_sl$food_purch_value + nf), NA_real_)
H_sl$own_prod_share <- H_sl$own_prod_item_share
ann <- (H_sl$food_purch_value + nf) * 365 / H_sl$diary_days
H_sl$log_cons_pae <- ifelse(is.finite(ann) & ann > 0 & is.finite(H_sl$hhsize) & H_sl$hhsize > 0, log(ann / H_sl$hhsize), NA_real_)
# 2011 district names -> GADM 4.1 (Western Area split into Rural and Urban in both)
sl_alias <- c(tonkilili = "Tonkolili", westernother = "Western Rural", westernurban = "Western Urban")
sp_sl <- spine[spine$country == "SierraLeone", ]
H_sl$Admin2 <- ifelse(kk(H_sl$district) %in% names(sl_alias), unname(sl_alias[kk(H_sl$district)]), H_sl$district)
H_sl$Admin2 <- sp_sl$Admin2[match(kk(H_sl$Admin2), kk(sp_sl$Admin2))]; H_sl$Admin1 <- sp_sl$Admin1[match(H_sl$Admin2, sp_sl$Admin2)]
stopifnot(all(!is.na(H_sl$Admin2)))
H_sl$country <- "SierraLeone"
cat(sprintf("  %d households in %d districts (diary %.0f days median); weighted food (purchase) share %.2f, own-production (item) share %.2f, HDDS %.1f, any ASF %.2f\n",
            nrow(H_sl), n_distinct(H_sl$Admin2), median(H_sl$diary_days), wm(H_sl$food_share, H_sl$w), wm(H_sl$own_prod_share, H_sl$w), wm(H_sl$hdds, H_sl$w), wm(H_sl$any_asf, H_sl$w)))

# ── GHANA GLSS7 2016/17 ──────────────────────────────────────────────────────
cat("\n[Ghana GLSS7 2016/17]\n")
g7 <- haven::read_dta("data/LSMS/g7aggregates_hhlevel.dta") |> as.data.frame()
H_gh <- data.frame(country = "Ghana", id = as.character(seq_len(nrow(g7))), w = num(g7$WTA_S), region = as.character(haven::as_factor(g7$region)),
                   food_share = ifelse(num(g7$HHEXP_N) > 0, num(g7$TOTFDAL) / num(g7$HHEXP_N), NA_real_),
                   own_prod_share = ifelse(num(g7$TOTFOOD) > 0, num(g7$FD_P) / num(g7$TOTFOOD), NA_real_),
                   log_cons_pae = log(pmax(num(g7$padq_hh_R), 1)), stringsAsFactors = FALSE)
XW <- read.csv("data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv", stringsAsFactors = FALSE)
sp_gh <- spine[spine$country == "Ghana", ]
r10 <- XW$adm1_paper_10[match(kk(sp_gh$Admin1), kk(XW$admin1_16))]; r10[is.na(r10)] <- sp_gh$Admin1[is.na(r10)]
cat(sprintf("  %d households in %d regions; weighted food share %.2f, own-production %.2f\n", nrow(H_gh), n_distinct(H_gh$region), wm(H_gh$food_share, H_gh$w), wm(H_gh$own_prod_share, H_gh$w)))

# ── aggregate to the spine ───────────────────────────────────────────────────
# consumption per adult equivalent is in each country's own currency and year,
# so it is centred on the household-weighted national mean: what enters the
# pooled set is log consumption RELATIVE to the country, not a currency level
centre_cons <- function(H) { if (any(is.finite(H$log_cons_pae))) H$log_cons_pae_rel <- H$log_cons_pae - wm(H$log_cons_pae, H$w) else H$log_cons_pae_rel <- NA_real_; H }
H_mw <- centre_cons(H_mw); H_gm <- centre_cons(H_gm); H_sl <- centre_cons(H_sl); H_gh <- centre_cons(H_gh)
IND <- c("food_share", "own_prod_share", "log_cons_pae_rel", "hdds", "any_asf", "any_fish", "any_meat", "any_eggs", "any_dairy", "any_pulses_nuts",
         "any_fruit", "any_veg", "any_dgl", "any_vita_fv", "asf_purchase_share")
agg <- function(H, keys) {
  for (v in IND) if (!v %in% names(H)) H[[v]] <- NA_real_
  H |> filter(!is.na(.data[[keys[length(keys)]]])) |> group_by(across(all_of(keys))) |>
    summarise(n_hh = n(), across(all_of(IND), ~ wm(as.numeric(.x), w)), .groups = "drop") |> as.data.frame()
}
# Malawi: TA where >= MIN_HH households, else the district (Admin1) mean
mw_ta <- agg(H_mw, c("Admin1", "Admin2")); mw_d1 <- agg(H_mw, "Admin1")
sp_mw <- spine[spine$country == "Malawi", ]
B_mw <- left_join(sp_mw, mw_ta, by = c("Admin1", "Admin2"))
use_d1 <- is.na(B_mw$n_hh) | B_mw$n_hh < MIN_HH
for (v in IND) B_mw[[v]][use_d1] <- mw_d1[[v]][match(B_mw$Admin1[use_d1], mw_d1$Admin1)]
B_mw$hces_level <- ifelse(use_d1, "admin1_broadcast", "admin2"); B_mw$n_hh[use_d1] <- mw_d1$n_hh[match(B_mw$Admin1[use_d1], mw_d1$Admin1)]
cat(sprintf("\nMalawi: %d TAs from their own households, %d from the district mean\n", sum(!use_d1), sum(use_d1)))
# Gambia: district (Admin-2) means; fall back to the division (Admin1) mean
gm_a2 <- agg(H_gm, c("Admin1", "Admin2")); gm_a1 <- agg(H_gm, "Admin1")
B_gm <- left_join(sp_gm, gm_a2, by = c("Admin1", "Admin2")); use_a1 <- is.na(B_gm$n_hh) | B_gm$n_hh < MIN_HH
for (v in IND) B_gm[[v]][use_a1] <- gm_a1[[v]][match(B_gm$Admin1[use_a1], gm_a1$Admin1)]
B_gm$hces_level <- ifelse(use_a1, "admin1_broadcast", "admin2"); B_gm$n_hh[use_a1] <- gm_a1$n_hh[match(B_gm$Admin1[use_a1], gm_a1$Admin1)]
cat(sprintf("Gambia: %d districts from their own households, %d from the division mean\n", sum(!use_a1), sum(use_a1)))
# Sierra Leone: districts
sl_a2 <- agg(H_sl, c("Admin1", "Admin2")); sp_sl <- spine[spine$country == "SierraLeone", ]
B_sl <- left_join(sp_sl, sl_a2, by = c("Admin1", "Admin2")); B_sl$hces_level <- "admin2"
cat(sprintf("Sierra Leone: %d of %d districts\n", sum(!is.na(B_sl$n_hh)), nrow(B_sl)))
# Ghana: regions broadcast
gh_r <- agg(H_gh, "region"); B_gh <- sp_gh; j <- match(kk(r10), kk(gh_r$region))
for (v in c("n_hh", IND)) B_gh[[v]] <- gh_r[[v]][j]
B_gh$hces_level <- "admin1_broadcast"
cat(sprintf("Ghana: %d of %d districts mapped to a GLSS7 region\n", sum(!is.na(j)), nrow(B_gh)))

OUT <- bind_rows(B_mw, B_gm, B_sl, B_gh)
names(OUT)[names(OUT) %in% IND] <- paste0("hces_", IND)
OUT <- OUT |> rename(hces_n_hh = n_hh)
OUT <- OUT[match(paste(spine$country, spine$Admin1, spine$Admin2), paste(OUT$country, OUT$Admin1, OUT$Admin2)), ]
write.csv(OUT, file.path(HDIR, "predictors_admin2_hces.csv"), row.names = FALSE)
for (H in list(H_mw, H_gm, H_sl, H_gh)) write.csv(H, file.path(HDIR, sprintf("hces_household_%s.csv", H$country[1])), row.names = FALSE)

# ── metadata ─────────────────────────────────────────────────────────────────
basis <- "Malawi IHS4 2016-17 (7-day recall, TA by EA point-in-polygon, district mean where a TA has < 15 households); The Gambia IHS 2015/16 (7-day recall, district by name; food share is PURCHASE-based: 7-day food purchases x 52 over food + annualised non-food purchases, no consumption aggregate exists); Sierra Leone SLIHS 2011 (diary purchases over five visits, own consumption per visit, district by code; food share PURCHASE-based like The Gambia; consumption level from annualised purchases per capita); clusters geolocated with the MICS 2017 GPS they share; diary item codes unlabelled so no recall indicators); Ghana GLSS7 2016/17 (expenditure aggregates, 10 regions broadcast through the 16-to-10 crosswalk)."
md <- data.frame(column = c(paste0("hces_", IND), "hces_n_hh"), source = "HCES microdata (IHS4, IHS 2015/16, SLIHS 2011, GLSS7)",
                 domain = "Household diet and consumption (HCES)", subnational = TRUE, stringsAsFactors = FALSE)
md$assumption <- paste0(c(
  "Food (purchased + own production + gifts) share of total household consumption",
  "Own-production share of food: by VALUE for Ghana, by ITEM (share of consumed items with any own-production quantity) for Malawi, The Gambia and Sierra Leone",
  "log consumption per adult equivalent (per capita for Malawi, the aggregate's basis; per capita annualised PURCHASES for Sierra Leone, no aggregate), centred on the household-weighted national mean because the aggregates are in different currencies; NA for The Gambia (no aggregate)",
  "Mean count of the 12 HDDS food groups consumed in the recall (7-day: Malawi, The Gambia; ~25-day diary: Sierra Leone)",
  "Share of households consuming any meat, fish, eggs or dairy in the recall (Malawi, The Gambia, Sierra Leone)",
  "Share consuming any fish (recall countries)", "Share consuming any meat/poultry (recall countries)", "Share consuming any eggs (recall countries)",
  "Share consuming any milk/dairy (recall countries)", "Share consuming any pulses, nuts or seeds (recall countries)",
  "Share consuming any fruit (recall countries)", "Share consuming any vegetables (recall countries)",
  "Share consuming any dark green leafy vegetables (recall countries)", "Share consuming any vitamin-A-rich fruit or vegetable (recall countries)",
  "Animal-source share of food PURCHASE value (recall countries)", "Households behind the area estimate"), ". ", basis)
md$n_countries <- vapply(md$column, function(v) sum(tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))), na.rm = TRUE), 0L)
md$countries <- vapply(md$column, function(v) { s <- tapply(OUT[[v]], OUT$country, function(z) any(is.finite(z))); paste(names(s)[!is.na(s) & s], collapse = ";") }, "")
md$completeness <- round(vapply(md$column, function(v) mean(is.finite(OUT[[v]])), 0), 3)
write.csv(md, file.path(HDIR, "predictors_admin2_hces_metadata.csv"), row.names = FALSE)
cat("\n=== HCES block ===\n"); print(md[, c("column", "n_countries", "countries", "completeness")], row.names = FALSE)
cat("\nweighted country means:\n"); print(OUT |> group_by(country) |> summarise(across(c(hces_food_share, hces_own_prod_share, hces_hdds, hces_any_asf, hces_any_dgl), ~ round(mean(.x, na.rm = TRUE), 2))), n = 10)
cat("DONE\n")
