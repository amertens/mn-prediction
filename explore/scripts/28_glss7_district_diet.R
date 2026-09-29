# =============================================================================
# SUPERSEDED by explore/scripts/29_glss7_block_classify.R.
#
# This script's district table is WRONG and is kept only as the "before" case
# that 29 compares against. GLSS7 labels section-9b items by BRAND, and the
# project's classify_item() - written against the generic labels used by
# Malawi, The Gambia and Sierra Leone - sent 270 of 484 real foods to "misc".
# Script 29 fixes it using GLSS7's own item-code block structure. Use
# explore/out/29_glss7_district_diet.csv, not 28's.
# =============================================================================
# =============================================================================
# explore/scripts/28_glss7_district_diet.R   [HC-06, task 3]
#
# GLSS7 food-group indicators at Ghana district level, using the project's OWN
# classify_item() / hh_from_items() so Ghana stays pooled-comparable with
# Malawi, The Gambia and Sierra Leone.
#
# Input is the streamed section-9b table from script 25 (528,678 household x
# item rows, 13,924 households, 484 items).
#
# THE ONE THING MISSING is the district NAME. GLSS7 labels every geographic
# variable - REGION, loc2, loc5, loc7, ez - EXCEPT `district`, which is a
# disclosure control in the public release, not a lost file. The codes are
# region*100 + the district's index in the official GSS ordering within region.
# The 216 districts are the set created 28 June 2012 (NOT a 2010 census set,
# which had 170); per-region counts 22/20/16/25/26/30/27/26/13/11 = 216. The
# ordering is not alphabetical - see HC-07. So this writes the
# finished district-level table keyed by CODE, and the join to the project's
# 260-district spine waits on one external artefact: the GSS 2010 district
# list in official per-region order.
#
# It also measures what that artefact is worth: the share of district-level
# variance in each indicator that is WITHIN region, i.e. exactly the signal a
# region-only broadcast throws away.
#
#   Rscript explore/scripts/28_glss7_district_diet.R
# -> explore/out/28_glss7_district_diet.csv   214 district codes
#    explore/out/28_glss7_region_diet.csv     10 regions
#    explore/out/28_glss7_variance_within_region.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(haven)})
ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction/"

# ── the project's own classifier, taken verbatim from the production builder ──
want <- c("classify_item", "dgl_item", "vita_item", "asf_groups", "hh_from_items")
ex   <- parse(file.path(ROOT, "scripts/covariates/build_hces_diet_block.R"))
got  <- character()
for (e in ex) {
  if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% want) {
    eval(e, envir = globalenv()); got <- c(got, as.character(e[[2]]))
  }
}
stopifnot(setequal(got, want))
message("reused from build_hces_diet_block.R: ", paste(got, collapse = ", "))

# ── items ────────────────────────────────────────────────────────────────────
IT  <- read.csv(file.path(ROOT, "explore/out/25_glss7_hh_item.csv"),
                stringsAsFactors = FALSE)
LAB <- read.csv(file.path(ROOT, "explore/out/25_glss7_item_labels.csv"),
                stringsAsFactors = FALSE)
IT  <- left_join(IT, LAB, by = "code")
message(nrow(IT), " household-item rows | ", dplyr::n_distinct(IT$hid),
        " households | ", dplyr::n_distinct(IT$code), " items | ",
        sum(is.na(IT$item)), " rows with no label")

LAB$group <- vapply(LAB$item, classify_item, "", USE.NAMES = FALSE)
LAB$dgl   <- dgl_item(LAB$item)
LAB$vita  <- vita_item(LAB$item)
cat("\n== items per food group ==\n"); print(table(LAB$group))
cat("\ndark green leafy items:", sum(LAB$dgl),
    "| vitamin-A-rich fruit/veg items:", sum(LAB$vita), "\n")
IT <- left_join(IT, LAB[, c("code", "group", "dgl", "vita")], by = "code")
IT$group[is.na(IT$group)] <- "misc"

# section 9b is a purchase diary with no own-production quantity column, so
# own_prod_item_share is not identified here and is dropped downstream
HH <- hh_from_items(IT, id = IT$hid, group = IT$group, consumed = 1L,
                    value_purch = IT$value, own_qty = NA_real_,
                    dgl = IT$dgl, vita = IT$vita)
HH$own_prod_item_share <- NULL
message("\nhousehold indicators: ", nrow(HH), " households, ",
        ncol(HH) - 1, " columns")

# ── geography and weights ────────────────────────────────────────────────────
S0 <- haven::read_dta(file.path(ROOT, "data/LSMS/GHA_2017/g7sec0.dta")) |>
  transmute(id = as.character(hid), clust, region = as.integer(region),
            district = as.integer(district), w = as.numeric(WTA_S),
            hhsize = as.numeric(hhsize))
H <- inner_join(HH, S0, by = "id")
message("joined to geography: ", nrow(H), " of ", nrow(HH), " households | ",
        dplyr::n_distinct(H$district), " district codes | ",
        dplyr::n_distinct(H$clust), " EAs")

IND <- setdiff(names(HH), "id")
wmean <- function(x, w) {
  x <- as.numeric(x); ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok])
}
agg <- function(d, key) {
  d |> group_by(across(all_of(key))) |>
    summarise(n_hh = n(), n_ea = dplyr::n_distinct(clust),
              across(all_of(IND), ~ wmean(.x, w)), .groups = "drop")
}
DIS <- agg(H, c("region", "district"))
REG <- agg(H, "region")

cat("\n== households per district code ==\n")
print(summary(DIS$n_hh)); cat("EAs per district code: median ",
    stats::median(DIS$n_ea), " (min ", min(DIS$n_ea), ")\n", sep = "")

# ── what the missing crosswalk is worth ──────────────────────────────────────
# share of district-level variance that lies WITHIN region: the part a
# region-only broadcast (what Ghana has today) cannot represent
VW <- do.call(rbind, lapply(IND, function(v) {
  x <- DIS[[v]]; ok <- is.finite(x)
  if (sum(ok) < 20 || stats::var(x[ok]) == 0) return(NULL)
  fit <- stats::lm(x[ok] ~ factor(DIS$region[ok]))
  data.frame(indicator = v, sd_district = stats::sd(x[ok]),
             within_region_share = 1 - summary(fit)$r.squared)
}))
VW <- VW[order(-VW$within_region_share), ]
cat("\n== share of district variance that is WITHIN region ==\n")
print(VW, row.names = FALSE, digits = 3)
cat("\nmedian within-region share:",
    signif(stats::median(VW$within_region_share), 3), "\n")

dir.create(file.path(ROOT, "explore/out"), showWarnings = FALSE)
write.csv(DIS, file.path(ROOT, "explore/out/28_glss7_district_diet.csv"), row.names = FALSE)
write.csv(REG, file.path(ROOT, "explore/out/28_glss7_region_diet.csv"), row.names = FALSE)
write.csv(VW,  file.path(ROOT, "explore/out/28_glss7_variance_within_region.csv"), row.names = FALSE)
write.csv(LAB[, c("code", "item", "group", "dgl", "vita")],
          file.path(ROOT, "explore/out/28_glss7_item_classified.csv"), row.names = FALSE)
cat("\nwrote district (", nrow(DIS), "), region (", nrow(REG),
    "), variance and item-classification tables\n", sep = "")
