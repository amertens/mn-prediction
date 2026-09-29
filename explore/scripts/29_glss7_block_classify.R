# =============================================================================
# explore/scripts/29_glss7_block_classify.R   [HC-06, task 3 - corrected]
#
# Script 28 built the Ghana district diet table but got it WRONG, and this
# fixes it. GLSS7 labels section-9b items by BRAND - "Geisha", "Titus",
# "John West" are canned fish; "Peak", "Nido", "Cowbell" are milk; "Gino" is
# tomato paste - while the project's classify_item() was written against the
# generic labels used by Malawi, The Gambia and Sierra Leone. It therefore
# dumped 270 of 484 real foods into "misc", understating every indicator.
#
# The fix uses GLSS7's own structure rather than a new regex. Item codes run in
# contiguous blocks separated by numeric gaps, one food per block, each ending
# in "Other <food>" (10-23 imported rice, 58-62 bread, 95-101 beef, ...). So:
#
#   1. segment items into blocks wherever the code sequence gaps
#   2. take each block's food group by majority vote of the items the project
#      classifier already resolves on their own generic labels
#   3. re-run the PROJECT's classify_item(label, code, block) passing that
#      block - its final line is `if (!is.na(block)) return(block)`, which is
#      exactly the fallback this is for
#
# So the brands inherit their block's food, and no project logic is edited;
# Ghana stays pooled-comparable with the other three countries.
#
# Blocks with no self-resolving member stay "misc" and are reported, with the
# expenditure share they carry, rather than guessed at.
#
# ONE GUARD, added after auditing the first run. A block vote propagates the
# classifier's false positives as well as its hits: "Football game" matches
# \bgame\b and made a recreation block meat, "Duck soap" made a soap block
# meat, "Insecticides" matches \binsect and made a pesticide block meat, and
# "Other fruit drink" made a juice block fruit. Every one of these sits beyond
# block 31. GLSS7's diary is ordered food first, then condiments (block 32,
# Maggi cube), salt, coffee, beverages, spirits, tobacco and finally non-food,
# so block 31 (spices, codes 523-529) is the last food block. Voting is
# therefore restricted to blocks 1-31 - one structural cut rather than a
# growing veto list, and it coincides with the project's own convention that
# condiments, stimulants and beverages are "misc" whatever plant they come
# from. Items beyond block 31 keep whatever classify_item() gives them on
# their own label, so "Cooked rice and stew" is still cereals.
LAST_FOOD_BLOCK <- 31L
#
#   Rscript explore/scripts/29_glss7_block_classify.R
# -> explore/out/29_glss7_item_classified.csv   485 items, block-resolved
#    explore/out/29_glss7_district_diet.csv     214 district codes
#    explore/out/29_glss7_region_diet.csv       10 regions
#    explore/out/29_glss7_variance_within_region.csv
#    explore/out/29_glss7_blocks.csv            the audit trail
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(haven)})
ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction/"

want <- c("classify_item", "dgl_item", "vita_item", "asf_groups", "hh_from_items")
ex   <- parse(file.path(ROOT, "scripts/covariates/build_hces_diet_block.R"))
for (e in ex) if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
                  is.name(e[[2]]) && as.character(e[[2]]) %in% want) eval(e, envir = globalenv())
message("reused verbatim: ", paste(want, collapse = ", "))

LAB <- read.csv(file.path(ROOT, "explore/out/25_glss7_item_labels.csv"),
                stringsAsFactors = FALSE)
LAB$item <- trimws(LAB$item)
LAB <- LAB[order(LAB$code), ]

# -- 1. segment into blocks on gaps in the code sequence ---------------------
LAB$block_id <- cumsum(c(1L, as.integer(diff(LAB$code) > 1L)))
# -- 2. each block's group = majority vote of the self-resolving items -------
LAB$g0 <- vapply(LAB$item, classify_item, "", USE.NAMES = FALSE)
blk <- LAB |> group_by(block_id) |>
  summarise(n_items = n(), n_self = sum(g0 != "misc"),
            group = { t <- table(g0[g0 != "misc"])
                      if (!length(t) || block_id[1] > LAST_FOOD_BLOCK) NA_character_
                      else names(t)[which.max(t)] },
            evidence = paste(utils::head(item[g0 != "misc"], 3), collapse = " | "),
            first_item = item[1], .groups = "drop")
message("blocks: ", nrow(blk), " | resolved ", sum(!is.na(blk$group)),
        " | unresolved ", sum(is.na(blk$group)))

# -- 3. re-run the project classifier WITH the block ------------------------
LAB <- left_join(LAB, blk[, c("block_id", "group")], by = "block_id")
LAB$group_new <- mapply(classify_item, LAB$item, LAB$code, LAB$group,
                        USE.NAMES = FALSE)
LAB$dgl  <- dgl_item(LAB$item)
LAB$vita <- vita_item(LAB$item)

lev <- sort(unique(c(LAB$g0, LAB$group_new)))
cat("\n== items per food group: before -> after ==\n")
print(cbind(before = table(factor(LAB$g0, levels = lev)),
            after  = table(factor(LAB$group_new, levels = lev))))
cat("\nbrands rescued from misc by their block:",
    sum(LAB$g0 == "misc" & LAB$group_new != "misc"), "\n")

cat("\n== the ten blocks that rescued the most items ==\n")
res <- LAB |> filter(g0 == "misc", group_new != "misc") |>
  count(block_id, group_new, name = "rescued") |>
  left_join(blk[, c("block_id", "n_items", "evidence")], by = "block_id") |>
  arrange(-rescued)
print(utils::head(as.data.frame(res), 10), row.names = FALSE)

cat("\n== food blocks left unresolved (no self-classifying member) ==\n")
print(as.data.frame(blk[is.na(blk$group) & blk$block_id <= LAST_FOOD_BLOCK,
                        c("block_id", "n_items", "n_self", "first_item")]),
      row.names = FALSE)
cat("\nblocks beyond", LAST_FOOD_BLOCK, "(condiments, beverages, non-food) not voted:",
    sum(blk$block_id > LAST_FOOD_BLOCK), "\n")

# -- rebuild the household and district tables -------------------------------
IT <- read.csv(file.path(ROOT, "explore/out/25_glss7_hh_item.csv"), stringsAsFactors = FALSE)
IT <- left_join(IT, LAB[, c("code", "group_new", "dgl", "vita")], by = "code")
IT$group_new[is.na(IT$group_new)] <- "misc"
unres <- sum(IT$value[IT$group_new == "misc"], na.rm = TRUE) / sum(IT$value, na.rm = TRUE)
cat(sprintf("\nshare of all diary expenditure still in misc: %.1f%%\n", 100 * unres))

HH <- hh_from_items(IT, id = IT$hid, group = IT$group_new, consumed = 1L,
                    value_purch = IT$value, own_qty = NA_real_,
                    dgl = IT$dgl, vita = IT$vita)
HH$own_prod_item_share <- NULL

S0 <- haven::read_dta(file.path(ROOT, "data/LSMS/GHA_2017/g7sec0.dta")) |>
  transmute(id = as.character(hid), clust, region = as.integer(region),
            district = as.integer(district), w = as.numeric(WTA_S))
H <- inner_join(HH, S0, by = "id")
IND <- setdiff(names(HH), "id")
wmean <- function(x, w) { x <- as.numeric(x); ok <- is.finite(x) & is.finite(w) & w > 0
                          if (!any(ok)) NA_real_ else sum(x[ok] * w[ok]) / sum(w[ok]) }
agg <- function(d, key) d |> group_by(across(all_of(key))) |>
  summarise(n_hh = n(), n_ea = dplyr::n_distinct(clust),
            across(all_of(IND), ~ wmean(.x, w)), .groups = "drop")
DIS <- agg(H, c("region", "district")); REG <- agg(H, "region")

cat("\n== district means: script 28 (broken) vs 29 (block-resolved) ==\n")
OLD <- read.csv(file.path(ROOT, "explore/out/28_glss7_district_diet.csv"))
cmp <- do.call(rbind, lapply(intersect(IND, names(OLD)), function(v)
  data.frame(indicator = v, before = mean(OLD[[v]], na.rm = TRUE),
             after = mean(DIS[[v]], na.rm = TRUE))))
cmp$change <- cmp$after - cmp$before
print(cmp, row.names = FALSE, digits = 3)

VW <- do.call(rbind, lapply(IND, function(v) {
  x <- DIS[[v]]; ok <- is.finite(x)
  if (sum(ok) < 20 || stats::var(x[ok]) == 0) return(NULL)
  fit <- stats::lm(x[ok] ~ factor(DIS$region[ok]))
  data.frame(indicator = v, sd_district = stats::sd(x[ok]),
             within_region_share = 1 - summary(fit)$r.squared) }))
VW <- VW[order(-VW$within_region_share), ]
cat("\n== share of district variance that is WITHIN region ==\n")
print(VW, row.names = FALSE, digits = 3)
cat("\nmedian within-region share:", signif(stats::median(VW$within_region_share), 3), "\n")

O <- file.path(ROOT, "explore/out")
write.csv(LAB[, c("code", "item", "block_id", "g0", "group_new", "dgl", "vita")],
          file.path(O, "29_glss7_item_classified.csv"), row.names = FALSE)
write.csv(blk, file.path(O, "29_glss7_blocks.csv"), row.names = FALSE)
write.csv(DIS, file.path(O, "29_glss7_district_diet.csv"), row.names = FALSE)
write.csv(REG, file.path(O, "29_glss7_region_diet.csv"), row.names = FALSE)
write.csv(VW,  file.path(O, "29_glss7_variance_within_region.csv"), row.names = FALSE)
cat("\nwrote item, block, district (", nrow(DIS), "), region (", nrow(REG),
    ") and variance tables\n", sep = "")
