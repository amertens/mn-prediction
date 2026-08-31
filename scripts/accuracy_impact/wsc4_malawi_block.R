# =============================================================================
# scripts/accuracy_impact/wsc4_malawi_block.R
#
# WS-C4. Admit Malawi to the individual-level arms.
#
# Malawi's merged dataset carries ONE `gw_` column against 430 to 599 in the
# other three countries, so its questionnaire arm has never been a questionnaire
# arm: WS-C1 correctly refused to score it and recorded `not_computed`.
#
# The micronutrient survey was nested in the 2015-16 DHS. The README says to
# merge on MCLUSTER, MNUMBER and M01; MCLUSTER does not exist under that name in
# the merged file, but `gw_cnum` is it. Verified before writing anything:
#
#   (gw_cnum, mnumber, m01) -> (hv001, hv002, hvidx)   3097 of 3099 rows
#
# 105 of the DHS's 850 clusters, which is the MNS subsample.
#
#   Rscript scripts/accuracy_impact/wsc4_malawi_block.R
# -> data/IPD/Malawi/Malawi_dhs_questionnaire_block.rds   (sidecar, additive)
# -> results/tables/wsc4_malawi_block_manifest.csv
#
# ADDITIVE. The merged dataset is not touched. A sidecar plus
# attach_malawi_questionnaire() keeps the provenance visible: these columns come
# from a different instrument than the other three countries' `gw_` blocks, and
# a reader who forgets that will over-read the comparison.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
source(here("R", "data_prep.R"))

CACHE <- "C:/Users/andre/AppData/Local/andre/rdhs/Cache/datasets"
OUT   <- here("data", "IPD", "Malawi", "Malawi_dhs_questionnaire_block.rds")
MAN   <- here("results", "tables", "wsc4_malawi_block_manifest.csv")
MIN_COMPLETE <- 0.80   # same threshold WS-F7 settled on, applied WITHIN the
                       # population a recode covers, not across the whole file
MAX_LEVELS   <- 30L    # beyond this a code is an identifier, not an item
MIN_MATCH    <- 0.50   # roster coverage floor, reported separately from
                       # completeness because they fail for different reasons

zz <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

mns <- readRDS(here("data", "IPD", "Malawi", "Malawi_merged_dataset.rds"))
stopifnot(all(c("gw_cnum", "mnumber", "m01") %in% names(mns)))
mns_key <- data.frame(cnum = zz(mns$gw_cnum), hh = zz(mns$mnumber),
                      ln = zz(mns$m01), row = seq_len(nrow(mns)))
# The MNS holds two populations in one file. IR covers only the women and KR
# only the children, so a match rate or a completeness computed over all 3099
# rows judges those recodes against respondents they were never going to reach.
# The first version of this script did exactly that and dropped both blocks --
# which are the ones carrying dietary recall, vitamin A capsule receipt,
# deworming and iron in pregnancy, i.e. every nutrition-proximal item the
# household block does not have.
pop <- as.character(mns$population)
is_child <- grepl("child", pop, ignore.case = TRUE)
is_woman <- !is_child & !is.na(pop)
cat(sprintf("[wsc4] Malawi MNS: %d rows (%d children, %d women), %d clusters\n",
            nrow(mns), sum(is_child), sum(is_woman),
            length(unique(mns_key$cnum))))

# ---------------------------------------------------------------------------
# Which DHS columns are eligible at all.
#
# Three exclusions, in order of how badly each would corrupt the arm:
#   1. the measurement bands, via the guard added for exactly this reason
#   2. identifiers, weights and design variables, which are not questionnaire
#      responses and would let the model recover the district directly
#   3. free-text and date fields, which carry no comparable information
# ---------------------------------------------------------------------------
ID_PAT <- paste0("^(hhid|caseid|bidx|midx|h[vw]idx|hhidx|hv00[0-9]|hv10[0-5]|",
                 "v00[0-9]|v10[0-5]|b16|idx|awfact|hv02[12345]|v00[678]|v01[678])")
DESIGN_PAT <- "wgt|weight|psu|strata|stratum|_id$|line|^sdist|^shdist|^v024$|^hv024$"
DATE_PAT   <- "cmc|^v01[1-7]$|date"

eligible <- function(nm) {
  keep <- !grepl(ID_PAT, nm) & !grepl(DESIGN_PAT, nm, ignore.case = TRUE) &
    !grepl(DATE_PAT, nm, ignore.case = TRUE)
  # The measurement guard is applied to the PREFIXED name, because that is what
  # it is scoped to and what the pipeline will see.
  keep & !is_biomarker_column(paste0("gw_", nm))
}

# ---------------------------------------------------------------------------
# The four recodes, each with the key it joins on.
#
# KR is birth-level: `b16` is the child's line number in the household roster,
# which is what m01 is. A birth that is not in the roster (died, or living
# elsewhere) has b16 = 0 and cannot match, which is correct.
# ---------------------------------------------------------------------------
SPECS <- list(
  list(file = "MWHR7ADT", keys = c("hv001", "hv002"),          tag = "hh", pop = "all"),
  list(file = "MWPR7ADT", keys = c("hv001", "hv002", "hvidx"), tag = "pr", pop = "all"),
  list(file = "MWIR7ADT", keys = c("v001", "v002", "v003"),    tag = "ir", pop = "women"),
  list(file = "MWKR7ADT", keys = c("v001", "v002", "b16"),     tag = "kr", pop = "child")
)

blocks <- list(); manifest <- list()
for (sp in SPECS) {
  p <- file.path(CACHE, paste0(sp$file, ".rds"))
  if (!file.exists(p)) { cat(sprintf("[wsc4] %s ABSENT, skipped\n", sp$file)); next }
  d <- readRDS(p)
  if (!all(sp$keys %in% names(d))) {
    cat(sprintf("[wsc4] %s missing keys, skipped\n", sp$file)); next
  }
  k <- lapply(sp$keys, function(z) zz(d[[z]]))
  jk <- if (length(k) == 3) paste(k[[1]], k[[2]], k[[3]], sep = "_") else
    paste(k[[1]], k[[2]], sep = "_")
  mk <- if (length(k) == 3) paste(mns_key$cnum, mns_key$hh, mns_key$ln, sep = "_") else
    paste(mns_key$cnum, mns_key$hh, sep = "_")

  cand <- setdiff(names(d), sp$keys)
  cand <- cand[eligible(cand)]
  if (!length(cand)) next

  # Position of each MNS row in this recode. First match only: the household
  # recode has one row per household and the person recodes one per person, so
  # a duplicate here would be a data error rather than a legitimate expansion.
  pos <- match(mk, jk)
  # The rows this recode is ELIGIBLE to reach. Judging IR against all 3099
  # rows counts every child as a miss and drops the block at 26 percent.
  in_pop <- switch(sp$pop, child = is_child, women = is_woman,
                   rep(TRUE, nrow(mns_key)))
  hit <- !is.na(pos) & in_pop
  matched <- sum(hit)
  denom <- sum(in_pop)
  match_rate <- if (denom > 0) matched / denom else 0
  if (matched < MIN_MATCH * denom) {
    cat(sprintf("[wsc4] %-9s matched only %d/%d in-population rows, dropped\n",
                sp$file, matched, denom))
    next
  }

  kept <- character(0)
  blk <- list()
  for (nm in cand) {
    v <- d[[nm]]
    if (inherits(v, "haven_labelled")) v <- zz(v)
    if (is.character(v) || inherits(v, "Date")) next
    if (!is.numeric(v)) v <- suppressWarnings(as.numeric(v))
    x <- v[pos]
    x[!in_pop] <- NA_real_          # a woman has no KR row and vice versa
    # Completeness means ITEM non-response, so it is measured among the rows
    # that actually matched. Measuring it across all in-population rows
    # conflates it with roster coverage, and because IR reaches only 76 percent
    # of the MNS women, a threshold of 0.80 then exceeded the attainable
    # maximum and silently rejected every column in the block.
    if (mean(is.finite(x[hit])) < MIN_COMPLETE) next
    u <- unique(x[is.finite(x)])
    if (length(u) < 2L || length(u) > MAX_LEVELS) next
    blk[[nm]] <- x; kept <- c(kept, nm)
  }
  # Never silent. The first version rejected every IR and KR column and said
  # nothing, which reads as "the recode was not there" rather than "the
  # threshold was unreachable". A block that yields nothing must say so.
  if (!length(kept)) {
    cat(sprintf("[wsc4] %-9s pop=%-5s matched %4d/%4d (%.0f%%) but 0 of %d candidates cleared completeness %.2f\n",
                sp$file, sp$pop, matched, denom, 100 * match_rate,
                length(cand), MIN_COMPLETE))
    next
  }
  b <- as.data.frame(blk, stringsAsFactors = FALSE)
  names(b) <- paste0("gw_", names(b))
  blocks[[sp$tag]] <- b
  manifest[[sp$tag]] <- data.frame(
    recode = sp$file, tag = sp$tag, column = names(b),
    class = biomarker_column_class(names(b)),
    population = sp$pop, matched_rows = matched, population_rows = denom,
    # Two distinct quantities, reported separately on purpose.
    match_rate = round(match_rate, 4),
    complete_in_matched = round(vapply(b, function(z)
      mean(is.finite(z[hit])), numeric(1)), 4),
    complete_in_pop = round(vapply(b, function(z)
      mean(is.finite(z[in_pop])), numeric(1)), 4),
    n_levels = vapply(b, function(z) length(unique(z[is.finite(z)])), integer(1)),
    stringsAsFactors = FALSE)
  cat(sprintf("[wsc4] %-9s pop=%-5s matched %4d/%4d (%.0f%%), %4d of %4d candidates kept\n",
              sp$file, sp$pop, matched, denom, 100 * match_rate,
              length(kept), length(cand)))
}

if (!length(blocks)) stop("[wsc4] no block built")
block <- suppressMessages(dplyr::bind_cols(blocks))
# A column can appear in more than one recode (v106 is in IR and KR). Keep the
# first and record the drop rather than silently suffixing.
dup <- duplicated(sub("\\.\\.\\.[0-9]+$", "", names(block)))
names(block) <- sub("\\.\\.\\.[0-9]+$", "", names(block))
if (any(dup)) {
  cat(sprintf("[wsc4] %d duplicate columns across recodes, first kept\n", sum(dup)))
  block <- block[, !dup, drop = FALSE]
}

# ---------------------------------------------------------------------------
# The guard, asserted on the built object rather than on a pattern.
# ---------------------------------------------------------------------------
cl <- biomarker_column_class(names(block))
if (any(cl != "questionnaire")) {
  bad <- names(block)[cl != "questionnaire"]
  stop(sprintf("[wsc4] %d measurement columns reached the block: %s",
               length(bad), paste(head(bad, 10), collapse = ", ")))
}
stopifnot(nrow(block) == nrow(mns))
cat(sprintf("\n[wsc4] block: %d rows x %d columns, all questionnaire-class\n",
            nrow(block), ncol(block)))

attr(block, "wsc4_provenance") <- paste0(
  "Malawi DHS 2015-16 standard recodes, joined to the micronutrient survey on ",
  "(gw_cnum, mnumber, m01) -> (hv001, hv002, hvidx). These are DHS ",
  "questionnaire items, NOT the instrument the other three countries' gw_ ",
  "blocks come from. Built by scripts/accuracy_impact/wsc4_malawi_block.R.")
saveRDS(block, OUT)
readr::write_csv(dplyr::bind_rows(manifest), MAN)
cat(sprintf("-> %s\n-> %s\n",
            file.path("data", "IPD", "Malawi", basename(OUT)),
            file.path("results", "tables", basename(MAN))))
