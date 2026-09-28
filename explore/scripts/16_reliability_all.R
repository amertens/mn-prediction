# =============================================================================
# explore/scripts/16_reliability_all.R   [probe RL-01]
#
# HOW MUCH DISTRICT SIGNAL DOES EACH BIOMARKER ACTUALLY CARRY?
#
# MA-01 found, as a by-product, that district split-half reliability differs
# enormously by biomarker: Malawi vitamin A (RBP) 0.350 against B12 0.696 and
# zinc 0.703. That is a MEASUREMENT-side explanation for a pattern the project
# has been trying to explain with covariates - vitamin A cells are the weak ones
# and B12 the strongest cell on the record. A predictor cannot recover district
# signal the biomarker does not reliably carry.
#
# This extends it to every country x outcome the pipeline defines, using the
# project's own config, loader and population masks so the populations match
# the targets exactly.
#
# TWO RELIABILITIES, AND THE DIFFERENCE BETWEEN THEM IS THE POINT
#   respondent split  split respondents at random within district. This is what
#                     MA-01 reported. Where a district is ONE cluster it counts
#                     the cluster effect as district signal, which CE-01 showed
#                     is a real inflation (16-27% of the within ceiling on
#                     multi-cluster units) and which matters a lot here: the
#                     project records 57% of Gambia, 83% of Ghana and 85% of
#                     Malawi districts as single-cluster.
#   cluster split     split whole CLUSTERS within district, for districts with
#                     at least two. No cluster effect can be counted as district
#                     signal. This is the honest number, on the subset of
#                     districts that can support it.
#
# Both are reported at half sample and Spearman-Brown corrected to full sample,
# r_full = 2 r_half / (1 + r_half), so they are comparable to a reliability
# coefficient rather than to a half-sample correlation.
#
#   Rscript explore/scripts/16_reliability_all.R
#   EXP_SPLITS=400 to tighten the Monte-Carlo error
# -> explore/out/16_reliability_all.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
EXP_ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"
setwd(EXP_ROOT)
source("explore/R/harness.R")
source("R/config.R")
source("R/data_prep.R")

NSPLIT <- as.integer(Sys.getenv("EXP_SPLITS", "300"))
set.seed(20260928L)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))
sb <- function(r) if (!is.finite(r) || r <= -1) NA_real_ else 2 * r / (1 + r)

#' Split-half reliability of district means.
#' @param unit vector defining the units that are split (respondent id, or
#'   cluster id for the cluster split). Whole units move together.
split_half <- function(y, district, unit, nsplit = NSPLIT, min_units = 2L,
                       min_districts = 8L, method = "spearman", wt = NULL) {
  d <- data.frame(y = y, district = as.character(district), unit = as.character(unit),
                  wt = if (is.null(wt)) 1 else wt)
  d$wt[!is.finite(d$wt) | d$wt <= 0] <- NA
  d$wt[is.na(d$wt)] <- stats::median(d$wt, na.rm = TRUE)
  if (!all(is.finite(d$wt))) d$wt <- 1
  d <- d[is.finite(d$y) & !is.na(d$district) & !is.na(d$unit) & is.finite(d$wt), ]
  if (!nrow(d)) return(c(r = NA_real_, districts = NA_real_))
  # keep districts that have enough distinct units to be split at all
  nu <- tapply(d$unit, d$district, function(z) length(unique(z)))
  keep <- names(nu)[nu >= min_units]
  d <- d[d$district %in% keep, ]
  if (dplyr::n_distinct(d$district) < min_districts) return(c(r = NA_real_, districts = dplyr::n_distinct(d$district)))

  key <- paste(d$district, d$unit)
  ud <- unique(data.frame(district = d$district, unit = d$unit, key = key,
                          stringsAsFactors = FALSE))
  out <- rep(NA_real_, nsplit)
  for (b in seq_len(nsplit)) {
    half <- stats::ave(seq_len(nrow(ud)), ud$district,
                       FUN = function(i) sample(rep(c(1L, 2L), length.out = length(i))))
    hof <- stats::setNames(half, ud$key)
    h <- hof[key]
    wm <- function(k) tapply(seq_len(sum(k)), d$district[k], function(i)
      stats::weighted.mean(d$y[k][i], d$wt[k][i]))
    a1 <- wm(h == 1L); a2 <- wm(h == 2L)
    k <- intersect(names(a1), names(a2))
    if (length(k) < min_districts) next
    out[b] <- suppressWarnings(stats::cor(a1[k], a2[k], method = method))
  }
  c(r = mean(out, na.rm = TRUE), districts = dplyr::n_distinct(d$district))
}

CFG <- get_country_configs()
rows <- list()

for (cn in names(CFG)) {
  cc <- CFG[[cn]]
  dat <- tryCatch(load_merged_data(cc$data_path), error = function(e) {
    message("  !! load failed for ", cc$country, ": ", conditionMessage(e)); NULL })
  if (is.null(dat)) next
  a2 <- cc$admin2_col; psu <- cc$psu_col
  if (!all(c(a2, psu) %in% names(dat))) {
    message("  !! ", cc$country, ": missing ", a2, " or ", psu); next }

  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    ycol <- oc$continuous
    if (is.null(ycol) || !ycol %in% names(dat)) {
      message("  -- ", cc$country, " ", on, ": no continuous column (", ycol, ")"); next }
    keep <- outcome_population_mask(dat, cc, oc, label = "[RL-01]")
    d <- dat[keep, , drop = FALSE]
    y <- num(d[[ycol]])
    # match the pipeline's scale: already-log columns are left alone, strictly
    # positive concentrations are logged, anything else is used as supplied
    is_log <- grepl("log", ycol, ignore.case = TRUE)
    if (!is_log && all(y[is.finite(y)] > 0)) y <- log(y)
    ok <- is.finite(y) & !is.na(d[[a2]]) & !is.na(d[[psu]])
    if (sum(ok) < 100) {
      message("  -- ", cc$country, " ", on, ": only ", sum(ok), " usable rows"); next }
    wcol <- cc$weight_col
    wv <- if (!is.null(wcol) && wcol %in% names(d)) num(d[[wcol]])[ok] else rep(1, sum(ok))
    y <- y[ok]; dis <- as.character(d[[a2]][ok]); cl <- as.character(d[[psu]][ok])

    ncl <- tapply(cl, dis, function(z) length(unique(z)))
    nre <- tapply(y, dis, length)
    resp <- split_half(y, dis, seq_along(y), min_units = 4L)
    clus <- split_half(y, dis, cl, min_units = 2L)
    respw <- split_half(y, dis, seq_along(y), min_units = 4L, wt = wv)

    # THE SAME CELL ON THE PREVALENCE SCALE. WS1a (R/reliability_empirical.R)
    # computes reliability of district PREVALENCE from the binary outcome with
    # Pearson; the block above is the continuous LEVEL with Spearman. They are
    # different quantities, and the gap between them prices what dichotomising
    # at a cutoff costs in district signal - which matters because the level
    # target beats the prevalence target throughout the project.
    bcol <- oc$binary
    bin <- if (!is.null(bcol) && bcol %in% names(d)) num(d[[bcol]])[ok] else {
      cut <- oc$cutoff
      if (is.null(cut)) rep(NA_real_, length(y)) else {
        yy <- if (identical(oc$cutoff_scale, "log") || is_log) y else exp(y)
        if (identical(oc$cutoff_dir, "less")) as.numeric(yy < cut) else as.numeric(yy > cut)
      }
    }
    bin[!is.finite(bin)] <- NA
    resp_b <- if (sum(is.finite(bin)) > 100 && length(unique(stats::na.omit(bin))) > 1)
      split_half(bin, dis, seq_along(bin), min_units = 4L, method = "pearson") else
      c(r = NA_real_, districts = NA_real_)

    rows[[paste(cc$country, on)]] <- data.frame(
      country = cc$country, outcome = on, column = ycol,
      n = length(y), districts = dplyr::n_distinct(dis),
      med_resp_per_district = stats::median(nre),
      med_clusters_per_district = stats::median(ncl),
      pct_single_cluster = round(100 * mean(ncl < 2), 1),
      r_half_respondent = unname(resp["r"]),
      r_full_respondent = sb(unname(resp["r"])),
      districts_resp = unname(resp["districts"]),
      r_half_cluster = unname(clus["r"]),
      r_full_cluster = sb(unname(clus["r"])),
      districts_cluster = unname(clus["districts"]),
      r_half_weighted = unname(respw["r"]),
      r_full_weighted = sb(unname(respw["r"])),
      r_half_prev = unname(resp_b["r"]),
      r_full_prev = sb(unname(resp_b["r"])),
      prev_mean = mean(bin, na.rm = TRUE),
      stringsAsFactors = FALSE)
    message(sprintf("  %-12s %-14s n=%5d  districts=%3d  resp r=%.3f  cluster r=%s",
                    cc$country, on, length(y), dplyr::n_distinct(dis),
                    unname(resp["r"]),
                    ifelse(is.finite(clus["r"]), sprintf("%.3f", clus["r"]), "  na")))
  }
}

R <- dplyr::bind_rows(rows)
R$inflation <- round(R$r_full_respondent - R$r_full_cluster, 3)
exp_write(R, "16_reliability_all")

cat("\n== district reliability by cell (Spearman-Brown corrected to full sample) ==\n")
p <- R[order(-R$r_full_respondent),
       c("country", "outcome", "n", "districts", "med_resp_per_district",
         "pct_single_cluster", "r_full_respondent", "r_full_cluster", "inflation")]
print(p, row.names = FALSE, digits = 3)

cat("\n== by NUTRIENT, pooled over countries ==\n")
R$nutrient <- sub("^(child|women)_", "", R$outcome)
agg <- aggregate(cbind(r_full_respondent, r_full_cluster) ~ nutrient, data = R,
                 FUN = function(z) round(mean(z, na.rm = TRUE), 3))
agg$cells <- aggregate(r_full_respondent ~ nutrient, data = R, FUN = length)$r_full_respondent
print(agg[order(-agg$r_full_respondent), ], row.names = FALSE)

cat("\n== by COUNTRY ==\n")
agg2 <- aggregate(cbind(r_full_respondent, r_full_cluster, pct_single_cluster) ~ country,
                  data = R, FUN = function(z) round(mean(z, na.rm = TRUE), 3))
print(agg2[order(-agg2$r_full_respondent), ], row.names = FALSE)

cat("\n== how much of the respondent-split reliability is CLUSTER effect? ==\n")
hv <- R[is.finite(R$r_full_cluster), ]
cat(sprintf("cells with a cluster-split estimate: %d of %d\n", nrow(hv), nrow(R)))
if (nrow(hv)) cat(sprintf("mean inflation (respondent minus cluster split): %+.3f; inflated in %d of %d cells\n",
    mean(hv$inflation, na.rm = TRUE), sum(hv$inflation > 0, na.rm = TRUE), nrow(hv)))

cat("\n== does reliability explain which cells the model predicts well? ==\n")
BM <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/benchmarks_v2_cells.csv"),
               stringsAsFactors = FALSE)
b <- BM[BM$estimand == "infill" & BM$target == "level" & BM$arm == "domain_index",
        c("country", "outcome", "spearman")]
names(b)[3] <- "index_infill_spearman"
m <- merge(R, b, by = c("country", "outcome"))
if (nrow(m) > 4) {
  ct <- suppressWarnings(stats::cor.test(m$r_full_respondent, m$index_infill_spearman,
                                         method = "spearman"))
  cat(sprintf("spearman(district reliability, index in-fill accuracy) = %.3f  p = %.4f  over %d cells\n",
              ct$estimate, ct$p.value, nrow(m)))
  print(m[order(-m$r_full_respondent),
          c("country", "outcome", "r_full_respondent", "index_infill_spearman")],
        row.names = FALSE, digits = 3)
}

cat("
== WHAT DOES DICHOTOMISING COST? continuous level vs binary prevalence ==
")
cat("   (both respondent-split, Spearman-Brown corrected; level uses Spearman,
")
cat("    prevalence uses Pearson on the binary, matching WS1a's definition)

")
D <- R[is.finite(R$r_full_prev), c("country","outcome","prev_mean",
                                   "r_full_respondent","r_full_prev")]
D$cost <- round(D$r_full_respondent - D$r_full_prev, 3)
print(D[order(-D$cost), ], row.names = FALSE, digits = 3)
cat(sprintf("
mean cost of dichotomising: %+.3f  (level more reliable in %d of %d cells)
",
    mean(D$cost, na.rm = TRUE), sum(D$cost > 0, na.rm = TRUE), sum(is.finite(D$cost))))

cat("
== does the reliability gap track the ACCURACY gap between the two targets? ==
")
bl <- BM[BM$estimand == "infill" & BM$arm == "domain_index",
         c("country","outcome","target","spearman")]
bw <- reshape(bl, idvar = c("country","outcome"), timevar = "target", direction = "wide")
names(bw) <- sub("^spearman[.]", "acc_", names(bw))
mm <- merge(D, bw, by = c("country","outcome"))
if (nrow(mm) > 4) {
  mm$acc_gap <- mm$acc_level - mm$acc_prev
  ct <- suppressWarnings(cor.test(mm$cost, mm$acc_gap, method = "spearman"))
  cat(sprintf("spearman(reliability gap, accuracy gap) = %.3f  p = %.4f  over %d cells
",
      ct$estimate, ct$p.value, nrow(mm)))
  cat(sprintf("mean accuracy gap (level minus prevalence) = %+.3f
",
      mean(mm$acc_gap, na.rm = TRUE)))
}

cat("
== ROBUSTNESS: unweighted vs SURVEY-WEIGHTED district means ==
")
cat("   (the pipeline's targets are survey-weighted; the table above is not)

")
Rw <- R[is.finite(R$r_full_weighted), c("country","outcome","r_full_respondent","r_full_weighted")]
Rw$diff <- round(Rw$r_full_weighted - Rw$r_full_respondent, 3)
print(Rw[order(-abs(Rw$diff)), ], row.names = FALSE, digits = 3)
cat(sprintf("
mean |difference| %.4f ; spearman between the two orderings %.3f over %d cells
",
    mean(abs(Rw$diff), na.rm = TRUE),
    suppressWarnings(cor(Rw$r_full_respondent, Rw$r_full_weighted, method = "spearman")),
    nrow(Rw)))
BM2 <- read.csv(file.path(EXP_ROOT, "results/tables/protocol_v2/benchmarks_v2_cells.csv"),
                stringsAsFactors = FALSE)
b2 <- BM2[BM2$estimand == "infill" & BM2$target == "level" & BM2$arm == "domain_index",
          c("country","outcome","spearman")]
names(b2)[3] <- "acc"
m2 <- merge(Rw, b2, by = c("country","outcome"))
if (nrow(m2) > 4) {
  c1 <- suppressWarnings(cor.test(m2$r_full_respondent, m2$acc, method = "spearman"))
  c2 <- suppressWarnings(cor.test(m2$r_full_weighted,   m2$acc, method = "spearman"))
  cat(sprintf("reliability-predicts-accuracy: unweighted %.3f (p=%.4f) | weighted %.3f (p=%.4f), n=%d
",
      c1$estimate, c1$p.value, c2$estimate, c2$p.value, nrow(m2)))
}
