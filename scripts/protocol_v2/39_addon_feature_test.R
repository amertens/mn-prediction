# =============================================================================
# scripts/protocol_v2/39_addon_feature_test.R   [AD-xx]
#
# DOES AN ADD-ON FEATURE BLOCK EARN ITS PLACE? A reusable test.
#
# Any new block of district columns (time-matched food prices, time-matched
# climate, satellite embeddings, ...) is scored the same way, so the answers
# are comparable and nothing enters the vocabulary on a hunch:
#   scan     within-country Spearman of each new column with the biomarker
#            level, per outcome; how many countries agree on the sign
#   in-fill  5-fold district CV (replicated), zero-tuning domain index on
#            base vocabulary / base + add-on / add-on alone / (base with a
#            named block replaced by the add-on) -- identical folds, so the
#            comparison is paired
#   LOCO     leave-one-country-out at the district rung and the regional
#            tier, the same sets plus climate+soil with and without the add-on
# Columns absent for a country are rank-normalised as usual (dropped under
# 70% coverage within that country); for the pooled LOCO matrix a country's
# missing add-on columns are imputed to 0 (the rank-normal median) and the
# imputation is reported. The add-on becomes its own domain, so the domain
# PCs and the index treat it as one block.
#
# Environment:
#   ADDON_FILE    csv with country, Admin1, Admin2 and numeric columns
#   ADDON_DOMAIN  domain label (first 12 characters must not collide with an
#                 existing domain's first 12 characters)
#   ADDON_TAG     output suffix
#   ADDON_DROP    optional regex of existing columns to remove in a 'replace' set
#   AD_REPS       in-fill replicates (default 10)
#
#   ADDON_FILE=... ADDON_DOMAIN=... ADDON_TAG=... Rscript scripts/protocol_v2/39_addon_feature_test.R
# -> results/tables/protocol_v2/addon_<tag>_{scan,infill,loco}.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"
TAG <- Sys.getenv("ADDON_TAG", "addon"); AFILE <- Sys.getenv("ADDON_FILE"); ADOM <- Sys.getenv("ADDON_DOMAIN", "Add-on block")
ADROP <- Sys.getenv("ADDON_DROP", ""); REPS <- as.integer(Sys.getenv("AD_REPS", "10")); MIN_TRAIN <- 20L; set.seed(20260904L)
stopifnot(nzchar(AFILE), file.exists(AFILE))
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
BND <- readRDS("dashboard/data/admin2_boundaries.rds")
COUNTRIES <- c(gambia = "Gambia", ghana = "Ghana", malawi = "Malawi", sierraleone = "SierraLeone")
CENT <- do.call(rbind, lapply(names(COUNTRIES), function(lc) {
  b <- BND[[lc]]; xy <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(b))))
  data.frame(country = COUNTRIES[[lc]], Admin1 = as.character(sf::st_drop_geometry(b)$Admin1),
             Admin2 = as.character(sf::st_drop_geometry(b)$Admin2), lon = xy[, 1], lat = xy[, 2], stringsAsFactors = FALSE) }))
A <- read.csv(AFILE, check.names = FALSE, stringsAsFactors = FALSE); stopifnot(all(c("country", "Admin1", "Admin2") %in% names(A)))
ACOLS <- setdiff(names(A), c("country", "Admin1", "Admin2")); ACOLS <- ACOLS[vapply(ACOLS, function(cc) is.numeric(A[[cc]]), TRUE)]
ACOLS <- setdiff(ACOLS, names(S))
A <- A[!duplicated(A[, c("country", "Admin1", "Admin2")]), c("country", "Admin1", "Admin2", ACOLS)]
S2 <- left_join(S, A, by = c("country", "Admin1", "Admin2"))
jr <- S2 |> group_by(country) |> summarise(n = dplyr::n(), joined = sum(is.finite(.data[[ACOLS[1]]])), .groups = "drop")
cat(sprintf("add-on '%s': %d columns | join rate by country: %s\n", TAG, length(ACOLS), paste(sprintf("%s %d/%d", jr$country, jr$joined, jr$n), collapse = ", ")))
domain_of <- c(domain_of, stats::setNames(rep(ADOM, length(ACOLS)), ACOLS))
domains <- sort(unique(stats::na.omit(domain_of[c(PREDS, ACOLS)])))
pref <- make.names(substr(domains, 1, 12)); if (any(duplicated(pref))) stop("domain prefix collision: ", paste(domains[duplicated(pref) | duplicated(pref, fromLast = TRUE)], collapse = " | "))
prefix_of <- stats::setNames(domains, pref); col_domain <- function(cols) unname(prefix_of[sub("_PC[0-9]+$", "", cols)])
CS <- c("Climate and weather", "Soil characteristics")
PRED_DROP <- if (nzchar(ADROP)) grep(ADROP, PREDS, value = TRUE) else character(0)
SETS <- list(base = PREDS, plus = c(PREDS, ACOLS), addon_only = ACOLS)
if (length(PRED_DROP)) SETS$replace <- c(setdiff(PREDS, PRED_DROP), ACOLS)
cat("sets:", paste(sprintf("%s(%d)", names(SETS), lengths(SETS)), collapse = " "), if (length(PRED_DROP)) sprintf("| replace drops %d columns matching '%s'", length(PRED_DROP), ADROP) else "", "\n")
ALLCOLS <- unique(c(PREDS, ACOLS))

build <- function(cn, on, tier = "admin2") {
  t <- TG[TG$country == cn & TG$outcome == on, ]; t <- t[is.finite(t$y_level) & is.finite(t$y_prev) & is.finite(t$n_eff_cont) & t$n_eff_cont > 0, ]
  if (nrow(t) < 12) return(NULL)
  if (tier == "admin1") {
    a <- t |> group_by(Admin1) |> summarise(y = stats::weighted.mean(y_level, n_eff_cont), yp = stats::weighted.mean(y_prev, n_eff_cont), w = sum(n_eff_cont), .groups = "drop")
    x <- S2[S2$country == cn, c("Admin1", ALLCOLS)] |> group_by(Admin1) |> summarise(across(all_of(ALLCOLS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    c1 <- CENT[CENT$country == cn, ] |> group_by(Admin1) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    m <- a |> inner_join(x, by = "Admin1") |> inner_join(c1, by = "Admin1"); region <- m$Admin1
  } else {
    m <- t[, c("Admin1", "Admin2", "y_level", "y_prev", "n_eff_cont")]; names(m)[3:5] <- c("y", "yp", "w")
    m <- m |> inner_join(S2[S2$country == cn, c("Admin1", "Admin2", ALLCOLS)], by = c("Admin1", "Admin2")) |>
      inner_join(CENT[CENT$country == cn, c("Admin1", "Admin2", "lon", "lat")], by = c("Admin1", "Admin2")); region <- m$Admin1
  }
  m <- m[is.finite(m$y) & is.finite(m$lon), ]; if (nrow(m) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, ALLCOLS, drop = FALSE])); if (ncol(Xr) < 20) return(NULL)
  list(country = cn, n = nrow(m), y = m$y, yp = m$yp, w = m$w, X = Xr, region = region, lon = m$lon, lat = m$lat, addon_kept = intersect(ACOLS, colnames(Xr)))
}
run_index <- function(tr, te, Y, D, aux) { p <- tryCatch(ARMS_V2[["domain_index"]](tr, te, Y, NULL, D, aux), error = function(e) NULL); if (is.null(p) || length(p) != length(te)) rep(NA_real_, length(te)) else p }
sub_D <- function(Dm, cols_set) { keep <- which(col_domain(colnames(Dm)) %in% unique(stats::na.omit(domain_of[cols_set]))); Dm[, keep, drop = FALSE] }

# ── scan ─────────────────────────────────────────────────────────────────────
scan <- list()
for (cn in COUNTRIES) for (on in unique(TG$outcome[TG$country == cn])) { t <- TG[TG$country == cn & TG$outcome == on, ]
  m <- inner_join(t[is.finite(t$y_level), c("Admin1", "Admin2", "y_level")], A[A$country == cn, ], by = c("Admin1", "Admin2")); if (nrow(m) < 8) next
  for (cc in ACOLS) { x <- m[[cc]]; ok <- is.finite(x) & is.finite(m$y_level); if (sum(ok) < 8 || stats::sd(x[ok]) == 0) next
    scan[[length(scan) + 1L]] <- data.frame(column = cc, outcome = on, country = cn, n = sum(ok), rho = suppressWarnings(stats::cor(x[ok], m$y_level[ok], method = "spearman")), stringsAsFactors = FALSE) } }
SC <- bind_rows(scan)
if (nrow(SC)) { SC$z <- atanh(pmax(pmin(SC$rho, 0.999), -0.999)) * sqrt(pmax(SC$n - 3, 1))
  SCS <- SC |> group_by(column, outcome) |> summarise(countries = dplyr::n(), mean_rho = mean(rho), same_sign = max(sum(rho > 0), sum(rho < 0)), meta_z = sum(z) / sqrt(dplyr::n()), .groups = "drop") |> arrange(desc(abs(meta_z)))
  write.csv(SC, file.path(OUTDIR, sprintf("addon_%s_scan.csv", TAG)), row.names = FALSE) }

# ── in-fill ──────────────────────────────────────────────────────────────────
rows <- list(); cells <- TG |> distinct(country, outcome) |> filter(country %in% COUNTRIES)
for (i in seq_len(nrow(cells))) { cn <- cells$country[i]; on <- cells$outcome[i]
  cl <- tryCatch(build(cn, on), error = function(e) NULL); if (is.null(cl)) next
  for (target in c("level", "prev")) { Y <- if (target == "level") cl$y else .v2_logit(cl$yp); ynat <- if (target == "level") cl$y else cl$yp
    aux <- list(lon = cl$lon, lat = cl$lat, Admin1 = cl$region, y_nat = Y)
    for (r in seq_len(REPS)) { folds <- make_folds_v2("kfold_district", cl$n, k = 5, rep_id = r)
      for (sname in names(SETS)) { cols <- intersect(SETS[[sname]], colnames(cl$X)); if (length(cols) < 1) next
        pred <- rep(NA_real_, cl$n)
        for (f in unique(folds)) { te <- which(folds == f); tr <- which(folds != f); if (length(tr) < 12) next
          Dm <- domain_representation_v2(cl$X[, cols, drop = FALSE], domain_of, sign_rows = tr); if (!ncol(Dm)) next
          pred[te] <- run_index(tr, te, Y, Dm, aux) }
        s <- score_v2(ynat, pred, cl$w, scale = if (target == "prev") "prev" else "level")
        rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, target = target, rep = r, set = sname, n_cols = length(cols), addon_cols = length(intersect(cols, ACOLS)), spearman = s$spearman, stringsAsFactors = FALSE) } } }
  cat("infill done", cn, on, "| add-on columns kept:", length(cl$addon_kept), "\n") }
IF <- bind_rows(rows); write.csv(IF, file.path(OUTDIR, sprintf("addon_%s_infill.csv", TAG)), row.names = FALSE)

# ── LOCO ─────────────────────────────────────────────────────────────────────
lrows <- list()
for (tier in c("admin2", "admin1")) for (target in c("level", "prev")) for (on in unique(TG$outcome)) {
  cl <- list(); for (cn in COUNTRIES) { z <- tryCatch(build(cn, on, tier), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) next
  base_common <- Reduce(intersect, lapply(cl, function(z) setdiff(colnames(z$X), ACOLS))); if (length(base_common) < 20) next
  addon_union <- unique(unlist(lapply(cl, function(z) intersect(ACOLS, colnames(z$X))))); if (!length(addon_union)) next
  cols_all <- c(base_common, addon_union)
  Xm <- do.call(rbind, lapply(cl, function(z) { M <- matrix(0, nrow = z$n, ncol = length(cols_all), dimnames = list(NULL, cols_all)); have <- intersect(cols_all, colnames(z$X)); M[, have] <- z$X[, have]; M }))
  imputed <- vapply(cl, function(z) length(setdiff(addon_union, colnames(z$X))), 0L)
  yfun <- function(z) if (target == "level") z$y else .v2_logit(z$yp)
  Y <- unlist(lapply(cl, function(z) as.numeric(scale(yfun(z))))); ynat <- unlist(lapply(cl, function(z) if (target == "level") z$y else z$yp)); wv <- unlist(lapply(cl, function(z) z$w))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L)); aux <- list(Admin1 = paste(ctry, unlist(lapply(cl, function(z) z$region))), y_nat = Y)
  LSETS <- list(base = base_common, plus = cols_all, addon_only = addon_union, cs = base_common[domain_of[base_common] %in% CS], cs_plus = c(base_common[domain_of[base_common] %in% CS], addon_union))
  if (length(PRED_DROP)) LSETS$replace <- c(setdiff(base_common, PRED_DROP), addon_union)
  for (h in names(cl)) { te <- which(ctry == h); tr <- which(ctry != h); if (length(tr) < MIN_TRAIN) next
    for (sname in names(LSETS)) { cols <- LSETS[[sname]]; if (length(cols) < 1) next
      Dm <- domain_representation_v2(Xm[, cols, drop = FALSE], domain_of, sign_rows = tr); if (!ncol(Dm)) next
      p <- run_index(tr, te, Y, Dm, aux); s <- score_v2(ynat[te], p, wv[te], scale = if (target == "prev") "prev" else "level")
      lrows[[length(lrows) + 1L]] <- data.frame(tier = tier, target = target, outcome = on, heldout = h, set = sname, n_cols = length(cols), addon_imputed_heldout = unname(imputed[h]), spearman = s$spearman, stringsAsFactors = FALSE) } }
  cat("loco done", tier, target, on, "| add-on columns imputed:", paste(names(imputed), imputed, collapse = " "), "\n")
}
LO <- bind_rows(lrows); write.csv(LO, file.path(OUTDIR, sprintf("addon_%s_loco.csv", TAG)), row.names = FALSE)

# ── report ───────────────────────────────────────────────────────────────────
cat(sprintf("\n===== AD [%s]: add-on block '%s' (%d columns) =====\n", TAG, ADOM, length(ACOLS)))
if (nrow(SC)) { cat("\n-- scan: strongest column x outcome associations with the biomarker level (meta z over countries; + = more deficiency) --\n")
  print(as.data.frame(head(SCS |> mutate(across(c(mean_rho, meta_z), ~ round(.x, 3))), 15)), row.names = FALSE)
  cat(sprintf("columns with |meta z| > 2 in any outcome: %d of %d; column x outcome pairs with all countries agreeing in sign (>= 3 countries): %d of %d\n",
              dplyr::n_distinct(SCS$column[abs(SCS$meta_z) > 2]), length(ACOLS), sum(SCS$same_sign == SCS$countries & SCS$countries >= 3), sum(SCS$countries >= 3))) }
if (nrow(IF)) { CELL <- IF |> group_by(country, outcome, target, set) |> summarise(spearman = median(spearman, na.rm = TRUE), .groups = "drop")
  cat("\n-- in-fill (5-fold district CV, replicated): mean Spearman over cells --\n")
  print(as.data.frame(CELL |> group_by(target, set) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3), median_rho = round(median(spearman, na.rm = TRUE), 3), .groups = "drop")), row.names = FALSE)
  for (tg in unique(CELL$target)) { W <- pivot_wider(CELL[CELL$target == tg, ], names_from = set, values_from = spearman)
    for (sname in setdiff(names(SETS), "base")) if (sname %in% names(W)) { d <- W[[sname]] - W$base
      cat(sprintf("  %-5s %-10s vs base: better in %2d of %2d | median %+.3f | mean %+.3f\n", tg, sname, sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE), mean(d, na.rm = TRUE))) } } }
if (nrow(LO)) { cat("\n-- LOCO transport: mean Spearman over cells --\n")
  print(as.data.frame(LO |> group_by(tier, target, set) |> summarise(cells = dplyr::n(), mean_rho = round(mean(spearman, na.rm = TRUE), 3), positive = sum(spearman > 0, na.rm = TRUE), .groups = "drop")), row.names = FALSE)
  for (tier in unique(LO$tier)) for (tg in unique(LO$target)) { W <- pivot_wider(LO[LO$tier == tier & LO$target == tg, c("outcome", "heldout", "set", "spearman")], names_from = set, values_from = spearman)
    for (pair in list(c("plus", "base"), c("cs_plus", "cs"), c("addon_only", "base"), c("replace", "base"))) if (all(pair %in% names(W))) { d <- W[[pair[1]]] - W[[pair[2]]]
      cat(sprintf("  %-6s %-5s %-10s vs %-5s: better in %2d of %2d | median %+.3f | mean %+.3f\n", tier, tg, pair[1], pair[2], sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE), mean(d, na.rm = TRUE))) } } }
cat("\nDONE\n")
