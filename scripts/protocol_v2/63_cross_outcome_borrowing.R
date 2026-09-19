# =============================================================================
# scripts/protocol_v2/63_cross_outcome_borrowing.R   [XO-01, 2026-09-19]
#
# CAN A COUNTRY BORROW A BIOMARKER IT DID NOT MEASURE FROM THE ONES IT DID?
#
# Leave-one-country-out transport (estimand C), with the held-out country's
# OTHER survey biomarkers offered to the index as predictors while its target
# outcome stays hidden. This is the deployment case "the survey measured iron
# but not folate": the outcome-to-outcome relationship is learned on the
# training countries that measured both and applied through the target
# country's own measured biomarkers.
#
# The block is built in memory from targets_v2.csv (the other outcomes'
# district y_level, rank-normalised within country like every predictor) and
# is NEVER written to predictors_admin2_shared.csv: it is same-survey by
# construction, which the leakage policy exists to keep out of the headline
# vocabulary. Each borrowed biomarker is its own one-column domain, so the
# zero-tuning index gives it its own Fisher-z weight.
#
# Arms (each run through domain_index and domain_enet), per cell:
#   base            covariates only  (GUARD 1: reproduces civ_transport_guards.csv)
#   same_nutrient   + the other population of the same nutrient (child <-> women)
#   other_xpop      + other nutrients measured in the OTHER population
#   other_samepop   + other nutrients measured in the SAME population (same
#                     blood sample: shares person and cluster noise, so the
#                     gain here is an upper bound, not a deployment claim)
#   all_other       + every other biomarker in the block
#   block_only_*    the borrowed biomarkers alone, no covariates
#
# Pre-registered (docs/superpowers/specs/2026-09-19-cross-outcome-borrowing-design.md):
# cross-nutrient delta within +-0.05; same-nutrient cross-population +0.05 to
# +0.10; folate / B12 cells gain nothing.
#
#   V2_DOMAIN_SET=cs|cs_top5|all   Rscript scripts/protocol_v2/63_cross_outcome_borrowing.R
# -> results/tables/protocol_v2/cross_outcome_loco_<set>.csv       cell x arm x model
#    results/tables/protocol_v2/cross_outcome_summary_<set>.csv    per outcome-class x arm
#    results/tables/protocol_v2/cross_outcome_weights_<set>.csv    the index weight each
#                                                                  borrowed biomarker got
#    results/figures/protocol_v2/fig_cross_outcome_delta_<set>.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
set.seed(20260919L)

P2   <- "results/tables/protocol_v2"
FDIR <- "results/figures/protocol_v2"; dir.create(FDIR, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")

TG <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")

DOMAIN_SET <- Sys.getenv("V2_DOMAIN_SET", "cs")
DOMAIN_SETS <- list(
  cs      = c("Climate and weather", "Soil characteristics"),
  cs_top5 = c("Climate and weather", "Soil characteristics", "Anaemia and haemoglobin",
              "Agricultural production, land use", "Infection and inflammation burden"),
  all     = NULL)
if (!DOMAIN_SET %in% names(DOMAIN_SETS)) stop("V2_DOMAIN_SET must be one of: ", paste(names(DOMAIN_SETS), collapse = ", "))
cand <- if (is.null(DOMAIN_SETS[[DOMAIN_SET]])) MD$column else MD$column[MD$domain %in% DOMAIN_SETS[[DOMAIN_SET]]]
PREDS <- drop_near_outcome_v2(intersect(cand, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
OUTCOMES  <- c("child_iron", "child_vitA", "women_iron", "women_vitA", "women_folate", "women_b12")
cat(sprintf("domain set %s: %d covariate columns\n", DOMAIN_SET, length(PREDS)))

# ── the biomarker block: every outcome's district y_level, one column each ──
BLK <- TG |> filter(is.finite(y_level)) |>
  transmute(country, Admin1, Admin2, col = paste0("svy_", outcome), y_level) |>
  pivot_wider(names_from = col, values_from = y_level)
BLOCK_COLS <- setdiff(names(BLK), c("country", "Admin1", "Admin2"))
# each borrowed biomarker is its own domain; the label's first 12 characters
# must stay unique (build_domain_pcs_v2 keys PCs on that prefix)
svy_domain <- function(col) paste0("svy ", sub("^svy_", "", col))
for (b in BLOCK_COLS) domain_of[b] <- svy_domain(b)
stopifnot(!anyDuplicated(substr(svy_domain(BLOCK_COLS), 1, 12)))
# the PC-axis name build_domain_pcs_v2 gives each one-column biomarker domain, mapped back to the outcome
PC_OF <- stats::setNames(sub("^svy_", "", BLOCK_COLS), paste0(make.names(substr(svy_domain(BLOCK_COLS), 1, 12)), "_PC1"))
# GUARD 2: the block never coincides with a covariate column
stopifnot(!any(BLOCK_COLS %in% PREDS))

nutrient_of   <- function(o) sub("^(child|women)_", "", o)
population_of <- function(o) sub("_.*$", "", o)

# ── cells ────────────────────────────────────────────────────────────────────
build_cell <- function(cn, on, target) {
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff"  else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t[, c("Admin1", "Admin2", ycol, wcol)],
                  S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  m <- left_join(m, BLK[BLK$country == cn, c("Admin1", "Admin2", BLOCK_COLS)], by = c("Admin1", "Admin2"))
  blk <- BLOCK_COLS[colSums(is.finite(as.matrix(m[, BLOCK_COLS, drop = FALSE]))) >= 12]   # GUARD 3
  blk <- setdiff(blk, paste0("svy_", on))                                                   # never the target itself
  Xr <- prep_predictors_v2(as.matrix(m[, c(PREDS, blk), drop = FALSE]))
  if (sum(colnames(Xr) %in% PREDS) < 10) return(NULL)
  list(country = cn, n = nrow(m), y_nat = m[[ycol]],
       y_mod = if (target == "prev") .v2_logit(m[[ycol]]) else m[[ycol]],
       X = Xr, w = m[[wcol]], block = intersect(blk, colnames(Xr)))
}

# ── the arms: which columns each one may see ────────────────────────────────
arm_columns <- function(on, block, covariates) {
  nut <- nutrient_of(on); pop <- population_of(on)
  b_out <- sub("^svy_", "", block)
  same_nut <- block[nutrient_of(b_out) == nut & population_of(b_out) != pop]
  oth_x    <- block[nutrient_of(b_out) != nut & population_of(b_out) != pop]
  oth_s    <- block[nutrient_of(b_out) != nut & population_of(b_out) == pop]
  list(
    base                    = list(cols = covariates,                 shares_pop = FALSE),
    same_nutrient           = list(cols = c(covariates, same_nut),    shares_pop = FALSE),
    other_xpop              = list(cols = c(covariates, oth_x),       shares_pop = FALSE),
    other_samepop           = list(cols = c(covariates, oth_s),       shares_pop = TRUE),
    all_other               = list(cols = c(covariates, block),       shares_pop = NA),
    block_only_same_nutrient = list(cols = same_nut,                  shares_pop = FALSE),
    block_only_all          = list(cols = block,                      shares_pop = NA))
}

MODELS <- c(domain_index = "arm_domain_index_v2", domain_enet = "arm_domain_enet_v2")

rows <- list(); wrows <- list()
for (target in c("level", "prev")) for (on in OUTCOMES) {
  cl <- list()
  for (cn in COUNTRIES) { z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL); if (!is.null(z)) cl[[cn]] <- z }
  if (length(cl) < 3) { cat(sprintf("  %-5s %-13s skipped: %d countries\n", target, on, length(cl))); next }
  common_cov <- Reduce(intersect, lapply(cl, function(z) intersect(colnames(z$X), PREDS)))
  common_blk <- Reduce(intersect, lapply(cl, function(z) z$block))   # identical block in every fold
  if (length(common_cov) < 10) next
  Xm   <- do.call(rbind, lapply(cl, function(z) z$X[, c(common_cov, common_blk), drop = FALSE]))
  ctry <- rep(names(cl), vapply(cl, function(z) z$n, 0L))
  Y    <- unlist(lapply(cl, function(z) as.numeric(scale(z$y_mod))))
  ynat <- unlist(lapply(cl, function(z) z$y_nat))
  wv   <- unlist(lapply(cl, function(z) z$w))
  arms <- arm_columns(on, common_blk, common_cov)
  for (a in names(arms)) {
    cols <- arms[[a]]$cols
    if (!length(cols)) next
    if (grepl("^block_only", a) && identical(cols, common_cov)) next
    if (!grepl("^block_only", a) && a != "base" && !length(setdiff(cols, common_cov))) next   # nothing to add
    Xa <- Xm[, cols, drop = FALSE]
    for (mod in names(MODELS)) {
      fn <- get(MODELS[[mod]])
      for (cn in names(cl)) {
        te <- which(ctry == cn); tr <- which(ctry != cn)
        Dm <- domain_representation_v2(Xa, domain_of, sign_rows = tr)
        p  <- tryCatch(fn(tr, te, Y, Xa, Dm, NULL), error = function(e) rep(NA_real_, length(te)))
        s  <- score_v2(ynat[te], p, wv[te], scale = if (target == "prev") "prev" else "level")
        rows[[length(rows) + 1L]] <- data.frame(
          domain_set = DOMAIN_SET, target = target, outcome = on, heldout = cn, arm = a, model = mod,
          n_units = length(te), n_cov = length(intersect(cols, common_cov)), n_block = length(setdiff(cols, common_cov)),
          block_cols = paste(sub("^svy_", "", setdiff(cols, common_cov)), collapse = ";"),
          shares_population = arms[[a]]$shares_pop,
          spearman = s$spearman, pearson = s$pearson, topk = s$topk, stringsAsFactors = FALSE)
        # the weight the index gave each borrowed biomarker, on the training rows of this fold
        if (mod == "domain_index" && length(setdiff(cols, common_cov))) {
          z <- .index_weights_v2(Dm[tr, , drop = FALSE], Y[tr]); sdD <- apply(Dm[tr, , drop = FALSE], 2, sd)
          bcols <- intersect(colnames(Dm), names(PC_OF))
          for (bc in bcols) wrows[[length(wrows) + 1L]] <- data.frame(
            domain_set = DOMAIN_SET, target = target, outcome = on, heldout = cn, arm = a,
            biomarker = PC_OF[[bc]], z = z[[bc]], sd_axis = sdD[[bc]],
            share_abs_contrib = abs(z[[bc]] * sdD[[bc]]) / sum(abs(z * sdD)), stringsAsFactors = FALSE)
        }
      }
    }
  }
  cat(sprintf("  %-5s %-13s %d countries, %d covariates, block: %s\n", target, on, length(cl), length(common_cov),
              paste(sub("^svy_", "", common_blk), collapse = " ")))
}
R <- bind_rows(rows)
R <- R |> group_by(domain_set, target, outcome, heldout, model) |>
  mutate(delta_vs_base = spearman - spearman[arm == "base"][1],
         topk_delta_vs_base = topk - topk[arm == "base"][1]) |> ungroup()
write.csv(R, file.path(P2, sprintf("cross_outcome_loco_%s.csv", DOMAIN_SET)), row.names = FALSE)
W <- bind_rows(wrows); write.csv(W, file.path(P2, sprintf("cross_outcome_weights_%s.csv", DOMAIN_SET)), row.names = FALSE)

# ── GUARD 1: base on cs reproduces the published four-country transport ──────
if (DOMAIN_SET == "cs") {
  G <- read.csv("results/tables/policy_deck/civ_transport_guards.csv", stringsAsFactors = FALSE) |>
    filter(arm == "4 countries only") |> select(outcome, heldout, ref = spearman)
  chk <- R |> filter(target == "level", arm == "base", model == "domain_index") |>
    inner_join(G, by = c("outcome", "heldout")) |> mutate(diff = abs(spearman - ref))
  cat(sprintf("\nGUARD 1: base vs civ_transport_guards.csv on %d cells, max |diff| = %.5f\n", nrow(chk), max(chk$diff)))
  if (max(chk$diff) > 5e-4) { print(chk); stop("GUARD 1 failed: base arm does not reproduce the published transport") }
}

# ── summary ──────────────────────────────────────────────────────────────────
R <- R |> mutate(outcome_class = ifelse(nutrient_of(outcome) %in% c("iron", "vitA"), "iron / vitA", "folate / B12"))
SUMM <- R |> filter(arm != "base") |> group_by(domain_set, target, model, outcome_class, arm, shares_population) |>
  summarise(cells = n(), base_mean = mean(spearman - delta_vs_base, na.rm = TRUE),
            arm_mean = mean(spearman, na.rm = TRUE), mean_delta = mean(delta_vs_base, na.rm = TRUE),
            median_delta = median(delta_vs_base, na.rm = TRUE), cells_up = sum(delta_vs_base > 0, na.rm = TRUE),
            mean_topk_delta = mean(topk_delta_vs_base, na.rm = TRUE), .groups = "drop") |>
  arrange(target, model, outcome_class, arm)
write.csv(SUMM, file.path(P2, sprintf("cross_outcome_summary_%s.csv", DOMAIN_SET)), row.names = FALSE)

cat(sprintf("\n=== XO-01  domain set %s  (Spearman, leave-one-country-out; delta = arm - base) ===\n", DOMAIN_SET))
for (tg in c("level", "prev")) for (mod in names(MODELS)) {
  cat(sprintf("\n--- target %s, model %s ---\n", tg, mod))
  print(as.data.frame(SUMM |> filter(target == tg, model == mod) |>
    select(outcome_class, arm, cells, base_mean, arm_mean, mean_delta, median_delta, cells_up)), digits = 3, row.names = FALSE)
}
cat("\nper-cell deltas, level target, domain_index:\n")
print(as.data.frame(R |> filter(target == "level", model == "domain_index", !grepl("^block_only", arm), arm != "base") |>
  select(outcome, heldout, arm, spearman, delta_vs_base) |>
  pivot_wider(names_from = arm, values_from = c(spearman, delta_vs_base))), digits = 2, row.names = FALSE)
if (nrow(W)) {
  cat("\nindex weight of each borrowed biomarker (all_other arm, level, mean over folds):\n")
  print(as.data.frame(W |> filter(target == "level", arm == "all_other") |> group_by(outcome, biomarker) |>
    summarise(z = mean(z), share = mean(share_abs_contrib), .groups = "drop") |>
    pivot_wider(names_from = biomarker, values_from = c(z, share))), digits = 2, row.names = FALSE)
}

# ── figure: delta vs base per cell, level target ─────────────────────────────
FG <- R |> filter(target == "level", !grepl("^block_only", arm), arm != "base") |>
  mutate(arm = factor(arm, levels = c("same_nutrient", "other_xpop", "other_samepop", "all_other"),
                      labels = c("+ same nutrient,\nother population", "+ other nutrients,\nother population",
                                 "+ other nutrients,\nsame population\n(shares blood sample)", "+ all other\nbiomarkers")),
         model = ifelse(model == "domain_index", "Zero-tuning index", "Domain elastic net"))
p <- ggplot(FG, aes(x = outcome, y = delta_vs_base, colour = heldout)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_hline(yintercept = c(-0.05, 0.05), linetype = "dotted", colour = "grey70") +
  geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
  facet_grid(model ~ arm) + coord_flip() +
  labs(title = sprintf("XO-01: borrowing other survey biomarkers under leave-one-country-out transport (%s covariates)", DOMAIN_SET),
       subtitle = "Change in Spearman against the covariates-only index, per held-out country; dotted lines at the pre-registered +/- 0.05",
       x = NULL, y = "Spearman, arm minus base", colour = "Held-out country") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
ggsave(file.path(FDIR, sprintf("fig_cross_outcome_delta_%s.png", DOMAIN_SET)), p, width = 13, height = 6.5, dpi = 150, bg = "white")
cat("\nDONE\n")
