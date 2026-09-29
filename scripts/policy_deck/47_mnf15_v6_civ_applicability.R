# =============================================================================
# scripts/policy_deck/47_mnf15_v6_civ_applicability.R
#
# "Is Cote d'Ivoire like the places we learned from?" (MNF15 v6 talk)
#
# Area of applicability (Meyer & Pebesma 2021, https://arxiv.org/abs/2005.07939;
# CAST tutorial cast04-AOA-tutorial) of the climate + soil index that ranks Cote
# d'Ivoire's 33 districts (scripts/policy_deck/04_civ_climate_soil_prediction.R).
# It asks only whether each district's climate and soil sit inside the range of
# the surveyed districts the model learned from. It says nothing about level
# offsets, assays or accuracy.
#
# Columns: exactly the climate + soil columns the Cote d'Ivoire model uses:
#   the "Climate and weather" + "Soil characteristics" columns of the shared set,
#   through drop_near_outcome_v2() at tiers open,survey_public, then
#   prep_predictors_v2() per surveyed-country cell and on Cote d'Ivoire (70%
#   coverage, non-constant), and the intersection over the pool, as loco_pass()
#   of script 04 builds it. Checked identical over the six outcomes.
#
# Method (fixed before looking at results; base R, CAST not installed):
#   scale     pooled RAW values (not the within-country ranks, where every country
#             overlaps by construction); each column standardised by the training
#             rows' mean and SD; any missing value -> training median
#   DI        distance to the nearest training district / mean pairwise distance
#             among training districts (unweighted Euclidean, all columns equal)
#   threshold country-blocked CV on the training set: each training district's DI
#             to the nearest training district in ANOTHER country; threshold =
#             min(Q3 + 1.5 IQR, max), as CAST's .di_threshold() (aoa-helpers.R)
#   inside    DI <= threshold (CAST)
#   rows      Cote d'Ivoire: training = all surveyed districts of the four countries
#             each surveyed country held out: training = the other three countries,
#             threshold from their own country-blocked CV
#
# Outputs (new files only):
#   results/figures/mnf15_v6/a6_civ_applicability.png
#   results/tables/policy_deck/civ_applicability.csv
#     one row per district: country (CoteDIvoire = trained on all four; a surveyed
#     country = that country held out), district, DI, threshold, inside
#
#   Rscript scripts/policy_deck/47_mnf15_v6_civ_applicability.R
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
if (!nzchar(Sys.getenv("V2_PREDICTOR_TIERS"))) Sys.setenv(V2_PREDICTOR_TIERS = "open,survey_public")
source("R/protocol_v2.R")

P2  <- "results/tables/protocol_v2"; PD <- "results/tables/policy_deck"
OUT <- "results/figures/mnf15_v6"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
chk <- function(ok, msg) if (!isTRUE(ok)) stop("CHECK FAILED: ", msg)

TG  <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
CIV <- readRDS("results/transportability/civ_canonical_admin2_full_v2.rds")
chk(is.data.frame(CIV) && nrow(CIV) == 33 && all(c("Admin1", "Admin2") %in% names(CIV)), "CIV database: 33 districts with Admin1/Admin2")

CS_DOMAINS <- c("Climate and weather", "Soil characteristics")
COUNTRIES  <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
OUTCOMES   <- c("child_iron", "child_vitA", "women_iron", "women_vitA", "women_folate", "women_b12")
PREDS <- drop_near_outcome_v2(intersect(MD$column[MD$domain %in% CS_DOMAINS], names(S)), MD)
chk(all(PREDS %in% names(CIV)), "every climate + soil column exists in the CIV database")
cat(sprintf("climate + soil vocabulary after the predictor policy: %d columns\n", length(PREDS)))

# ── 1. the columns the Cote d'Ivoire model uses (script 04, loco_pass with CIV) ──
cell_rows <- function(cn, on) {          # as build_cell() of script 04 (level target)
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t$y_level) & is.finite(t$n_eff_cont), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t[, c("Admin1", "Admin2")], S[S$country == cn, c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 10) return(NULL)
  list(cols = colnames(Xr), keys = paste(m$Admin1, m$Admin2, sep = "|"))
}
civ_cols <- colnames(prep_predictors_v2(as.matrix(CIV[, PREDS])))
common_by_outcome <- list(); common4_by_outcome <- list(); KEYS <- list()
for (on in OUTCOMES) {
  cl <- list()
  for (cn in COUNTRIES) {
    z <- cell_rows(cn, on)
    if (!is.null(z)) { cl[[cn]] <- z$cols; KEYS[[cn]] <- union(KEYS[[cn]], z$keys) }
  }
  chk(length(cl) >= 3, paste("at least three surveyed countries for", on))
  common_by_outcome[[on]]  <- Reduce(intersect, c(cl, list(CIV = civ_cols)))
  common4_by_outcome[[on]] <- Reduce(intersect, cl)
}
chk(length(unique(lapply(common_by_outcome, sort))) == 1, "the CIV model's column set is the same for every outcome")
COLS <- common_by_outcome[[1]]
COLS4 <- common4_by_outcome[[1]]
chk(length(unique(lapply(common4_by_outcome, sort))) == 1, "the four-country column set is the same for every outcome")
cat(sprintf("columns the Cote d'Ivoire model uses: %d (four-country model without CIV: %d)\n", length(COLS), length(COLS4)))
G <- read.csv(file.path(PD, "civ_transport_guards.csv"), stringsAsFactors = FALSE)
chk(all(G$n_common[G$arm == "with CIV in the vocabulary"] == length(COLS)), "column count matches civ_transport_guards.csv (with CIV)")
chk(all(G$n_common[G$arm == "4 countries only"] == length(COLS4)), "column count matches civ_transport_guards.csv (4 countries only)")

# training rows: every surveyed district (union over outcomes), raw values
TR <- bind_rows(lapply(COUNTRIES, function(cn) {
  s <- S[S$country == cn, ]
  s[paste(s$Admin1, s$Admin2, sep = "|") %in% KEYS[[cn]], c("country", "Admin1", "Admin2", PREDS)]
}))
chk(all(table(TR$country)[COUNTRIES] == c(30, 75, 87, 14)), "surveyed districts: Gambia 30, Ghana 75, Malawi 87, Sierra Leone 14")

# why each column left: coverage < 70% or constant, per surveyed country and in CIV
why <- function(x) {
  if (mean(is.finite(x)) < 0.70) sprintf("%d of %d districts have a value", sum(is.finite(x)), length(x))
  else if (length(unique(x[is.finite(x)])) < 2) sprintf("constant (%s) across its districts", format(unique(x[is.finite(x)])))
  else NA_character_
}
dropped <- setdiff(PREDS, COLS)
DR <- bind_rows(lapply(dropped, function(j) {
  w <- c(CIV = why(CIV[[j]]), sapply(COUNTRIES, function(cn) why(TR[[j]][TR$country == cn])))
  data.frame(column = j, where = names(w)[!is.na(w)], reason = w[!is.na(w)], row.names = NULL)
}))
cat(sprintf("\n=== %d columns dropped by the model's own rules ===\n", length(dropped)))
print(DR, row.names = FALSE)

# SoilGrids in CIV: coverage and scale against the surveyed districts
sg <- grep("^soilgrids_", PREDS, value = TRUE)
SG <- data.frame(column = sg,
                 civ_present = sapply(sg, function(j) sum(is.finite(CIV[[j]]))),
                 civ_median = sapply(sg, function(j) median(CIV[[j]], na.rm = TRUE)),
                 train_median = sapply(sg, function(j) median(TR[[j]], na.rm = TRUE)),
                 train_min = sapply(sg, function(j) min(TR[[j]], na.rm = TRUE)),
                 train_max = sapply(sg, function(j) max(TR[[j]], na.rm = TRUE)), row.names = NULL)
SG$ratio <- SG$civ_median / SG$train_median
cat("\n=== CIV SoilGrids columns: coverage and scale ===\n"); print(SG, digits = 4, row.names = FALSE)
chk(all(!sg %in% COLS), "the 70% rule drops every SoilGrids column")

# scale guard on the kept columns: raw values must be on one scale across countries
SC <- data.frame(column = COLS,
                 civ_med = sapply(COLS, function(j) median(CIV[[j]], na.rm = TRUE)),
                 lo = sapply(COLS, function(j) min(TR[[j]], na.rm = TRUE)),
                 hi = sapply(COLS, function(j) max(TR[[j]], na.rm = TRUE)), row.names = NULL)
SC$civ_med_inside <- SC$civ_med >= SC$lo & SC$civ_med <= SC$hi
cat(sprintf("\nscale guard: CIV median inside the surveyed range for %d of %d kept columns\n", sum(SC$civ_med_inside), nrow(SC)))
if (any(!SC$civ_med_inside)) print(SC[!SC$civ_med_inside, ], row.names = FALSE)
cat(sprintf("missing values in the kept columns: training %d, CIV %d\n",
            sum(!is.finite(as.matrix(TR[, COLS]))), sum(!is.finite(as.matrix(CIV[, COLS])))))

# ── 2. dissimilarity index and threshold (CAST, unweighted) ─────────────────────
aoa_base <- function(Xtr, grp, Xnew) {
  Xtr <- as.matrix(Xtr); Xnew <- as.matrix(Xnew)
  med <- apply(Xtr, 2, median, na.rm = TRUE)
  imp <- function(X) { for (j in seq_len(ncol(X))) X[!is.finite(X[, j]), j] <- med[j]; X }
  Xtr <- imp(Xtr); Xnew <- imp(Xnew)
  mu <- colMeans(Xtr); s <- apply(Xtr, 2, sd)
  chk(all(s > 0), "every column varies in the training rows")
  Ztr  <- sweep(sweep(Xtr, 2, mu), 2, s, "/")
  Znew <- sweep(sweep(Xnew, 2, mu), 2, s, "/")
  D <- as.matrix(dist(Ztr)); diag(D) <- NA
  dbar <- mean(D, na.rm = TRUE)                                  # mean pairwise distance
  Dcv <- D; Dcv[outer(grp, grp, "==")] <- NA                     # country-blocked CV
  train_di <- apply(Dcv, 1, min, na.rm = TRUE) / dbar
  thres <- min(unname(quantile(train_di, 0.75) + 1.5 * IQR(train_di)), max(train_di))
  thres_obs <- unname(grDevices::boxplot.stats(unname(train_di))$stats[5])   # sensitivity only: older whisker convention
  tZ <- t(Ztr)
  di <- apply(Znew, 1, function(z) sqrt(min(colSums((tZ - z)^2)))) / dbar
  list(DI = di, threshold = thres, threshold_obs = thres_obs, train_di = train_di, dbar = dbar)
}

run_all <- function(cols) {
  res <- list()
  a <- aoa_base(TR[, cols], TR$country, CIV[, cols])
  res[["CoteDIvoire"]] <- data.frame(country = "CoteDIvoire", Admin1 = CIV$Admin1, Admin2 = CIV$Admin2,
                                     DI = a$DI, threshold = a$threshold, threshold_obs = a$threshold_obs, dbar = a$dbar)
  for (cn in COUNTRIES) {
    tr <- TR$country != cn
    a <- aoa_base(TR[tr, cols], TR$country[tr], TR[!tr, cols])
    res[[cn]] <- data.frame(country = cn, Admin1 = TR$Admin1[!tr], Admin2 = TR$Admin2[!tr],
                            DI = a$DI, threshold = a$threshold, threshold_obs = a$threshold_obs, dbar = a$dbar)
  }
  out <- bind_rows(res); out$inside <- out$DI <= out$threshold
  out
}
A <- run_all(COLS)
SUMM <- A |> group_by(country) |>
  summarise(n = n(), inside = sum(inside), threshold = first(threshold), mean_pairwise = first(dbar),
            median_ratio = median(DI / threshold), max_ratio = max(DI / threshold), .groups = "drop") |>
  mutate(share_inside = inside / n)
cat(sprintf("\n=== area of applicability, %d columns ===\n", length(COLS)))
print(as.data.frame(SUMM), digits = 3, row.names = FALSE)

civ <- A[A$country == "CoteDIvoire", ]
cat("\nCote d'Ivoire districts nearest the edge (DI / threshold):\n")
print(head(civ[order(-civ$DI), c("Admin1", "Admin2", "DI", "threshold")] |> mutate(ratio = DI / threshold), 5), digits = 3, row.names = FALSE)

# sensitivity (printed only) 1: the threshold as the largest training DI below the
# fence (boxplot.stats, the older whisker convention) instead of CAST's min(fence, max)
cat("\nsensitivity: threshold = largest training DI below Q3 + 1.5 IQR (boxplot.stats)\n")
print(as.data.frame(A |> group_by(country) |>
  summarise(n = n(), inside_cast = sum(inside), inside_whisker_obs = sum(DI <= threshold_obs),
            threshold_cast = first(threshold), threshold_whisker_obs = first(threshold_obs), .groups = "drop")),
  digits = 3, row.names = FALSE)

# sensitivity 2: held-out rows on the four-country model's own columns (these
# include the 7 SoilGrids columns, so the CIV row is not meaningful and is omitted)
A4 <- run_all(COLS4)
cat(sprintf("\nsensitivity: held-out rows on the %d columns of the four-country model (no CIV in the vocabulary)\n", length(COLS4)))
print(as.data.frame(A4 |> filter(country != "CoteDIvoire") |> group_by(country) |>
  summarise(n = n(), inside = sum(inside), threshold = first(threshold), .groups = "drop")),
  digits = 3, row.names = FALSE)

# ── 3. held-out accuracy next to the inside share (no calibration possible) ────
PP <- read.csv(file.path(PD, "v3_percell_pairs.csv"), stringsAsFactors = FALSE) |>
  filter(estimand == "country", arm == "domain_index") |>
  group_by(country) |> summarise(outcomes = n(), pairs_all_domain_index = mean(pairs), .groups = "drop")
GS <- G |> filter(arm == "with CIV in the vocabulary") |>
  group_by(country = heldout) |> summarise(spearman_climate_soil = mean(spearman), .groups = "drop")
CAL <- SUMM |> filter(country != "CoteDIvoire") |> select(country, n, inside, share_inside) |>
  left_join(PP, by = "country") |> left_join(GS, by = "country") |> arrange(desc(share_inside))
cat("\n=== held out: share inside vs accuracy when held out ===\n")
cat("pairs = share of district pairs in the survey's order, v3_percell_pairs.csv (country, domain_index = ALL-domain index)\n")
cat("spearman_climate_soil = the map's own model, civ_transport_guards.csv (with CIV in the vocabulary)\n")
print(as.data.frame(CAL), digits = 3, row.names = FALSE)

# ── 4. outputs ──────────────────────────────────────────────────────────────────
dup <- A |> count(country, Admin2) |> filter(n > 1)
OUTCSV <- A |> transmute(country, district = Admin2, DI = signif(DI, 6), threshold = signif(threshold, 6), inside)
if (nrow(dup)) OUTCSV <- A |> transmute(country, admin1 = Admin1, district = Admin2, DI = signif(DI, 6), threshold = signif(threshold, 6), inside)
write.csv(OUTCSV, file.path(PD, "civ_applicability.csv"), row.names = FALSE)
cat("\nwrote", file.path(PD, "civ_applicability.csv"), nrow(OUTCSV), "rows\n")

TEAL <- "#0F7B8A"; ORANGE <- "#B45309"; GUIDE <- "#9AA0A6"; INK <- "#1A1A1A"; SOFT <- "#5F6368"
NAME <- c(CoteDIvoire = "C\u00f4te d'Ivoire", Gambia = "The Gambia", Ghana = "Ghana", Malawi = "Malawi", SierraLeone = "Sierra Leone")
held <- CAL$country                                         # most inside first
ord  <- c("CoteDIvoire", held)
ypos <- setNames(rev(seq_along(ord)), ord)
ypos["CoteDIvoire"] <- ypos["CoteDIvoire"] + 0.35           # a little air under Cote d'Ivoire

PL <- A |> mutate(ratio = DI / threshold, y = ypos[country])
set.seed(47L)
PL <- PL |> group_by(country) |> mutate(yj = y + (runif(n()) - 0.5) * 0.42) |> ungroup()
RT <- SUMM |> mutate(y = ypos[country], txt = sprintf("%d of %d districts inside", inside, n))
xmax <- max(PL$ratio) * 1.07
xlab <- -0.035 * xmax                                       # row labels, left of the panel (clip off)
civ_y <- unname(ypos["CoteDIvoire"])
LB <- data.frame(y = ypos[held], name = NAME[held])

p <- ggplot(PL) +
  annotate("segment", x = 1, xend = 1, y = min(ypos) - 0.35, yend = civ_y + 0.72, colour = GUIDE, linewidth = 0.9) +
  annotate("text", x = 1, y = civ_y + 0.82, label = "edge of the training range", colour = SOFT, size = 5.2, vjust = 0) +
  annotate("text", x = 0.97, y = civ_y + 0.52, label = "inside", colour = TEAL, size = 4.8, hjust = 1, fontface = "bold") +
  annotate("text", x = 1.03, y = civ_y + 0.52, label = "outside", colour = ORANGE, size = 4.8, hjust = 0, fontface = "bold") +
  annotate("segment", x = xlab * 6, xend = xmax, y = civ_y - 0.66, yend = civ_y - 0.66, colour = "grey85", linewidth = 0.4) +
  geom_point(aes(ratio, yj, fill = inside), shape = 21, colour = "white", stroke = 0.35, size = 3.1, alpha = 0.9) +
  annotate("text", x = xlab, y = civ_y, label = NAME[["CoteDIvoire"]], hjust = 1, size = 5.8, fontface = "bold", colour = INK) +
  geom_text(data = LB, aes(x = xlab, y = y + 0.13, label = name), hjust = 1, size = 5.5, colour = INK) +
  geom_text(data = LB, aes(x = xlab, y = y - 0.17), label = "when held out", hjust = 1, size = 4.3, colour = SOFT) +
  geom_text(data = RT, aes(x = Inf, y = y, label = txt), hjust = -0.06, size = 5.2, colour = INK) +
  scale_fill_manual(values = c(`TRUE` = TEAL, `FALSE` = ORANGE), guide = "none") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = scales::breaks_pretty(6)) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(min(ypos) - 0.35, civ_y + 1.12), expand = FALSE, clip = "off") +
  labs(x = "How unlike the training districts (1 = edge of the range the model learned from)", y = NULL) +
  theme_minimal(base_size = 17) +
  theme(text = element_text(colour = INK),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 13.5, colour = "grey30"),
        axis.title.x = element_text(size = 14.5, colour = "grey20", margin = margin(t = 8)),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
        plot.margin = margin(6, 205, 6, 140))
ggsave(file.path(OUT, "a6_civ_applicability.png"), p, width = 11, height = 4.6, dpi = 220, bg = "white")
cat("wrote", file.path(OUT, "a6_civ_applicability.png"), "\n")
