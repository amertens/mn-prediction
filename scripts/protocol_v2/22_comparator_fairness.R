# =============================================================================
# scripts/protocol_v2/22_comparator_fairness.R   [CF-01]
#
# HOW MUCH OF "REGIONAL AVERAGES ARE NO BETTER THAN CHANCE" IS THE JACKKNIFE?
#
# The in-fill survey baseline gives each district its region's mean computed
# WITHOUT the district's own respondents (jackknife), so the comparison with
# the model is information-symmetric. But regions here hold 3-5 surveyed
# districts. Removing the worst district from a three-district mean leaves the
# average of its two milder neighbours, so the most deficient districts are
# handed the LOWEST anchors: in Malawi child iron the correlation between a
# district's prevalence and its anchor falls from +0.63 (with self) to +0.19
# (jackknifed), and burden capture from 32% to 13% -- below random.
#
# That is not a bug in the jackknife; it is what a survey that missed the
# district would actually know. But the SIZE of the effect at 3-5 units per
# region means the baseline is also carrying a small-group artefact, and the
# NCE sentence "regional averages are no better than chance" rests on it.
#
# THIS SCRIPT scores four regional-mean comparators on identical data, in-fill:
#   full      region mean INCLUDING the district (what a surveyed district's
#             region average really is; leaks the district's own answer)
#   jk        jackknifed (the protocol-v2 baseline)
#   jk_shrunk jackknifed mean shrunk toward the national mean by k/(k+1),
#             k = number of OTHER districts in the region -- so a region with
#             two other districts is trusted 2/3, with five others 5/6. A
#             cheap empirical-Bayes-style repair of the small-group artefact
#             that never sees the district's own answer.
#   split     the region's other districts split at random in half; anchor =
#             mean of one half (averaged over 20 splits). Also never sees the
#             district; noisier than jk by construction.
# against the model's saved cell-level results (domain_index, in-fill), on
# Spearman and burden capture. If jk_shrunk moves the baseline substantially
# toward `full`, the "no better than chance" wording should be scoped or
# softened; if it does not, the jackknife baseline stands as written.
#
#   Rscript scripts/protocol_v2/22_comparator_fairness.R
# -> results/tables/protocol_v2/comparator_fairness.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
OUTDIR <- "results/tables/protocol_v2"; TOPFRAC <- 0.20; NSPLIT <- 20L
set.seed(20260903L)

TG  <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
POP <- readRDS("dashboard/data/admin2_population.rds")
POP$country <- gsub(" ", "", POP$country)  # FIX 2026-09-04: the file spells "Sierra Leone" with a space; without this the join silently dropped the country
M   <- read.csv(file.path(OUTDIR, "nce_targeting_metrics.csv"))
B   <- read.csv(file.path(OUTDIR, "benchmarks_v2_cells.csv"))
pop_for <- function(on) if (grepl("^child", on)) "pop_child" else "pop_women"

capture <- function(y, pop, score) {
  ok <- is.finite(y) & is.finite(pop) & is.finite(score); if (sum(ok) < 5) return(NA_real_)
  y <- y[ok]; pop <- pop[ok]; score <- score[ok]; b <- y * pop
  k <- max(1L, round(TOPFRAC * length(y))); sel <- order(score, decreasing = TRUE)[seq_len(k)]
  sum(b[sel]) / sum(b)
}
sp <- function(y, s) { ok <- is.finite(y) & is.finite(s)
  if (sum(ok) < 5 || stats::sd(s[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(y[ok], s[ok], method = "spearman")) }

rows <- list()
cells <- unique(M[M$estimand == "infill" & M$arm == "domain_index", c("country", "outcome")])
for (i in seq_len(nrow(cells))) {
  cn <- cells$country[i]; on <- cells$outcome[i]
  t <- TG[TG$country == cn & TG$outcome == on, ]
  pp <- POP[POP$country == cn, c("Admin2", pop_for(on))]; names(pp)[2] <- "pop"
  m <- merge(t, pp, by = "Admin2")
  m <- m[is.finite(m$y_prev) & is.finite(m$n_eff) & m$n_eff > 0 & is.finite(m$pop) & m$pop > 0, ]
  n <- nrow(m); if (n < 12) next
  reg <- as.character(m$Admin1); y <- m$y_prev; w <- m$n_eff
  rs <- tapply(y * w, reg, sum); rn <- tapply(w, reg, sum); cnt <- table(reg)
  nat <- stats::weighted.mean(y, w)

  full <- as.numeric((rs / rn)[reg])
  denom <- rn[reg] - w
  jk <- as.numeric(ifelse(denom > 0, (rs[reg] - y * w) / denom, nat))
  k_other <- as.numeric(cnt[reg]) - 1
  shrink <- k_other / (k_other + 1)
  jk_shrunk <- shrink * jk + (1 - shrink) * nat

  split <- rep(0, n)
  for (s in seq_len(NSPLIT)) {
    anc <- rep(nat, n)
    for (r in unique(reg)) {
      idx <- which(reg == r)
      for (j in idx) {
        others <- setdiff(idx, j); if (!length(others)) next
        half <- if (length(others) == 1) others else sample(others, max(1L, floor(length(others) / 2)))
        anc[j] <- stats::weighted.mean(y[half], w[half])
      }
    }
    split <- split + anc / NSPLIT
  }

  model_sp  <- mean(B$spearman[B$country == cn & B$outcome == on & B$estimand == "infill" &
                                 B$target == "prev" & B$arm == "domain_index"], na.rm = TRUE)
  model_cap <- mean(M$capture_top20[M$country == cn & M$outcome == on & M$estimand == "infill" &
                                      M$arm == "domain_index"], na.rm = TRUE)
  for (a in c("full", "jk", "jk_shrunk", "split")) {
    s <- get(a)
    rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, n_units = n,
      units_per_region = round(n / length(cnt), 1), arm = a,
      spearman = sp(y, s), capture = capture(y, m$pop, s), stringsAsFactors = FALSE)
  }
  rows[[length(rows) + 1L]] <- data.frame(country = cn, outcome = on, n_units = n,
    units_per_region = round(n / length(cnt), 1), arm = "model_domain_index",
    spearman = model_sp, capture = model_cap, stringsAsFactors = FALSE)
}
R <- bind_rows(rows)
write.csv(R, file.path(OUTDIR, "comparator_fairness.csv"), row.names = FALSE)

cat("===== CF-01: regional-mean comparators vs the model, in-fill, ", length(unique(paste(R$country, R$outcome))), " cells =====\n", sep = "")
cat("(burden capture: random = 0.20)\n\n")
S <- R |> group_by(arm) |> summarise(cells = dplyr::n(),
        spearman = round(mean(spearman, na.rm = TRUE), 3),
        capture = round(mean(capture, na.rm = TRUE), 3), .groups = "drop") |> arrange(desc(capture))
print(as.data.frame(S), row.names = FALSE)
W <- tidyr::pivot_wider(R[, c("country", "outcome", "arm", "capture")], names_from = arm, values_from = capture)
hh <- function(a, b) { d <- W[[a]] - W[[b]]
  sprintf("capture  %-18s vs %-18s better in %2d of %2d | median %+.3f", a, b,
          sum(d > 0, na.rm = TRUE), sum(is.finite(d)), median(d, na.rm = TRUE)) }
cat("\n", hh("model_domain_index", "jk"), "\n", hh("model_domain_index", "jk_shrunk"), "\n",
    hh("model_domain_index", "split"), "\n", hh("model_domain_index", "full"), "\n",
    hh("jk_shrunk", "jk"), "\n", hh("full", "jk"), "\n")
cat("\nby units-per-region band (capture):\n")
R$band <- cut(R$units_per_region, c(0, 3.5, 5, 100), labels = c("<=3.5", "3.5-5", ">5"))
print(as.data.frame(R |> group_by(band, arm) |>
  summarise(capture = round(mean(capture, na.rm = TRUE), 3), .groups = "drop") |>
  tidyr::pivot_wider(names_from = arm, values_from = capture)), row.names = FALSE)
cat("\nDONE\n")
