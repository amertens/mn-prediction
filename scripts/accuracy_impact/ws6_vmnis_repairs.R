# =============================================================================
# scripts/accuracy_impact/ws6_vmnis_repairs.R
#
# WS6. Repairs to the national (VMNIS) layer.
#
# WS6a  The sampling term in national_noise_ceiling() implies an effective
#       sample size of order ten. Find out why, correct it, and report the old
#       and new component tables side by side.
# WS6b  Two panels return sd_resid exactly 0.000, which means the lmer fit is at
#       a variance boundary and r_max is not trustworthy there.
# WS6d  Report signed bias and MAE separately in the composition table, and mark
#       the near-degenerate cells.
#
# WHAT THE STATED NUMBER IMPLIES
# ------------------------------
# sd_sampling 0.816 on the logit scale at a prevalence near 25 percent
# (source: results/tables/frozen_2026-09/national_vmnis_ceiling.csv, row 1).
# The delta-method variance of logit(p) is 1 / (n p (1-p)), so
#   0.816^2 = 1.5 / ((n - 1) x 0.25 x 0.75)  =>  n - 1 ~ 12.
# A national nutrition survey does not have twelve respondents. Something in the
# computation is not the sample size.
#
#   Rscript scripts/accuracy_impact/ws6_vmnis_repairs.R
# -> results/tables/vmnis_sampling_audit.csv
# -> results/tables/national_vmnis_ceiling_revised.csv
# -> results/tables/national_composition_revised.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
targets::tar_source(here("R"))
TDIR <- here("results", "tables")

nat <- vmnis_national()
adj_level <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(grepl("inflam", x, ignore.case = TRUE) ~ "inflammation_adjusted",
                   grepl("^none$", x, ignore.case = TRUE) ~ "unadjusted",
                   x == "" | grepl("not specified", x, ignore.case = TRUE) ~ "unspecified",
                   TRUE ~ "other")
}
panels <- unique(do.call(rbind, lapply(NAT_PANEL_MAP, function(x)
  data.frame(mn = x[1], pop = x[2], stringsAsFactors = FALSE))))

# ---------------------------------------------------------------------------
# WS6a. Where the sampling term comes from.
#
# Four candidate mechanisms, each measurable:
#   (1) v_bar is the ARITHMETIC MEAN of a reciprocal. The mean of 1/n is
#       dominated by the smallest n, not by the typical one.
#   (2) prevalence is clamped at 0.005, so p(1-p) can fall to 0.005, and a
#       single near-zero-prevalence survey contributes up to 301/(n-1).
#   (3) rows with a missing Samplesize drop out of v_bar but stay in the lmer
#       fit, so the two are computed on different sets.
#   (4) the ceiling is computed on UN-AGGREGATED VMNIS records while the LOCO
#       model is scored on country-year means, so they describe different
#       populations.
# ---------------------------------------------------------------------------
audit <- list()
for (i in seq_len(nrow(panels))) {
  mn <- panels$mn[i]; pp <- panels$pop[i]; lab <- paste(mn, "|", pp)
  d <- nat |> filter(mn_group == mn, pop == pp)
  d <- d[is.finite(d$prev) & !is.na(d$iso3c), , drop = FALSE]
  if (!nrow(d)) next
  n <- suppressWarnings(as.numeric(d$Samplesize))
  p_cl <- pmin(pmax(d$prev, 0.005), 0.995)
  v_s  <- 1.5 / (pmax(n, 2) - 1) / (p_cl * (1 - p_cl))
  vf   <- v_s[is.finite(v_s)]
  # The same term with the median instead of the mean, and with the clamp
  # loosened to a value a real survey can actually attain.
  p_lo <- pmin(pmax(d$prev, 0.02), 0.98)
  v_alt <- 1.5 / (pmax(n, 2) - 1) / (p_lo * (1 - p_lo))
  audit[[lab]] <- data.frame(
    panel = lab, rows = nrow(d),
    n_missing_samplesize = sum(!is.finite(n)),
    n_min = suppressWarnings(min(n, na.rm = TRUE)),
    n_median = stats::median(n, na.rm = TRUE),
    n_max = suppressWarnings(max(n, na.rm = TRUE)),
    n_under_50 = sum(n < 50, na.rm = TRUE),
    prev_under_2pct = sum(d$prev < 0.02, na.rm = TRUE),
    v_bar_mean = round(mean(vf, na.rm = TRUE), 4),
    v_bar_median = round(stats::median(vf, na.rm = TRUE), 4),
    sd_samp_from_mean = round(sqrt(mean(vf, na.rm = TRUE)), 4),
    sd_samp_from_median = round(sqrt(stats::median(vf, na.rm = TRUE)), 4),
    sd_samp_median_softclamp = round(sqrt(stats::median(v_alt[is.finite(v_alt)],
                                                        na.rm = TRUE)), 4),
    implied_n_from_mean = round(1.5 / mean(vf, na.rm = TRUE) /
                                  (0.25 * 0.75) + 1, 1),
    stringsAsFactors = FALSE)
}
aud <- bind_rows(audit)
readr::write_csv(aud, file.path(TDIR, "vmnis_sampling_audit.csv"))
cat("=== WS6a: the sampling term, decomposed ===\n")
print(as.data.frame(aud[, c("panel","rows","n_missing_samplesize","n_min",
                            "n_median","n_under_50","prev_under_2pct",
                            "sd_samp_from_mean","sd_samp_from_median",
                            "implied_n_from_mean")]), row.names = FALSE)

# ---------------------------------------------------------------------------
# The corrected ceiling.
#
# Three changes, each justified by the audit above:
#   - the sampling term uses the MEDIAN survey's variance, not the mean of a
#     reciprocal, so one tiny survey cannot set the error floor for a panel;
#   - the clamp moves from 0.005 to 0.02, above which the delta method is not
#     dominated by its own boundary;
#   - rows with no Samplesize are excluded from the lmer fit as well, so the
#     variance components and the sampling term describe the same rows.
# A boundary flag is carried so WS6b can mark unusable rows rather than
# quietly reporting them.
# ---------------------------------------------------------------------------
# national_noise_ceiling_v2() has been PROMOTED into R/national_vmnis.R as the
# definition of national_noise_ceiling(). The published behaviour is still
# reachable through its arguments, so the two can be compared in one run:
#   published: centre = "mean",   clamp = 0.005, require_n = FALSE
#   revised:   centre = "median", clamp = 0.02,  require_n = TRUE
national_noise_ceiling_published <- function(d, deff = 1.5)
  national_noise_ceiling(d, deff = deff, clamp = 0.005, centre = "mean",
                         require_n = FALSE)

old <- list(); new <- list()
for (i in seq_len(nrow(panels))) {
  mn <- panels$mn[i]; pp <- panels$pop[i]; lab <- paste(mn, "|", pp)
  d <- nat |> filter(mn_group == mn, pop == pp) |>
    mutate(method = paste(adj_level(Dataadjustedfor), trimws(as.character(Indicator))))
  o <- tryCatch(national_noise_ceiling_published(d), error = function(e) NULL)
  nw <- tryCatch(national_noise_ceiling(d), error = function(e) NULL)
  if (!is.null(o))  old[[lab]] <- cbind(panel = lab, version = "published", o)
  if (!is.null(nw)) new[[lab]] <- cbind(panel = lab, version = "revised", nw)
}
ot <- bind_rows(old); nt <- bind_rows(new)
both <- bind_rows(ot, nt)
readr::write_csv(both, file.path(TDIR, "national_vmnis_ceiling_revised.csv"))

cat("\n=== WS6a/WS6b: published against revised, side by side ===\n")
print(as.data.frame(both[order(both$panel, both$version),
  c("panel","version","surveys","sd_country","sd_method","sd_resid",
    "sd_sampling","r_max_report","r_max_standardised")]), row.names = FALSE)
cat("\n=== WS6b: which rows are usable? ===\n")
print(as.data.frame(nt[, c("panel","sd_resid","resid_at_boundary",
                           "sampling_exceeds_resid","usable")]), row.names = FALSE)

# ---------------------------------------------------------------------------
# WS6d. Signed bias and MAE reported separately, with the near-degenerate cells
# marked. The oracle removes the LEVEL error; it cannot remove a pattern error,
# and a cell whose true prevalence is 1.3 to 2.5 percent has almost no level to
# get wrong.
# ---------------------------------------------------------------------------
comp <- read.csv(file.path(TDIR, "national_composition.csv"), stringsAsFactors = FALSE)
lev  <- tryCatch(read.csv(file.path(TDIR, "national_composition_levels.csv"),
                          stringsAsFactors = FALSE), error = function(e) NULL)
comp$true_pp <- if (!is.null(lev) && "true_national_pp" %in% names(lev))
  lev$true_national_pp[match(paste(comp$outcome, tolower(comp$country)),
                             paste(lev$outcome, tolower(lev$country)))] else NA_real_
# "Near-degenerate" means there is almost no level to get wrong. The four
# women's vitamin A cells sit at 1.3 to 2.5 percent true prevalence, so an arm
# that removes level error can barely improve on one that does nothing, and the
# oracle's apparent weakness there is arithmetic rather than evidence.
comp$near_degenerate <- is.finite(comp$true_pp) & comp$true_pp < 5

s <- comp |> group_by(arm) |>
  summarise(cells = dplyr::n(),
            MAE_pp = round(mean(mae_pp), 2),
            mean_signed_bias_pp = round(mean(bias_pp), 2),
            mean_abs_bias_pp = round(mean(abs(bias_pp)), 2),
            .groups = "drop")
sd_ <- comp |> filter(!near_degenerate | is.na(near_degenerate)) |> group_by(arm) |>
  summarise(cells_excl_degenerate = dplyr::n(),
            MAE_pp_excl = round(mean(mae_pp), 2),
            mean_signed_bias_excl = round(mean(bias_pp), 2), .groups = "drop")
out <- dplyr::left_join(s, sd_, by = "arm")
readr::write_csv(out, file.path(TDIR, "national_composition_revised.csv"))
cat("\n=== WS6d: composition arms, signed bias separated from MAE ===\n")
print(as.data.frame(out), row.names = FALSE)
cat(sprintf("\nnear-degenerate cells (true national prevalence below 5 pp): %d of %d\n",
            sum(comp$near_degenerate, na.rm = TRUE) / length(unique(comp$arm)),
            nrow(comp) / length(unique(comp$arm))))
