# =============================================================================
# scripts/accuracy_impact/ws5_anchoring_budget.R
#
# WS5. How much survey does a regional anchor need?
#
# THE QUESTION, AND WHY WS2 CHANGES ITS SHAPE
# -------------------------------------------
# Open question 3 of docs/SESSION_FINDINGS_FOR_REVIEW.md asks how few and how
# thinly sampled the anchor regions can be before the Section 4 gain disappears.
# WS2 established that the Section 4 gain, as measured, does not survive a
# jackknife, and that a covariate-free arm assigning each district its region's
# design-based survey estimate outperforms the anchored covariate model on
# correlation, absolute error and bias.
#
# The design question therefore stops being "how much survey does the anchoring
# correction need" and becomes the more useful one: **how thinly can a region be
# sampled before its own mean stops being a good estimate for its districts?**
# That is the quantity a survey planner controls, and on WS2's evidence it is the
# quantity that carries the information.
#
# THE GRID
#   regions   how many of a country's Admin-1 regions are sampled at all.
#             Districts in an unsampled region fall back to the national mean
#             computed from the sampled regions.
#   fraction  what share of each sampled region's clusters is retained.
#
# Whole CLUSTERS are dropped, not individuals, because a survey planner buys
# clusters. Retaining a fraction of individuals within every cluster would price
# a design nobody can field.
#
# The anchor is always a jackknife: a district's regional mean is computed from
# the region's OTHER districts' retained clusters. Without that this script
# would measure the same circularity WS2 found in Section 4.
#
# Scored against the FULL survey's district estimates, which is the best
# available stand-in for the truth.
#
#   PROFILE=smoke   Ghana only, 5 replicates
#   Rscript scripts/accuracy_impact/ws5_anchoring_budget.R
# -> results/tables/anchoring_design_curve.csv
# -> results/figures/anchoring_design_curve.png
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full"); SUF <- if (PROFILE == "smoke") "_SMOKE" else ""
REPS <- as.integer(Sys.getenv("WS5_REPS", if (PROFILE == "smoke") "5" else "25"))
SEED <- 20260908L
FRACTIONS <- c(0.25, 0.5, 0.75, 1.0)
TDIR <- here("results", "tables"); FDIR <- here("results", "figures")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

cfgs <- get_country_configs()
if (PROFILE == "smoke") cfgs <- cfgs["Ghana"]
rows <- list()
set.seed(SEED)

for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", ocn), store = STORE),
                   error = function(e) NULL)
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", ocn), store = STORE),
                   error = function(e) NULL)
    if (is.null(od) || is.null(sv) || is.null(oc$binary)) next
    d <- od$data
    need <- c(cc$admin1_col, cc$admin2_col, cc$cluster_id, cc$weight_col, oc$binary)
    if (!all(need %in% names(d))) next
    y <- num(d[[oc$binary]]); w <- num(d[[cc$weight_col]])
    w[!is.finite(w) | w <= 0] <- 1
    a1 <- trimws(as.character(d[[cc$admin1_col]]))
    a2 <- trimws(as.character(d[[cc$admin2_col]]))
    cl <- as.character(d[[cc$cluster_id]])
    ok <- is.finite(y) & !is.na(a1) & !is.na(a2) & !is.na(cl)
    y <- y[ok]; w <- w[ok]; a1 <- a1[ok]; a2 <- a2[ok]; cl <- cl[ok]
    if (length(y) < 100) next
    regions <- sort(unique(a1))
    if (length(regions) < 3) next

    # The yardstick: the full survey's district and region estimates.
    truth_d <- sv[is.finite(sv$svy_prev), , drop = FALSE]
    key_d <- if (all(c("Admin1","Admin2") %in% names(truth_d)))
      paste(truth_d$Admin1, truth_d$Admin2) else truth_d$Admin2
    truth_a1 <- vapply(regions, function(r)
      stats::weighted.mean(y[a1 == r], w[a1 == r]), numeric(1))
    dist_reg <- tapply(a1, paste(a1, a2), function(z) z[1])

    n_reg_grid <- unique(pmax(2L, round(length(regions) *
                                          c(0.25, 0.5, 0.75, 1.0))))
    for (nr in n_reg_grid) for (fr in FRACTIONS) for (rep in seq_len(REPS)) {
      keep_reg <- sample(regions, nr)
      # retain a fraction of each kept region's clusters
      keepc <- rep(FALSE, length(y))
      for (r in keep_reg) {
        cu <- unique(cl[a1 == r])
        k <- max(1L, round(length(cu) * fr))
        sel <- sample(cu, k)
        keepc[a1 == r & cl %in% sel] <- TRUE
      }
      if (sum(keepc) < 30) next
      nat <- stats::weighted.mean(y[keepc], w[keepc])

      # Jackknife regional anchor: a district's own respondents are excluded.
      pred <- rep(NA_real_, nrow(truth_d))
      for (i in seq_len(nrow(truth_d))) {
        rg <- if ("Admin1" %in% names(truth_d)) truth_d$Admin1[i] else
          dist_reg[[paste(truth_d$Admin2[i])]] %||% NA
        dd <- truth_d$Admin2[i]
        if (is.na(rg) || !(rg %in% keep_reg)) { pred[i] <- nat; next }
        j <- which(keepc & a1 == rg & a2 != dd)
        pred[i] <- if (length(j) >= 5) stats::weighted.mean(y[j], w[j]) else nat
      }
      fin <- is.finite(pred) & is.finite(truth_d$svy_prev)
      if (sum(fin) < 5) next

      # Admin-1: the region's retained-cluster mean against its full-survey mean.
      pa <- vapply(regions, function(r) {
        if (!(r %in% keep_reg)) return(nat)
        j <- which(keepc & a1 == r)
        if (length(j) >= 5) stats::weighted.mean(y[j], w[j]) else nat
      }, numeric(1))

      rows[[length(rows) + 1L]] <- data.frame(
        country = cc$country, outcome = ocn, replicate = rep,
        n_regions_total = length(regions), n_regions_anchored = nr,
        fraction_clusters = fr,
        clusters_retained = length(unique(cl[keepc])),
        n_retained = sum(keepc),
        pct_survey_used = round(100 * sum(keepc) / length(y), 1),
        mae_admin2_pp = round(100 * mean(abs(pred[fin] - truth_d$svy_prev[fin])), 3),
        bias_admin2_pp = round(100 * mean(pred[fin] - truth_d$svy_prev[fin]), 3),
        mae_admin1_pp = round(100 * mean(abs(pa - truth_a1)), 3),
        bias_admin1_pp = round(100 * mean(pa - truth_a1), 3),
        stringsAsFactors = FALSE)
    }
    cat(sprintf("  [ok] %-13s %-13s regions=%d\n", cn, ocn, length(regions)))
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
readr::write_csv(res, file.path(TDIR, sprintf("anchoring_design_curve%s.csv", SUF)))

s <- res |> group_by(n_regions_anchored, fraction_clusters) |>
  summarise(settings = dplyr::n(),
            pct_survey = round(mean(pct_survey_used), 1),
            mae_a2 = round(mean(mae_admin2_pp), 3),
            mae_a2_lo = round(stats::quantile(mae_admin2_pp, .1), 3),
            mae_a2_hi = round(stats::quantile(mae_admin2_pp, .9), 3),
            bias_a2 = round(mean(bias_admin2_pp), 3),
            mae_a1 = round(mean(mae_admin1_pp), 3),
            .groups = "drop")
readr::write_csv(s, file.path(TDIR, sprintf("anchoring_design_summary%s.csv", SUF)))
cat("\n=== WS5: error against anchoring budget ===\n")
print(as.data.frame(s), row.names = FALSE)

# The reference point: every region anchored, every cluster retained.
full <- s[s$fraction_clusters == 1 & s$n_regions_anchored == max(s$n_regions_anchored), ]
if (nrow(full)) {
  ref <- full$mae_a2[1]
  tol <- 1.0   # percentage points
  ok <- s[s$mae_a2 <= ref + tol, ]
  ok <- ok[order(ok$pct_survey), ]
  cat(sprintf("\nfull-anchor Admin-2 MAE: %.3f pp\n", ref))
  cat(sprintf("smallest budget within %.1f pp of it: %.1f%% of the survey ",
              tol, ok$pct_survey[1]))
  cat(sprintf("(%d of %d regions, %.0f%% of clusters), MAE %.3f pp [%.3f, %.3f]\n",
              ok$n_regions_anchored[1], max(s$n_regions_anchored),
              100 * ok$fraction_clusters[1], ok$mae_a2[1],
              ok$mae_a2_lo[1], ok$mae_a2_hi[1]))
}

png(file.path(FDIR, sprintf("anchoring_design_curve%s.png", SUF)),
    width = 1100, height = 800, res = 130)
op <- par(mar = c(4.5, 4.5, 3, 1))
plot(s$pct_survey, s$mae_a2, type = "n",
     xlab = "Percent of the survey retained",
     ylab = "Admin-2 mean absolute error (pp)",
     main = "Error against anchoring budget")
frs <- sort(unique(s$fraction_clusters))
cols <- grDevices::hcl.colors(length(frs), "Dark 3")
for (i in seq_along(frs)) {
  z <- s[s$fraction_clusters == frs[i], ]
  z <- z[order(z$pct_survey), ]
  lines(z$pct_survey, z$mae_a2, col = cols[i], lwd = 2)
  points(z$pct_survey, z$mae_a2, col = cols[i], pch = 19)
}
if (nrow(full)) abline(h = full$mae_a2[1], lty = 2, col = "grey50")
legend("topright", bty = "n", lwd = 2, col = cols,
       legend = sprintf("%.0f%% of clusters", 100 * frs))
par(op); dev.off()
cat(sprintf("-> %s\n", file.path("results", "figures",
                                 sprintf("anchoring_design_curve%s.png", SUF))))
