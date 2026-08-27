# =============================================================================
# sandbox_parsimony/R/18_mixed_level_check.R
#
# QUESTION: is it appropriate to use Admin-1 for some countries and Admin-2 for
# others, so that every country contributes a similar NUMBER of units?
#
# It is testable, and it turns on two things.
#
#  (1) MAUP. Does the covariate-outcome relationship survive a change of
#      aggregation level? If the correlation between a predictor and deficiency
#      measured at Admin-1 differs from the same correlation at Admin-2 -- in
#      size or in sign -- then a pooled model fitted on mixed levels is fitting
#      a relationship that does not exist at either level. This is the modifiable
#      areal unit problem and it is the substantive objection.
#
#  (2) CONFOUNDING OF LEVEL WITH COUNTRY. Admin-1 units are ~3x better measured
#      than Admin-2 units. If Gambia contributes noisy Admin-2 rows while Ghana
#      contributes clean Admin-1 rows, then country, reliability and level are
#      perfectly confounded. A LOCO model would then look like it "fails to
#      transport" when it is only meeting a different noise level.
#
# Tested here:
#   A. per-variable correlation with the outcome at Admin-1 vs Admin-2, and the
#      regression slope, to see whether the relationship is stable across levels
#   B. reliability by level per country, to size the confounding
#   C. a pooled LOCO run three ways: all-Admin-2, all-Admin-1 (where possible),
#      and MIXED (Admin-1 for the countries with many regions, Admin-2 for the
#      small ones), so the practical cost is measured rather than argued
# =============================================================================
.ASSEMBLE_FNS_ONLY <- TRUE
source("sandbox_parsimony/R/00_assemble.R")
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
source("sandbox_parsimony/R/05a_loco_fns.R")
suppressPackageStartupMessages({library(dplyr)})

STORE <- "_targets_full/objects"
COUNTRIES <- c("gambia", "ghana", "sierraleone", "malawi", "tanzania")
CFG_KEY <- c(gambia = "Gambia", ghana = "Ghana", sierraleone = "SierraLeone",
             malawi = "Malawi", tanzania = "Tanzania")
WATER <- "^Lake |^Water|Reservoir|^Lac |^Sea$"
`%|%` <- function(a, b) if (is.null(a)) b else a

sink(stdout(), type = "message")
source("R/config.R")
cfg <- get_country_configs()
cov_ext <- readRDS("sandbox_parsimony/out/cov_ext.rds")
source("sandbox_parsimony/R/17_admin1_layer_fns.R")   # covariates_admin1(), outcome_at()

build <- function(ctry, oc, cc, level) {
  cv <- if (level == "Admin1") covariates_admin1(cov_ext[[ctry]]) else cov_ext[[ctry]]
  if (is.null(cv)) return(NULL)
  sv <- outcome_at(ctry, oc, cc, level); if (is.null(sv)) return(NULL)
  d <- dplyr::inner_join(sv, cv, by = "Admin2")
  d <- d[is.finite(d$svy_prev) & is.finite(d$n_svy) & is.finite(d$lon), , drop = FALSE]
  if (nrow(d) < 5) return(NULL)
  d$country <- ctry; d$level <- level
  d
}

# ---------------------------------------------------------------------------
# A. Does the covariate-outcome relationship survive the change of level?
# ---------------------------------------------------------------------------
maup <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  for (oc in cc$outcomes) {
    d1 <- build(ctry, oc, cc, "Admin1"); d2 <- build(ctry, oc, cc, "Admin2")
    if (is.null(d1) || is.null(d2) || nrow(d1) < 6 || nrow(d2) < 20) next
    cur <- intersect(curated_vars(names(d2)), names(d1))
    if (length(cur) < 5) next
    for (v in cur) {
      x1 <- suppressWarnings(as.numeric(d1[[v]])); x2 <- suppressWarnings(as.numeric(d2[[v]]))
      if (stats::sd(x1, na.rm = TRUE) < 1e-9 || stats::sd(x2, na.rm = TRUE) < 1e-9) next
      r1 <- suppressWarnings(cor(x1, d1$svy_prev, use = "complete.obs"))
      r2 <- suppressWarnings(cor(x2, d2$svy_prev, use = "complete.obs"))
      # slopes on standardised x so they are comparable across variables
      b1 <- suppressWarnings(stats::coef(stats::lm(d1$svy_prev ~ scale(x1)))[2])
      b2 <- suppressWarnings(stats::coef(stats::lm(d2$svy_prev ~ scale(x2)))[2])
      maup[[paste(ctry, oc$tag, v)]] <- data.frame(
        country = ctry, outcome = oc$tag, variable = v,
        n1 = nrow(d1), n2 = nrow(d2),
        r_admin1 = round(r1, 3), r_admin2 = round(r2, 3),
        slope_admin1_pp = round(100 * b1, 2), slope_admin2_pp = round(100 * b2, 2),
        sign_flip = is.finite(r1) && is.finite(r2) && sign(r1) != sign(r2),
        stringsAsFactors = FALSE)
    }
  }
}
mp <- bind_rows(maup)
write.csv(mp, "sandbox_parsimony/out/maup_check.csv", row.names = FALSE)

cat("\n=== A. Same covariate, same outcome, two aggregation levels ===\n")
ok <- is.finite(mp$r_admin1) & is.finite(mp$r_admin2)
cat(sprintf("  %d variable x country x outcome comparisons\n", sum(ok)))
cat(sprintf("  correlation between r(Admin-1) and r(Admin-2) across them: %.3f\n",
            cor(mp$r_admin1[ok], mp$r_admin2[ok])))
cat(sprintf("  sign flips between levels: %d of %d (%.0f%%)\n",
            sum(mp$sign_flip[ok]), sum(ok), 100 * mean(mp$sign_flip[ok])))
cat(sprintf("  median |r| at Admin-1 = %.3f, at Admin-2 = %.3f\n",
            median(abs(mp$r_admin1[ok])), median(abs(mp$r_admin2[ok]))))
cat(sprintf("  median |slope| per SD, Admin-1 = %.2f pp, Admin-2 = %.2f pp (ratio %.2f)\n",
            median(abs(mp$slope_admin1_pp[ok])), median(abs(mp$slope_admin2_pp[ok])),
            median(abs(mp$slope_admin1_pp[ok])) / median(abs(mp$slope_admin2_pp[ok]))))

cat("\n  per outcome:\n")
print(as.data.frame(mp[ok, ] |> group_by(outcome) |>
  summarise(comparisons = n(),
            r_between_levels = round(cor(r_admin1, r_admin2), 3),
            pct_sign_flip = round(100 * mean(sign_flip)),
            .groups = "drop")), row.names = FALSE)

# ---------------------------------------------------------------------------
# B. How badly would level be confounded with country?
# ---------------------------------------------------------------------------
rel <- list()
for (ctry in COUNTRIES) {
  cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
  for (oc in cc$outcomes) for (lv in c("Admin1", "Admin2")) {
    d <- build(ctry, oc, cc, lv); if (is.null(d)) next
    r <- reliability(d$svy_prev, d$n_svy, 1.5)
    rel[[paste(ctry, oc$tag, lv)]] <- data.frame(
      country = ctry, outcome = oc$tag, level = lv, n_units = nrow(d),
      median_n = stats::median(d$n_svy), lambda = round(r$lambda, 3),
      r_max = round(r$r_max, 3), stringsAsFactors = FALSE)
  }
}
rl <- bind_rows(rel)
write.csv(rl, "sandbox_parsimony/out/level_reliability_ext.csv", row.names = FALSE)

cat("\n=== B. Units and reliability by country and level ===\n")
print(as.data.frame(rl |> group_by(country, level) |>
  summarise(outcomes = n(), units = round(median(n_units)),
            median_n_per_unit = round(median(median_n)),
            mean_r_max = round(mean(r_max), 3), .groups = "drop") |>
  arrange(country, level)), row.names = FALSE)

# ---------------------------------------------------------------------------
# C. Pooled LOCO three ways
#    mixed = Admin-1 where a country has >= 15 regions, Admin-2 otherwise --
#    i.e. exactly the "keep the unit count comparable" proposal.
# ---------------------------------------------------------------------------
MIXED_RULE <- function(ctry, n1) if (is.finite(n1) && n1 >= 15) "Admin1" else "Admin2"

pool_for <- function(oc_tag, scheme) {
  frames <- list()
  for (ctry in COUNTRIES) {
    cc <- cfg[[CFG_KEY[[ctry]]]]; if (is.null(cc)) next
    oc <- Filter(function(o) o$tag == oc_tag, cc$outcomes)
    if (!length(oc)) next
    oc <- oc[[1]]
    d1 <- build(ctry, oc, cc, "Admin1")
    lv <- switch(scheme,
                 admin2 = "Admin2",
                 admin1 = "Admin1",
                 mixed  = MIXED_RULE(ctry, if (is.null(d1)) NA else nrow(d1)))
    d <- if (lv == "Admin1") d1 else build(ctry, oc, cc, "Admin2")
    if (is.null(d) || nrow(d) < 6) next
    frames[[ctry]] <- d
  }
  if (length(frames) < 3) return(NULL)
  covn <- function(df) {
    cv <- setdiff(names(df), c("Admin2", "Admin1", "svy_prev", "n_svy", "country",
                               "level", "lon", "lat", "n_admin2", ".w"))
    cv[vapply(cv, function(c) is.numeric(df[[c]]) && any(is.finite(df[[c]])), logical(1))]
  }
  common <- Reduce(intersect, lapply(frames, covn))
  if (length(common) < 10) return(NULL)
  list(data = bind_rows(lapply(frames, function(df)
         df[, c("country", "Admin2", "svy_prev", "n_svy", "lon", "lat", common)])),
       predictors = common)
}

loco_rows <- list()
for (oc_tag in c("child_vitA", "child_iron", "women_iron", "women_vitA")) {
  for (scheme in c("admin2", "admin1", "mixed")) {
    P <- pool_for(oc_tag, scheme); if (is.null(P)) next
    d <- P$data
    for (ho in unique(d$country)) {
      tr <- d[d$country != ho, ]; te <- d[d$country == ho, ]
      if (nrow(tr) < 20 || nrow(te) < 5) next
      pr <- tryCatch(loco_one(tr, te, P$predictors,
                              function(X, y, v) intersect(curated_vars(v), v),
                              alpha = 0, centering = "own", anchor = "train_mean"),
                     error = function(e) NULL)
      if (is.null(pr)) next
      m <- loco_metrics(pr, te); if (is.null(m)) next
      m$outcome <- oc_tag; m$scheme <- scheme; m$held_out <- ho
      m$n_train <- nrow(tr)
      loco_rows[[paste(oc_tag, scheme, ho)]] <- m
    }
  }
}
lc <- bind_rows(loco_rows)
write.csv(lc, "sandbox_parsimony/out/mixed_level_loco.csv", row.names = FALSE)

cat("\n=== C. Pooled LOCO under three unit schemes ===\n")
print(as.data.frame(lc |> group_by(scheme) |>
  summarise(cells = n(), countries = n_distinct(held_out),
            median_n_test = round(median(n_test)),
            rho = round(mean(spearman, na.rm = TRUE), 3),
            r = round(mean(pearson, na.rm = TRUE), 3),
            rmse_pp = round(mean(rmse_pp, na.rm = TRUE), 1), .groups = "drop")),
  row.names = FALSE)
cat("\n  NB the three schemes are NOT scored on the same targets -- an Admin-1\n")
cat("  prevalence is an easier thing to predict than an Admin-2 one. Read this\n")
cat("  as 'what does each product look like', not as a fair contest.\n")
message("\nSaved -> out/maup_check.csv, out/level_reliability_ext.csv, out/mixed_level_loco.csv")
