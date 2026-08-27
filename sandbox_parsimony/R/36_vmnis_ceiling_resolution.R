# =============================================================================
# sandbox_parsimony/R/36_vmnis_ceiling_resolution.R
#
# Reconcile the two disagreeing VMNIS ceilings.
#
# THE PROBLEM. FINDINGS section 23 reported r_max = 0.66 for vitamin A and
# section 24 reported 0.505 for the same panel. Nothing about the data changed
# between them -- only how repeat surveys were grouped before the within-country
# variance was computed:
#
#   section 23  one row per country-YEAR         -> method differences inside a
#                                                   country-year are AVERAGED
#                                                   AWAY, i.e. counted as signal
#   section 24  one row per country-year-METHOD  -> those differences are
#                                                   RETAINED, i.e. counted as error
#
# Neither is wrong; they answer different questions, and letting an arbitrary
# grouping decide is the actual defect. Note this is a VMNIS-only issue: the
# Admin-2 ceiling uses the binomial sampling variance and involves no repeat
# surveys at all.
#
# THE FIX. Make the estimand explicit and estimate the variance components
# instead of inferring error from a grouping choice. VMNIS records Samplesize,
# so the sampling term can be computed directly -- exactly as at Admin-2 --
# rather than bundled in with everything else:
#
#   logit(prev)_s = country_c + method_m + year trend + residual_s + sampling_s
#
#   sigma2_country   between-country signal -- what a national model predicts
#   sigma2_method    measurement heterogeneity (inflammation adjustment, biomarker)
#   sigma2_resid     real temporal change + everything unexplained
#   v_bar            mean binomial sampling variance, deff * p(1-p)/(n-1)
#
# Two ceilings then follow, for two DIFFERENT estimands:
#
#   predicting what a survey will REPORT (method unknown in advance):
#     lambda_raw = sigma2_country / (sigma2_country + sigma2_method + sigma2_resid + v_bar)
#
#   predicting the country's STANDARDISED prevalence (one fixed method):
#     lambda_std = sigma2_country / (sigma2_country + sigma2_resid + v_bar)
#
# The second is higher precisely because method variance has been moved out of
# the error term -- which is what section 25's outcome standardisation actually
# buys, now stated as a number rather than a grouping side effect.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(lme4)})

.logit <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
DEFF <- 1.5

nat <- readRDS("sandbox_parsimony/out/vmnis_national_rep.rds")

adj_level <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    grepl("inflam", x, ignore.case = TRUE) ~ "inflammation_adjusted",
    grepl("^none$", x, ignore.case = TRUE) ~ "unadjusted",
    x == "" | grepl("not specified", x, ignore.case = TRUE) ~ "unspecified",
    TRUE ~ "other")
}

dat <- nat |>
  mutate(iso3c = suppressWarnings(countrycode::countrycode(
           country, "country.name", "iso3c", warn = FALSE)),
         pop = as.character(Population),
         adj = adj_level(Dataadjustedfor),
         biomarker = ifelse(grepl("Retinol", as.character(Indicator)),
                            "retinol", "other"),
         method = paste(adj, biomarker, sep = "|"),
         n_svy = suppressWarnings(as.numeric(Samplesize))) |>
  filter(!is.na(iso3c), is.finite(prev))

PANELS <- list(c("Vitamin A", "Preschool-age children"),
               c("Zinc", "Preschool-age children"),
               c("Folate", "Non-pregnant women (NPW)"))

#' The two legacy estimators, reproduced so the reconciliation is checkable.
legacy_ceiling <- function(d, group_method) {
  g <- if (group_method) c("iso3c", "year", "adj", "biomarker") else c("iso3c", "year")
  dd <- d |> group_by(across(all_of(g))) |>
    summarise(prev = mean(prev), .groups = "drop") |>
    mutate(lg = .logit(prev))
  rep_c <- dd |> group_by(iso3c) |> filter(n() >= 2) |>
    summarise(sd_w = stats::sd(lg), .groups = "drop")
  if (!nrow(rep_c)) return(NA_real_)
  vw <- mean(rep_c$sd_w^2, na.rm = TRUE); vb <- stats::var(dd$lg)
  sqrt(max(0, 1 - vw / vb))
}

#' Variance-components decomposition with the sampling term computed directly.
vc_ceiling <- function(d) {
  d <- d |> mutate(lg = .logit(prev))
  # binomial sampling variance on the logit scale, delta method:
  #   Var(logit p) ~ deff / (n p (1-p))
  ok_n <- is.finite(d$n_svy) & d$n_svy > 10
  v_s <- rep(NA_real_, nrow(d))
  pp <- pmin(pmax(d$prev, 0.01), 0.99)
  v_s[ok_n] <- DEFF / (d$n_svy[ok_n] * pp[ok_n] * (1 - pp[ok_n]))
  v_bar <- mean(v_s, na.rm = TRUE)

  if (n_distinct(d$iso3c) < 10) return(NULL)
  m <- tryCatch(
    lme4::lmer(lg ~ year_c + (1 | iso3c) + (1 | method), data = d,
               control = lme4::lmerControl(calc.derivs = FALSE)),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  vc <- as.data.frame(lme4::VarCorr(m))
  gv <- function(g) { v <- vc$vcov[vc$grp == g]; if (length(v)) v[1] else 0 }
  s2_country <- gv("iso3c"); s2_method <- gv("method")
  s2_resid <- stats::sigma(m)^2
  # the residual absorbs sampling error too; subtract it so the components do
  # not double-count, floored at zero when the estimate goes negative
  s2_resid_only <- max(s2_resid - v_bar, 0)

  tot_raw <- s2_country + s2_method + s2_resid_only + v_bar
  tot_std <- s2_country + s2_resid_only + v_bar
  list(n = nrow(d), countries = n_distinct(d$iso3c),
       pct_with_n = round(100 * mean(ok_n), 0),
       s2_country = s2_country, s2_method = s2_method,
       s2_resid = s2_resid_only, v_bar = v_bar,
       r_max_raw = sqrt(s2_country / tot_raw),
       r_max_std = sqrt(s2_country / tot_std))
}

rows <- list()
for (pp in PANELS) {
  d <- dat |> filter(mn_group == pp[1], pop == pp[2])
  if (nrow(d) < 30) next
  d$year_c <- (d$year - mean(d$year, na.rm = TRUE)) / 10
  lab <- paste(pp, collapse = " / ")
  vc <- vc_ceiling(d)
  if (is.null(vc)) next
  rows[[lab]] <- data.frame(
    panel = lab, surveys = vc$n, countries = vc$countries,
    pct_with_samplesize = vc$pct_with_n,
    legacy_by_country_year = round(legacy_ceiling(d, FALSE), 3),
    legacy_by_country_year_method = round(legacy_ceiling(d, TRUE), 3),
    sd_country = round(sqrt(vc$s2_country), 3),
    sd_method = round(sqrt(vc$s2_method), 3),
    sd_resid = round(sqrt(vc$s2_resid), 3),
    sd_sampling = round(sqrt(vc$v_bar), 3),
    r_max_report = round(vc$r_max_raw, 3),
    r_max_standardised = round(vc$r_max_std, 3),
    stringsAsFactors = FALSE)
}

res <- bind_rows(rows)
write.csv(res, "sandbox_parsimony/out/vmnis_ceiling_resolution.csv", row.names = FALSE)

cat("\n=== Variance components of national logit prevalence ===\n")
print(as.data.frame(res |> select(panel, surveys, countries, pct_with_samplesize,
                                  sd_country, sd_method, sd_resid, sd_sampling)),
      row.names = FALSE)

cat("\n=== The two legacy numbers, and what replaces them ===\n")
cat("r_max_report       = ceiling for predicting what a survey WILL REPORT\n")
cat("                     (method unknown in advance, so method variance is error)\n")
cat("r_max_standardised = ceiling for predicting the country's prevalence on ONE\n")
cat("                     fixed method (method variance removed from the error)\n\n")
print(as.data.frame(res |> select(panel, legacy_by_country_year,
                                  legacy_by_country_year_method,
                                  r_max_report, r_max_standardised)),
      row.names = FALSE)

cat("\nReading: the legacy pair brackets the answer because the coarser grouping\n")
cat("silently counts method differences as signal and the finer one counts them\n")
cat("as error. The variance-components version does not have to choose -- it\n")
cat("reports a ceiling per estimand, and quotes the sampling term separately.\n")
message("\nSaved -> sandbox_parsimony/out/vmnis_ceiling_resolution.csv")
