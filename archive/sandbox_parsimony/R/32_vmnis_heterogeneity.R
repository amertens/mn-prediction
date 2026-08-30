# =============================================================================
# sandbox_parsimony/R/32_vmnis_heterogeneity.R
#
# Model the VMNIS survey heterogeneity explicitly, and rank the national
# predictors by importance.
#
# FINDINGS.md section 23 found the national model already at ~98% of its ceiling
# (r = 0.648 vs r_max = 0.66), where the ceiling is set by how much repeat
# national surveys of the SAME country disagree. Part of that disagreement is
# genuine temporal change, but part is measurement: VMNIS records the
# inflammation adjustment, the biomarker and the sample size per survey, and the
# raw means differ a lot by method --
#
#   Not Specified  n=290  25.3%   logit -1.46
#   Inflamation    n=177  21.4%   logit -1.67
#   None           n= 42  12.9%   logit -2.16
#
# CAUTION, and it drives the design below: that raw cross-tabulation confounds
# method with country and era -- countries that ran inflammation-adjusted
# surveys are not a random sample. A method effect is only identified WITHIN a
# country, or with country held constant. So the method term is estimated with
# country effects in the model, not from the marginal table.
#
# Three things are done here:
#   A. does adding the heterogeneity terms improve honest LOCO prediction?
#   B. does the CEILING rise once method variance is taken out of the
#      within-country disagreement that defines it?
#   C. permutation importance of the national predictors.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ranger); library(glmnet)})

`%|%` <- function(a, b) if (is.null(a)) b else a
.logit  <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.ilogit <- stats::plogis

nat <- readRDS("sandbox_parsimony/out/vmnis_national_rep.rds")
wdi <- readRDS("sandbox_parsimony/out/wdi_national.rds") |>
  filter(!is.na(iso3c), region != "Aggregates")

WDI_VARS <- c(gdp_pc = "NY.GDP.PCAP.PP.KD", undernourish = "SN.ITK.DEFC.ZS",
              food_prod_idx = "AG.PRD.FOOD.XD", cereal_yield = "AG.YLD.CREL.KG",
              ag_land_pct = "AG.LND.AGRI.ZS", urban_pct = "SP.URB.TOTL.IN.ZS",
              electricity = "EG.ELC.ACCS.ZS", water_basic = "SH.H2O.BASW.ZS",
              sanitation = "SH.STA.BASS.ZS", u5_mortality = "SH.DYN.MORT",
              life_exp = "SP.DYN.LE00.IN", dtp3 = "SH.IMM.IDPT",
              health_exp_pc = "SH.XPD.CHEX.PC.CD", fertility = "SP.DYN.TFRT.IN",
              female_school = "SE.SEC.ENRR.FE", pop_total = "SP.POP.TOTL",
              rural_pop_pct = "SP.RUR.TOTL.ZS")
V <- names(WDI_VARS)

fill_near <- function(df, vars, max_gap = 5L) {
  df |> arrange(iso3c, year) |> group_by(iso3c) |>
    mutate(across(all_of(vars), function(x) {
      if (all(is.na(x))) return(x)
      idx <- which(!is.na(x))
      vapply(seq_along(x), function(i) {
        if (!is.na(x[i])) return(x[i])
        j <- idx[which.min(abs(idx - i))]
        if (abs(j - i) <= max_gap) x[j] else NA_real_
      }, numeric(1))
    })) |> ungroup()
}
wdi <- fill_near(wdi, V)

# --- heterogeneity terms recorded by VMNIS ---------------------------------
# Collapsed to a small set of levels; the free-text Surveymethodology field is
# unusable as a factor (almost every value unique) and is dropped.
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
         biomarker = as.character(Indicator),
         log_n = log1p(suppressWarnings(as.numeric(Samplesize)))) |>
  filter(!is.na(iso3c)) |>
  inner_join(wdi |> select(iso3c, year, all_of(V)), by = c("iso3c", "year"))

# ===========================================================================
# B. Does the ceiling rise once method is accounted for?
#    r_max is defined by within-country disagreement across repeat surveys.
#    Recompute it after removing an estimated method effect, which is
#    identified WITHIN country (country fixed effects), not from the marginal.
# ===========================================================================
ceiling_for <- function(d, resid_method = FALSE) {
  d <- d |> mutate(lg = .logit(prev))
  if (resid_method) {
    # country + method fixed effects: the method contrast is identified only
    # from countries that ran surveys under more than one method.
    fit <- tryCatch(stats::lm(lg ~ factor(iso3c) + adj + biomarker, data = d),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- stats::coef(fit)
      meth <- cf[grepl("^adj|^biomarker", names(cf))]
      meth[!is.finite(meth)] <- 0
      mm <- stats::model.matrix(~ adj + biomarker, data = d)
      keep <- intersect(colnames(mm), names(meth))
      d$lg <- d$lg - as.numeric(mm[, keep, drop = FALSE] %*% meth[keep])
    }
  }
  rep_c <- d |> group_by(iso3c) |> filter(n() >= 2) |>
    summarise(sd_w = stats::sd(lg), .groups = "drop")
  vw <- mean(rep_c$sd_w^2, na.rm = TRUE); vb <- stats::var(d$lg)
  lam <- max(0, 1 - vw / vb)
  list(lambda = lam, r_max = sqrt(lam), n_rep = nrow(rep_c),
       within_sd = sqrt(vw), total_sd = sqrt(vb))
}

PANELS <- list(c("Vitamin A", "Preschool-age children"),
               c("Zinc", "Preschool-age children"),
               c("Folate", "Non-pregnant women (NPW)"))

cat("\n=== B. Ceiling before and after removing survey-method effects ===\n")
cat("method effects estimated WITHIN country (country fixed effects).\n\n")
cb <- list()
for (pp in PANELS) {
  d <- dat |> filter(mn_group == pp[1], pop == pp[2]) |>
    group_by(iso3c, year, adj, biomarker) |>
    summarise(prev = mean(prev), .groups = "drop")
  if (n_distinct(d$iso3c) < 10) next
  a <- ceiling_for(d, FALSE); b <- ceiling_for(d, TRUE)
  cb[[paste(pp, collapse = "/")]] <- data.frame(
    panel = paste(pp, collapse = " / "), surveys = nrow(d),
    countries = n_distinct(d$iso3c), repeated = a$n_rep,
    within_sd_raw = round(a$within_sd, 2), r_max_raw = round(a$r_max, 3),
    within_sd_adj = round(b$within_sd, 2), r_max_adj = round(b$r_max, 3),
    stringsAsFactors = FALSE)
}
cb <- bind_rows(cb); print(as.data.frame(cb), row.names = FALSE)

# ===========================================================================
# A. Does modelling heterogeneity improve honest LOCO prediction?
#    Three specifications, same folds, same learner:
#      wdi_only      national covariates + year
#      wdi_plus_het  + adjustment method, biomarker, log sample size
#      wdi_std       same fit, but PREDICTED at a standardised method
#                    (inflammation-adjusted retinol), i.e. "what would a
#                    BRINDA-style survey have found" -- the quantity a user
#                    actually wants, though it is scored against surveys that
#                    used whatever method they used, so it is penalised here.
# ===========================================================================
fit_panel <- function(d, label) {
  d <- d |> group_by(iso3c, year, adj, biomarker) |>
    summarise(prev = mean(prev), log_n = mean(log_n, na.rm = TRUE),
              across(all_of(V), mean), .groups = "drop") |> as.data.frame()
  d <- d[stats::complete.cases(d[, c("prev", V)]), , drop = FALSE]
  if (n_distinct(d$iso3c) < 12) return(NULL)
  d$log_n[!is.finite(d$log_n)] <- stats::median(d$log_n[is.finite(d$log_n)])

  mk <- function(with_het, std = FALSE) {
    X <- as.data.frame(d[, V, drop = FALSE])
    for (v in c("gdp_pc", "pop_total", "health_exp_pc", "cereal_yield"))
      X[[v]] <- log1p(pmax(X[[v]], 0))
    X$year <- d$year
    if (with_het) {
      a <- if (std) rep("inflammation_adjusted", nrow(d)) else d$adj
      b <- if (std) rep("Retinol (plasma or serum)", nrow(d)) else d$biomarker
      X$adj_infl <- as.numeric(a == "inflammation_adjusted")
      X$adj_unadj <- as.numeric(a == "unadjusted")
      X$adj_unspec <- as.numeric(a == "unspecified")
      X$bio_retinol <- as.numeric(b == "Retinol (plasma or serum)")
      X$log_n <- if (std) stats::median(d$log_n) else d$log_n
    }
    as.matrix(X)
  }

  y <- .logit(d$prev); ctys <- unique(d$iso3c)
  specs <- list(wdi_only = list(het = FALSE, std = FALSE),
                wdi_plus_het = list(het = TRUE, std = FALSE),
                wdi_std = list(het = TRUE, std = TRUE))
  out <- list()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    Xf <- mk(sp$het, FALSE)          # training always uses the TRUE method
    Xp <- mk(sp$het, sp$std)         # prediction may standardise it
    pr <- rep(NA_real_, nrow(d))
    for (ho in ctys) {
      te <- which(d$iso3c == ho); tr <- setdiff(seq_len(nrow(d)), te)
      if (length(tr) < 10) next
      rf <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = y[tr], Xf[tr, , drop = FALSE]),
                                    num.trees = 800, min.node.size = 5, seed = 1),
                     error = function(e) NULL)
      if (!is.null(rf))
        pr[te] <- .ilogit(stats::predict(
          rf, data = data.frame(Xp[te, , drop = FALSE]))$predictions)
    }
    ok <- is.finite(pr)
    if (sum(ok) < 10) next
    out[[nm]] <- data.frame(panel = label, spec = nm, n = sum(ok),
      mae_pp = round(100 * mean(abs(pr[ok] - d$prev[ok])), 2),
      pearson = round(suppressWarnings(cor(pr[ok], d$prev[ok])), 3),
      spearman = round(suppressWarnings(cor(pr[ok], d$prev[ok], method = "spearman")), 3),
      stringsAsFactors = FALSE)
  }
  list(res = bind_rows(out), data = d, y = y, X = mk(TRUE, FALSE))
}

cat("\n=== A. Does modelling survey heterogeneity help LOCO prediction? ===\n")
res <- list(); keep <- list()
for (pp in PANELS) {
  d <- dat |> filter(mn_group == pp[1], pop == pp[2])
  lab <- paste(pp, collapse = " / ")
  f <- tryCatch(fit_panel(d, lab), error = function(e) NULL)
  if (is.null(f) || !nrow(f$res)) next
  res[[lab]] <- f$res; keep[[lab]] <- f
}
rr <- bind_rows(res); print(as.data.frame(rr), row.names = FALSE)

# ===========================================================================
# C. Which national predictors matter?
# ===========================================================================
cat("\n=== C. Permutation importance of national predictors ===\n")
imp_rows <- list()
for (lab in names(keep)) {
  f <- keep[[lab]]
  df <- data.frame(y = f$y, f$X)
  m <- matrix(NA_real_, 10, ncol(f$X), dimnames = list(NULL, colnames(f$X)))
  for (i in 1:10) {
    rf <- tryCatch(ranger::ranger(y ~ ., data = df, num.trees = 1000,
                                  min.node.size = 5, importance = "permutation",
                                  seed = i), error = function(e) NULL)
    if (!is.null(rf)) m[i, names(rf$variable.importance)] <- rf$variable.importance
  }
  mi <- colMeans(m, na.rm = TRUE); mi[!is.finite(mi)] <- 0
  imp_rows[[lab]] <- data.frame(panel = lab, variable = names(mi),
                                imp = as.numeric(mi),
                                rel = as.numeric(mi) / max(mi[mi > 0]),
                                stringsAsFactors = FALSE)
}
imp <- bind_rows(imp_rows)
write.csv(imp, "sandbox_parsimony/out/national_var_importance.csv", row.names = FALSE)

cat("\nAveraged over panels (rel = importance relative to the top variable):\n")
print(as.data.frame(imp |> group_by(variable) |>
  summarise(panels = n(), mean_rel = round(mean(pmax(rel, 0)), 3), .groups = "drop") |>
  arrange(desc(mean_rel)) |> head(20)), row.names = FALSE)

cat("\nTop 5 per panel:\n")
print(as.data.frame(imp |> group_by(panel) |> slice_max(rel, n = 5) |>
  summarise(top5 = paste(variable, collapse = ", "), .groups = "drop")), row.names = FALSE)

write.csv(rr, "sandbox_parsimony/out/vmnis_heterogeneity.csv", row.names = FALSE)
write.csv(cb, "sandbox_parsimony/out/vmnis_ceiling_adjusted.csv", row.names = FALSE)
message("\nSaved -> out/vmnis_heterogeneity.csv, out/vmnis_ceiling_adjusted.csv, out/national_var_importance.csv")
