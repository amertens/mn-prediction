# =============================================================================
# sandbox_parsimony/R/31_national_model.R
#
# An HONEST national-prevalence model: no benchmarking, no national anchor, no
# survey from the target country. Held-out country never contributes an outcome.
#
# Why this design rather than aggregating the Admin-2 model. FINDINGS.md §19:
# the Admin-2 LOCO models recover national prevalence to ~5 pp mean absolute
# error and miss by 2-3x on iron, because they were built to resolve WITHIN-
# country pattern and are trained on 4 countries. Here the unit of observation
# IS the country-survey, so a 70-country VMNIS panel gives an order of magnitude
# more design points for the between-country relationship that national
# prediction actually depends on.
#
# Covariates: World Bank WDI national indicators, chosen a priori on the same
# causal story as the Admin-2 curated set -- diet quality and food supply,
# poverty and market access, health system, and demography. Lagged/interpolated
# to the survey year.
#
# Evaluation: LEAVE-ONE-COUNTRY-OUT. Every survey from the held-out country is
# removed, so the model has never seen that country's prevalence at any date.
# That is the situation the product is for.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(glmnet); library(ranger)})
source("sandbox_parsimony/R/37_national_covariates.R")

# Covariate source: "wdi" = the original 17 a-priori World Bank indicators,
# "panel" = the project's national country x year panel from mn-proxies,
# "both" = the union. Set with NAT_COV; NAT_COV_MIN sets the panel's minimum
# column coverage and NAT_SCREEN the per-fold prescreen width.
COV_SOURCE <- Sys.getenv("NAT_COV", "wdi")
SCREEN_K   <- as.integer(Sys.getenv("NAT_SCREEN", "150"))
KNN_K      <- as.integer(Sys.getenv("NAT_KNN", "5"))

nat <- readRDS("sandbox_parsimony/out/vmnis_national_rep.rds")

cov <- load_national_covariates(COV_SOURCE)
V        <- cov$vars
LOG_VARS <- cov$log_vars
message(sprintf("covariate source: %s | %d covariates | prescreen to %d per fold",
                cov$source, length(V), SCREEN_K))

# ---------------------------------------------------------------------------
# Join VMNIS to the covariates on ISO3 + survey year
# ---------------------------------------------------------------------------
nat <- nat |>
  mutate(iso3c = suppressWarnings(
    countrycode::countrycode(country, "country.name", "iso3c", warn = FALSE)),
    pop = as.character(Population))

n_cty_vmnis <- n_distinct(nat$iso3c[!is.na(nat$iso3c)])

dat <- nat |>
  filter(!is.na(iso3c)) |>
  inner_join(cov$df |> select(iso3c, year, all_of(V)),
             by = c("iso3c", "year"))

# NAT_COMMON=1 restricts to countries both sources cover, so wdi and panel are
# scored on identical folds and the contrast is about the covariates alone.
if (nzchar(Sys.getenv("NAT_COMMON"))) {
  keep_cty <- common_countries()
  dat <- dat |> filter(iso3c %in% keep_cty)
  message(sprintf("common support: restricted to %d countries", length(keep_cty)))
}

message(sprintf("VMNIS x %s joined: %d rows, %d countries (VMNIS had %d)",
                cov$source, nrow(dat), n_distinct(dat$iso3c), n_cty_vmnis))

# ---------------------------------------------------------------------------
# Model one micronutrient x population panel
# ---------------------------------------------------------------------------
.logit  <- function(p, eps = 0.005) stats::qlogis(pmin(pmax(p, eps), 1 - eps))
.ilogit <- stats::plogis

fit_loco_national <- function(d, vars, label, screen_k = SCREEN_K) {
  # Complete-case on the OUTCOME only. With ~1,600 candidate columns a
  # complete-case rule over the design matrix would drop nearly every survey,
  # so missing covariates are median-imputed inside each fold instead.
  d <- d[is.finite(d$prev), , drop = FALSE]
  if (!nrow(d)) return(NULL)

  .mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  # one row per country-year: average repeated indicators within a survey year
  d <- d |> group_by(iso3c, year) |>
    summarise(prev = mean(prev), n_svy = .mean(Samplesize),
              across(all_of(vars), .mean), .groups = "drop") |>
    as.data.frame()
  ctys <- unique(d$iso3c)
  if (length(ctys) < 12) return(NULL)

  X <- as.matrix(d[, vars, drop = FALSE])
  storage.mode(X) <- "double"
  X[!is.finite(X)] <- NA_real_
  keep <- colSums(!is.na(X)) > 0L &
    apply(X, 2, function(x) length(unique(x[!is.na(x)])) > 1L)
  X <- X[, keep, drop = FALSE]
  for (v in intersect(LOG_VARS, colnames(X)))
    X[, v] <- log1p(pmax(X[, v], 0))
  y <- .logit(d$prev)

  preds <- matrix(NA_real_, nrow(d), 3,
                  dimnames = list(NULL, c("null", "ridge", "rf")))
  n_used <- n_mi <- n_mi_sel <- integer(0)
  for (ho in ctys) {
    te <- which(d$iso3c == ho); tr <- setdiff(seq_len(nrow(d)), te)
    if (length(tr) < 10) next
    preds[te, "null"] <- .ilogit(mean(y[tr]))

    # Missingness indicators from the raw pattern, then KNN imputation with
    # training rows as the only donors, then screening on training rows. All
    # three inside the fold or the held-out country informs the design.
    #
    # Surveys are never dropped for missing covariates: with ~1,600 candidate
    # columns a complete-case rule would discard nearly all of them, and "this
    # country never reported this indicator" is signal in a DHS-derived panel,
    # not noise to be papered over with a median.
    MI  <- missing_indicators(X, tr)
    Xi  <- knn_impute_fold(X, tr, k = KNN_K)
    Xf  <- if (is.null(MI)) Xi else cbind(Xi, MI)
    n_mi <- c(n_mi, if (is.null(MI)) 0L else ncol(MI))

    sel <- prescreen_cols(Xf[tr, , drop = FALSE], y[tr], screen_k)
    Xtr <- Xf[tr, sel, drop = FALSE]; Xte <- Xf[te, sel, drop = FALSE]
    n_used <- c(n_used, ncol(Xtr))
    n_mi_sel <- c(n_mi_sel, sum(grepl("^mi_", colnames(Xtr))))

    cvm <- tryCatch(glmnet::cv.glmnet(Xtr, y[tr], alpha = 0, nfolds = 5),
                    error = function(e) NULL)
    if (!is.null(cvm))
      preds[te, "ridge"] <- .ilogit(as.numeric(
        stats::predict(cvm, newx = Xte, s = "lambda.min")))
    rf <- tryCatch(ranger::ranger(y ~ ., data = data.frame(y = y[tr], Xtr),
                                  num.trees = 800, min.node.size = 5, seed = 1),
                   error = function(e) NULL)
    if (!is.null(rf))
      preds[te, "rf"] <- .ilogit(stats::predict(
        rf, data = data.frame(Xte))$predictions)
  }

  out <- lapply(colnames(preds), function(k) {
    p <- preds[, k]; ok <- is.finite(p)
    if (sum(ok) < 10) return(NULL)
    data.frame(panel = label, source = cov$source, model = k,
               n_cov_pool = ncol(X),
               n_mi_added = if (length(n_mi)) round(mean(n_mi)) else NA_integer_,
               n_cov_used = if (length(n_used)) round(mean(n_used)) else NA_integer_,
               n_mi_selected = if (length(n_mi_sel)) round(mean(n_mi_sel)) else NA_integer_,
               n_surveys = sum(ok),
               n_countries = n_distinct(d$iso3c[ok]),
               mean_prev_pp = round(100 * mean(d$prev[ok]), 1),
               mae_pp = round(100 * mean(abs(p[ok] - d$prev[ok])), 2),
               rmse_pp = round(100 * sqrt(mean((p[ok] - d$prev[ok])^2)), 2),
               pearson = round(suppressWarnings(cor(p[ok], d$prev[ok])), 3),
               spearman = round(suppressWarnings(
                 cor(p[ok], d$prev[ok], method = "spearman")), 3),
               stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

PANELS <- list(
  c("Vitamin A", "Preschool-age children"),
  c("Zinc", "Preschool-age children"),
  c("Folate", "Non-pregnant women (NPW)"),
  c("Vitamin B12", "Non-pregnant women (NPW)"),
  c("Vitamin A", "School-age children"),
  c("Vitamin A", "Non-pregnant women (NPW)")
)

res <- list()
for (pp in PANELS) {
  d <- dat |> filter(mn_group == pp[1], pop == pp[2])
  lab <- paste(pp[1], "|", pp[2])
  r <- tryCatch(fit_loco_national(d, V, lab), error = function(e) NULL)
  if (is.null(r)) { message("  skipped (too few countries): ", lab); next }
  res[[lab]] <- r
  message(sprintf("  %-40s %d surveys / %d countries", lab,
                  r$n_surveys[1], r$n_countries[1]))
}

out <- bind_rows(res)
OUT_CSV <- sprintf("sandbox_parsimony/out/national_model_loco_%s%s.csv", cov$source,
                   if (nzchar(Sys.getenv("NAT_COMMON"))) "_common" else "")
write.csv(out, OUT_CSV, row.names = FALSE)

cat("\n=== Honest leave-one-COUNTRY-out national prediction (no benchmarking) ===\n")
cat("null = global mean of the training countries, the no-covariate comparator.\n\n")
print(as.data.frame(out |> arrange(panel, mae_pp)), row.names = FALSE)

cat("\n=== averaged over panels ===\n")
print(as.data.frame(out |> group_by(model) |>
  summarise(panels = n(), mae_pp = round(mean(mae_pp), 2),
            rmse_pp = round(mean(rmse_pp), 2),
            pearson = round(mean(pearson), 3),
            spearman = round(mean(spearman), 3), .groups = "drop") |>
  arrange(mae_pp)), row.names = FALSE)
message("Saved -> ", OUT_CSV)
