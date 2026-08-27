# =============================================================================
# sandbox_parsimony/R/14_decision_accuracy.R
#
# "Better than the competitors" is not the same as "accurate enough to use".
# This scores the winning model the way a programme would actually consume it:
#
#   who_acc        share of districts put in the right WHO severity band
#   who_acc_null   same, if you just painted the national mean on every district
#   top20_recall   of the worst 20% of districts by survey prevalence, what
#                  share does the model also place in its worst 20%
#   top20_lift     that recall divided by 0.20 (what random targeting gets)
#   mae_pp         mean absolute error in percentage points
#   mae_null_pp    same for the national-mean map
#
# The null is the honest comparator: a programme that has no model still knows
# the national prevalence. A model earns its keep only by beating that.
#
# Sampling noise in the OUTCOME is charged to the model here, exactly as it is
# in production reporting. Where the noise audit says a cell has no signal, the
# "true" district ranking is itself mostly noise, so top20_recall there is
# measuring nothing -- those cells are flagged, not dropped.
# =============================================================================
source("sandbox_parsimony/R/02_features.R")
source("sandbox_parsimony/R/03_core.R")
suppressPackageStartupMessages({library(dplyr); library(ranger)})

pooled_all <- readRDS("sandbox_parsimony/out/pooled_all.rds")
audit <- read.csv("sandbox_parsimony/out/noise_audit.csv", stringsAsFactors = FALSE)

WHO <- list(
  child_vitA   = c(0.02, 0.10, 0.20), women_vitA   = c(0.02, 0.10, 0.20),
  child_iron   = c(0.05, 0.20, 0.40), women_iron   = c(0.05, 0.20, 0.40),
  women_folate = c(0.05, 0.20, 0.40), women_b12    = c(0.05, 0.20, 0.40),
  child_zinc   = c(0.10, 0.20, 0.30), women_zinc   = c(0.10, 0.20, 0.30)
)
who_band <- function(p, oc) {
  th <- WHO[[oc]]; if (is.null(th)) return(rep(NA_character_, length(p)))
  cut(p, breaks = c(-Inf, th, Inf), labels = c("Low", "Mild", "Moderate", "Severe"),
      right = FALSE) |> as.character()
}

OUTCOMES <- names(pooled_all)
MIN_AREAS <- 25L
N_REP <- 10L

rows <- list()
for (oc in OUTCOMES) {
  P <- pooled_all[[oc]]; if (is.null(P)) next
  cur <- curated_vars(P$predictors)
  for (ctry in P$countries) {
    dat <- P$data[P$data$country == ctry, , drop = FALSE]
    dat <- dat[is.finite(dat$svy_prev) & is.finite(dat$n_svy), , drop = FALSE]
    if (nrow(dat) < MIN_AREAS) next

    n <- nrow(dat)
    acc <- list()
    for (rep_i in seq_len(N_REP)) {
      set.seed(700L + rep_i)
      fold <- sample(rep_len(1:5, n))
      oof <- rep(NA_real_, n); null_pred <- rep(NA_real_, n)
      for (f in 1:5) {
        te <- which(fold == f); tr <- setdiff(seq_len(n), te)
        pp <- prep_X(dat[tr, , drop = FALSE], dat[te, , drop = FALSE], P$predictors)
        sel <- intersect(cur, pp$vars)
        Xtr <- cbind(pp$Xtr[, sel, drop = FALSE], lon = dat$lon[tr], lat = dat$lat[tr])
        Xte <- cbind(pp$Xte[, sel, drop = FALSE], lon = dat$lon[te], lat = dat$lat[te])
        df <- data.frame(y = .logit(dat$svy_prev[tr]), Xtr, check.names = TRUE)
        rf <- tryCatch(ranger::ranger(y ~ ., data = df, num.trees = 1000,
                                      min.node.size = 5,
                                      case.weights = pmax(dat$n_svy[tr], 1),
                                      seed = rep_i), error = function(e) NULL)
        oof[te] <- if (is.null(rf)) mean(dat$svy_prev[tr]) else
          .ilogit(stats::predict(rf, data = data.frame(Xte, check.names = TRUE))$predictions)
        # the no-model comparator: the training-fold national mean, survey-weighted
        null_pred[te] <- stats::weighted.mean(dat$svy_prev[tr], pmax(dat$n_svy[tr], 1))
      }
      obs <- dat$svy_prev
      k <- max(2L, round(0.2 * n))
      worst_obs <- order(obs, decreasing = TRUE)[seq_len(k)]
      worst_mod <- order(oof, decreasing = TRUE)[seq_len(k)]
      acc[[rep_i]] <- data.frame(
        who_acc      = mean(who_band(oof, oc) == who_band(obs, oc), na.rm = TRUE),
        who_acc_null = mean(who_band(null_pred, oc) == who_band(obs, oc), na.rm = TRUE),
        top20_recall = length(intersect(worst_obs, worst_mod)) / k,
        mae_pp       = mean(abs(oof - obs)) * 100,
        mae_null_pp  = mean(abs(null_pred - obs)) * 100,
        spearman     = suppressWarnings(cor(oof, obs, method = "spearman")))
    }
    a <- bind_rows(acc)
    rmax <- audit$r_max_d15[audit$outcome == oc & audit$country == ctry]
    rows[[paste(oc, ctry)]] <- data.frame(
      outcome = oc, country = ctry, n_areas = n,
      r_max = if (length(rmax)) rmax else NA_real_,
      obs_prev_pp = round(100 * stats::weighted.mean(dat$svy_prev, pmax(dat$n_svy, 1)), 1),
      who_acc      = round(mean(a$who_acc), 3),
      who_acc_null = round(mean(a$who_acc_null), 3),
      who_gain     = round(mean(a$who_acc) - mean(a$who_acc_null), 3),
      top20_recall = round(mean(a$top20_recall), 3),
      top20_lift   = round(mean(a$top20_recall) / 0.2, 2),
      mae_pp       = round(mean(a$mae_pp), 2),
      mae_null_pp  = round(mean(a$mae_null_pp), 2),
      mae_gain_pp  = round(mean(a$mae_null_pp) - mean(a$mae_pp), 2),
      spearman     = round(mean(a$spearman), 3),
      stringsAsFactors = FALSE)
    message(sprintf("%-14s %-12s WHO %.2f (null %.2f)  top20 recall %.2f (lift %.1fx)  MAE %.1f vs null %.1f",
                    oc, ctry, mean(a$who_acc), mean(a$who_acc_null),
                    mean(a$top20_recall), mean(a$top20_recall)/0.2,
                    mean(a$mae_pp), mean(a$mae_null_pp)))
  }
}

res <- bind_rows(rows)
write.csv(res, "sandbox_parsimony/out/decision_accuracy.csv", row.names = FALSE)

cat("\n=== Is the map good enough to target with? (curated16 + lon/lat RF) ===\n")
print(as.data.frame(res |> arrange(desc(r_max))), row.names = FALSE)

cat("\n=== Averages, split by whether the cell has detectable signal ===\n")
print(as.data.frame(res |> mutate(grp = ifelse(r_max >= 0.35, "signal", "no signal")) |>
  group_by(grp) |> summarise(cells = n(),
    who_acc = round(mean(who_acc), 3), who_null = round(mean(who_acc_null), 3),
    top20_recall = round(mean(top20_recall), 3),
    mae_pp = round(mean(mae_pp), 2), mae_null_pp = round(mean(mae_null_pp), 2),
    spearman = round(mean(spearman), 3), .groups = "drop")), row.names = FALSE)
message("\nSaved -> sandbox_parsimony/out/decision_accuracy.csv")
