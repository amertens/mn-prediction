# =============================================================================
# scripts/covariates/14_relative_target.R
#
# DOES MODELLING PREVALENCE *RELATIVE TO THE NATIONAL LEVEL* FIX TRANSPORT?
#
# The project's central measured finding is that district RANKING transports
# across countries while the absolute LEVEL does not (mean |national bias|
# 10.6 pp, worst case 42 pp). The natural response is to stop asking the model
# to carry the level: fit on the district's deviation from its own country's
# national prevalence, and put the level back from a design-based national
# estimate at prediction time.
#
#   absolute target   y_i = logit(p_i)
#   relative target   y_i = logit(p_i) - logit(p_nat(country of i))
#
# WHY THIS IS A CROSS-COUNTRY TEST AND NOT A WITHIN-COUNTRY ONE
# -------------------------------------------------------------
# Within one country logit(p_nat) is a CONSTANT, so the relative target is the
# absolute target shifted by a constant. Every learner in the stack carries an
# intercept, so the fitted predictions are algebraically identical and the
# exercise is a no-op. (Verified numerically at the end of this script rather
# than asserted.) The offset only varies across countries -- and it varies a
# lot: child iron national prevalence runs 5.1% in Sierra Leone to 37.5% in
# The Gambia. So this is a leave-one-country-out test.
#
# THE FOUR ARMS, AND WHY THE FOURTH IS NOT OPTIONAL
# --------------------------------------------------
#   1 absolute            pooled fit on logit(p), predict the held-out country
#   2 relative + true     fit on deviations, reconstruct with the HELD-OUT
#                         country's own design-based national prevalence
#   3 relative + blind    reconstruct with the TRAINING countries' mean
#                         national prevalence -- uses nothing from the held-out
#                         country, so it isolates what centring alone buys
#   4 anchor only         predict EVERY district at the held-out country's
#                         national prevalence; no covariates at all
#
# Arm 2 is handed the held-out country's true national level, so its MAE will
# look excellent even if its deviations are pure noise. Arm 4 is the null that
# isolates that: arm 2 is only interesting where it beats arm 4. Reporting 2
# without 4 would be the same "looks strong because it was given the answer"
# error this project has already had to correct twice.
#
#   Rscript scripts/covariates/14_relative_target.R
# -> results/tables/relative_target_loco.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")

OUTCOMES  <- c("child_iron", "child_vitA", "women_iron", "women_vitA")
COUNTRIES <- c(Gambia = "gambia", Ghana = "ghana",
               SierraLeone = "sierraleone", Malawi = "malawi")
K_SCREEN <- 20L

lg  <- function(p) stats::qlogis(pmin(pmax(p, 0.002), 0.998))
inv <- stats::plogis

nat <- readr::read_csv(here("results", "tables", "national_estimates_all.csv"),
                       show_col_types = FALSE)
nat_of <- function(country_label, outcome) {
  v <- nat$obs_prev[nat$country == country_label & nat$outcome == outcome]
  if (length(v)) as.numeric(v[1]) else NA_real_
}

# ── Assemble one admin-2 frame per country x outcome ─────────────────────────
frames <- list()
for (lab in names(COUNTRIES)) {
  lc <- COUNTRIES[[lab]]
  cov <- tryCatch(targets::tar_read_raw(paste0("area_covariates_", lc), store = STORE),
                  error = function(e) NULL)
  if (is.null(cov)) { message("no covariates for ", lab); next }
  # area_covariates_* is a LIST(gee_admin2, polygons), not a data.frame. Reading
  # it straight into a join fails with "x and y must share the same src"; worse,
  # the earlier fallback to gee_admin2_* masked this by handing back a plain
  # frame, so the shape difference only surfaced once the targets were restored.
  if (is.list(cov) && !is.data.frame(cov) && "gee_admin2" %in% names(cov))
    cov <- cov$gee_admin2
  for (oc in OUTCOMES) {
    svy <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", lc, "_", oc), store = STORE),
                    error = function(e) NULL)
    if (is.null(svy) || !nrow(svy)) next
    by <- admin2_join_by(svy, cov)
    d <- dplyr::inner_join(svy[, c(by, "svy_prev", "n_svy")], cov, by = by)
    d <- d[is.finite(d$svy_prev), , drop = FALSE]
    if (nrow(d) < 8) next
    d$country <- lab; d$outcome <- oc
    frames[[paste(lab, oc)]] <- d
  }
}

covs_common <- Reduce(intersect, lapply(frames, function(d)
  names(d)[vapply(d, is.numeric, logical(1))]))
covs_common <- setdiff(covs_common, c("svy_prev", "n_svy"))
cat(sprintf("[rel] %d frames | %d common numeric covariates\n",
            length(frames), length(covs_common)))

prep_X <- function(d) {
  X <- as.matrix(d[, covs_common, drop = FALSE]); X[!is.finite(X)] <- NA
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  X
}

rows <- list()
for (oc in OUTCOMES) {
  fr <- frames[grepl(paste0(" ", oc, "$"), names(frames))]
  if (length(fr) < 2) next
  labs <- vapply(fr, function(d) d$country[1], character(1))
  for (ho in unique(labs)) {
    te <- dplyr::bind_rows(fr[labs == ho])
    tr <- dplyr::bind_rows(fr[labs != ho])
    if (!nrow(te) || nrow(tr) < 20) next
    p_nat_te <- nat_of(if (ho == "SierraLeone") "Sierra Leone" else ho, oc)
    tr_nat <- vapply(tr$country, function(cl)
      nat_of(if (cl == "SierraLeone") "Sierra Leone" else cl, oc), numeric(1))
    if (!is.finite(p_nat_te) || any(!is.finite(tr_nat))) next

    Xtr <- prep_X(tr); Xte <- prep_X(te)
    keep <- apply(Xtr, 2, function(z) stats::sd(z, na.rm = TRUE) > 0)
    Xtr <- Xtr[, keep, drop = FALSE]; Xte <- Xte[, keep, drop = FALSE]
    w_tr <- area_precision_weights(tr$svy_prev, tr$n_svy)

    fit_arm <- function(y) {
      sel <- .awsl_screen(Xtr, y, K_SCREEN)
      s <- tryCatch(.awsl_stack(Xtr[, sel, drop = FALSE], y, w_tr),
                    error = function(e) NULL)
      if (is.null(s)) return(rep(NA_real_, nrow(Xte)))
      .awsl_predict(s, Xte[, sel, drop = FALSE])
    }

    # arm 1 — absolute
    p_abs <- inv(fit_arm(lg(tr$svy_prev)))
    # arms 2/3 — relative: fit deviations, then re-anchor
    dev   <- fit_arm(lg(tr$svy_prev) - lg(tr_nat))
    p_rel_true  <- inv(lg(p_nat_te) + dev)
    p_rel_blind <- inv(lg(mean(tr_nat)) + dev)
    # arm 4 — anchor only, no covariates
    p_anchor <- rep(p_nat_te, nrow(te))

    met <- function(p, arm) {
      ok <- is.finite(p) & is.finite(te$svy_prev)
      data.frame(outcome = oc, held_out = ho, arm = arm, n_test = sum(ok),
        r = if (sum(ok) > 3 && stats::sd(p[ok]) > 0)
              round(suppressWarnings(stats::cor(te$svy_prev[ok], p[ok])), 3) else NA_real_,
        mae_pp  = round(100 * mean(abs(te$svy_prev[ok] - p[ok])), 2),
        bias_pp = round(100 * (mean(p[ok]) - mean(te$svy_prev[ok])), 2),
        stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1L]] <- rbind(
      met(p_abs, "1 absolute"), met(p_rel_true, "2 relative + true anchor"),
      met(p_rel_blind, "3 relative + blind anchor"), met(p_anchor, "4 anchor only"))
    cat(sprintf("  %-11s %-11s abs r=%.2f mae=%.1f | rel r=%.2f mae=%.1f | anchor mae=%.1f\n",
                oc, ho, met(p_abs,"")$r, met(p_abs,"")$mae_pp,
                met(p_rel_true,"")$r, met(p_rel_true,"")$mae_pp,
                met(p_anchor,"")$mae_pp))
  }
}

res <- dplyr::bind_rows(rows)
readr::write_csv(res, here("results", "tables", "relative_target_loco.csv"))

cat("\n================ RELATIVE vs ABSOLUTE TARGET, LOCO ====================\n")
print(res %>% group_by(arm) %>%
        summarise(cells = n(), mean_r = round(mean(r, na.rm = TRUE), 3),
                  med_r = round(stats::median(r, na.rm = TRUE), 3),
                  mean_mae = round(mean(mae_pp, na.rm = TRUE), 2),
                  mean_abs_bias = round(mean(abs(bias_pp), na.rm = TRUE), 2),
                  .groups = "drop") %>% as.data.frame(), row.names = FALSE)

w <- res %>% select(outcome, held_out, arm, mae_pp) %>%
  tidyr::pivot_wider(names_from = arm, values_from = mae_pp)
cat(sprintf("\nrelative(true anchor) beats absolute:    %d/%d cells\n",
            sum(w$`2 relative + true anchor` < w$`1 absolute`, na.rm = TRUE), nrow(w)))
cat(sprintf("relative(true anchor) beats ANCHOR ONLY: %d/%d cells  <- the honest test\n",
            sum(w$`2 relative + true anchor` < w$`4 anchor only`, na.rm = TRUE), nrow(w)))
cat(sprintf("relative(blind anchor) beats absolute:   %d/%d cells\n",
            sum(w$`3 relative + blind anchor` < w$`1 absolute`, na.rm = TRUE), nrow(w)))

# ── The within-country no-op, verified rather than asserted ──────────────────
cat("\n---- within-country check: relative target is a location shift ----\n")
d <- frames[[paste("Ghana", "child_iron")]]
if (!is.null(d)) {
  X <- prep_X(d); keep <- apply(X, 2, function(z) stats::sd(z) > 0)
  X <- X[, keep, drop = FALSE]
  wt <- area_precision_weights(d$svy_prev, d$n_svy)
  c0 <- lg(nat_of("Ghana", "child_iron"))
  f <- function(y) { sel <- .awsl_screen(X, y, K_SCREEN)
    .awsl_predict(.awsl_stack(X[, sel, drop = FALSE], y, wt), X[, sel, drop = FALSE]) }
  a <- f(lg(d$svy_prev)); b <- f(lg(d$svy_prev) - c0) + c0
  cat(sprintf("  Ghana child_iron, max |absolute - (relative + constant)| = %.3e\n",
              max(abs(a - b), na.rm = TRUE)))
}
cat("\n-> results/tables/relative_target_loco.csv\n")
