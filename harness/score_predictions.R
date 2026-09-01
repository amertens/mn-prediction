#!/usr/bin/env Rscript
# =============================================================================
# harness/score_predictions.R
#
# WS8c. Score a partner's prediction file against the locked held-out cells.
#
#   Rscript harness/score_predictions.R <predictions.csv> [out.csv]
#
# The prediction file must carry: country, outcome, Admin1, Admin2, pred
# (a proportion in [0, 1]).
#
# WHAT IT REFUSES, AND WHY
# ------------------------
# A file containing predictions for TRAINING cells is rejected outright rather
# than filtered. Section 2.3 measures a mean +0.182 gain in rank correlation
# from letting a covariate prescreen see the held-out regions, positive in all
# 20 measurable cells. A submission that has seen the test cells cannot be
# scored against one that has not, and quietly dropping the extra rows would
# hide that the submitter had them.
#
# WHAT IT REPORTS
#   r            Pearson correlation against the design-based district estimate
#   spearman     rank correlation
#   mae_pp       mean absolute error
#   bias_pp      SIGNED mean error, reported separately from mae_pp because the
#                two answer different questions and Section 6 conflated them
#   r_share      r divided by the EMPIRICAL reliability ceiling, suppressed
#                where that ceiling is below 0.05
#   topk_capture share of the truly worst k districts that the prediction also
#                ranks in its worst k
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: score_predictions.R <predictions.csv> [out.csv]")
HDIR <- dirname(normalizePath(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = FALSE))
if (!dir.exists(HDIR)) HDIR <- "harness"

pred <- utils::read.csv(args[1], stringsAsFactors = FALSE)
need <- c("country", "outcome", "Admin1", "Admin2", "pred")
miss <- setdiff(need, names(pred))
if (length(miss)) stop("prediction file is missing: ", paste(miss, collapse = ", "))

truth <- utils::read.csv(file.path(HDIR, "targets_heldout.csv"), stringsAsFactors = FALSE)
pub   <- utils::read.csv(file.path(HDIR, "targets_public.csv"), stringsAsFactors = FALSE)
ceil  <- tryCatch(utils::read.csv(file.path(HDIR, "ceilings.csv"), stringsAsFactors = FALSE),
                  error = function(e) NULL)

k <- function(d) paste(d$country, d$outcome, d$Admin1, d$Admin2)
kp <- k(pred); kt <- k(truth); ktr <- k(pub)

# REFUSAL 1: predictions for training cells.
leak <- intersect(kp, ktr)
if (length(leak)) {
  cat(sprintf("REFUSED: the file contains %d prediction(s) for TRAINING cells.\n",
              length(leak)))
  cat("The held-out set is listed in harness/heldout_cells.csv. Submit predictions\n")
  cat("for those cells only. Examples of the offending rows:\n")
  cat(paste0("  ", utils::head(leak, 5), collapse = "\n"), "\n")
  quit(status = 2)
}
# REFUSAL 2: nothing to score.
if (!length(intersect(kp, kt))) {
  cat("REFUSED: no prediction matches a held-out cell.\n"); quit(status = 2)
}
if (any(!is.finite(pred$pred)) || any(pred$pred < 0 | pred$pred > 1, na.rm = TRUE)) {
  cat("REFUSED: `pred` must be a finite proportion in [0, 1].\n"); quit(status = 2)
}

m <- dplyr::inner_join(truth, pred[, need], by = c("country","outcome","Admin1","Admin2"))
m <- m[is.finite(m$svy_prev) & is.finite(m$pred), , drop = FALSE]

topk <- function(obs, prd, frac = 0.25) {
  n <- length(obs); kk <- max(1L, round(n * frac))
  length(intersect(order(obs, decreasing = TRUE)[seq_len(kk)],
                   order(prd, decreasing = TRUE)[seq_len(kk)])) / kk
}
out <- m |> group_by(country, outcome) |>
  summarise(n_districts = dplyr::n(),
            r = round(suppressWarnings(stats::cor(svy_prev, pred)), 4),
            spearman = round(suppressWarnings(
              stats::cor(svy_prev, pred, method = "spearman")), 4),
            mae_pp = round(100 * mean(abs(pred - svy_prev)), 3),
            bias_pp = round(100 * mean(pred - svy_prev), 3),
            topk_capture = round(topk(svy_prev, pred), 3),
            .groups = "drop")
if (!is.null(ceil)) {
  kk2 <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
  out$r_max_emp <- ceil$r_max_emp[match(paste(kk2(out$country), out$outcome),
                                        paste(kk2(ceil$country), ceil$outcome))]
  out$r_share <- ifelse(is.finite(out$r_max_emp) & out$r_max_emp > 0.05,
                        round(out$r / out$r_max_emp, 3), NA_real_)
}
cat("=== scored on the locked held-out cells ===\n")
print(as.data.frame(out), row.names = FALSE)
cat(sprintf("\ncells scored: %d | median r: %.3f | mean MAE: %.2f pp | mean signed bias: %+.2f pp\n",
            nrow(out), stats::median(out$r, na.rm = TRUE),
            mean(out$mae_pp, na.rm = TRUE), mean(out$bias_pp, na.rm = TRUE)))
if ("r_share" %in% names(out))
  cat(sprintf("median r_share against the empirical ceiling: %.3f (%d of %d cells scorable)\n",
              stats::median(out$r_share, na.rm = TRUE),
              sum(is.finite(out$r_share)), nrow(out)))
if (length(args) > 1) { utils::write.csv(out, args[2], row.names = FALSE)
  cat(sprintf("-> %s\n", args[2])) }
