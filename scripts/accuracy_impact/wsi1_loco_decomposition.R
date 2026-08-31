# =============================================================================
# scripts/accuracy_impact/wsi1_loco_decomposition.R
#
# Score the anchored-level-plus-covariate-tilt decomposition under
# leave-one-country-out, sweeping the tilt so the ranking/error trade-off is
# reported rather than chosen.
#
#   Rscript scripts/accuracy_impact/wsi1_loco_decomposition.R
# -> results/tables/loco_decomposition.csv
# -> results/tables/loco_decomposition_summary.csv
#
# The blind anchor here is the LEAVE-ONE-COUNTRY-OUT MEAN of the training
# countries, which is the weakest honest anchor available: it is what you would
# use for a country with no national estimate at all. VMNIS supplies a better
# one where it has the country, and the gap between the two columns is the value
# of having a national estimate.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full")
TDIR  <- here("results", "tables")
TILTS <- c(0, 0.2, 0.35, 0.5, 0.75, 1)
SEED  <- 20260926L
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))

H <- suppressMessages(readr::read_csv(
  here("data", "covariates", "harmonized", "predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- setdiff(names(H), c("country", "Admin1", "Admin2"))
cfgs <- get_country_configs()

# ---------------------------------------------------------------------------
# Assemble one district table per outcome, pooled across countries.
# ---------------------------------------------------------------------------
gather_outcome <- function(ocn) {
  rows <- list()
  for (cn in names(cfgs)) {
    cc <- cfgs[[cn]]
    if (!ocn %in% names(cc$outcomes)) next
    sv <- tryCatch(targets::tar_read_raw(
      paste0("svy_admin2_", tolower(cn), "_", ocn), store = STORE),
      error = function(e) NULL)
    if (is.null(sv) || !nrow(sv)) next
    hc <- H[H$country == cn, , drop = FALSE]
    if (!nrow(hc)) next
    m <- dplyr::inner_join(sv, hc, by = admin2_join_by(sv, hc))
    if (nrow(m) < 8) next
    m$country <- cc$country
    rows[[cn]] <- m
  }
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

OUTCOMES <- unique(unlist(lapply(cfgs, function(z) names(z$outcomes))))
res <- list()

for (ocn in OUTCOMES) {
  d <- gather_outcome(ocn)
  if (is.null(d) || length(unique(d$country)) < 3) {
    cat(sprintf("[wsi1] %-13s skipped: needs 3+ countries\n", ocn)); next
  }
  # Covariates usable in EVERY country, so the transported model is not
  # silently different per fold.
  usable <- COVS[vapply(COVS, function(v) {
    x <- d[[v]]
    all(vapply(split(x, d$country), function(z)
      mean(is.finite(z)) > 0.8 && stats::sd(z, na.rm = TRUE) > 0, logical(1)))
  }, logical(1))]
  if (length(usable) < 5) {
    cat(sprintf("[wsi1] %-13s skipped: %d usable covariates\n", ocn, length(usable))); next
  }

  for (ho in unique(d$country)) {
    te <- d[d$country == ho, , drop = FALSE]
    tr <- d[d$country != ho, , drop = FALSE]
    if (nrow(te) < 8 || length(unique(tr$country)) < 2) next

    Xtr <- as.matrix(tr[, usable, drop = FALSE])
    Xte <- as.matrix(te[, usable, drop = FALSE])
    # Median-impute from the TRAINING rows only.
    for (j in seq_len(ncol(Xtr))) {
      mu <- stats::median(Xtr[, j], na.rm = TRUE)
      Xtr[!is.finite(Xtr[, j]), j] <- mu
      Xte[!is.finite(Xte[, j]), j] <- mu
    }
    set.seed(SEED)
    f <- tryCatch(.ds_fit(Xtr, tr$svy_prev, k_screen = min(20L, ncol(Xtr))),
                  error = function(e) NULL)
    if (is.null(f)) next
    pred <- tryCatch(.ds_predict(f, Xte), error = function(e) NULL)
    if (is.null(pred)) next

    # The two anchors.
    anchor_true  <- stats::weighted.mean(te$svy_prev, pmax(te$n_svy, 1), na.rm = TRUE)
    anchor_blind <- stats::weighted.mean(tr$svy_prev, pmax(tr$n_svy, 1), na.rm = TRUE)
    spread <- loco_training_spread(tr$svy_prev, tr$country)
    if (!is.finite(spread)) spread <- stats::sd(tr$svy_prev, na.rm = TRUE)

    sc <- loco_score_decomposition(te$svy_prev, pred, anchor_true, anchor_blind,
                                   spread, tilts = TILTS)
    if (is.null(sc)) next
    sc$outcome <- ocn; sc$held_out <- ho; sc$n_districts <- nrow(te)
    sc$n_cov <- length(usable); sc$spread_assumed <- round(spread, 4)
    sc$anchor_true_pp  <- round(100 * anchor_true, 2)
    sc$anchor_blind_pp <- round(100 * anchor_blind, 2)
    # The raw transported prediction, for reference: this is the `absolute` arm.
    o <- is.finite(te$svy_prev) & is.finite(pred)
    sc$absolute_r <- if (sum(o) > 3 && stats::sd(pred[o]) > 0)
      round(suppressWarnings(stats::cor(pred[o], te$svy_prev[o])), 4) else NA_real_
    sc$absolute_mae_pp <- round(100 * mean(abs(pred[o] - te$svy_prev[o])), 2)
    res[[paste(ocn, ho)]] <- sc
    cat(sprintf("[wsi1] %-13s hold out %-13s districts=%3d cov=%3d anchor_true=%5.1f blind=%5.1f\n",
                ocn, ho, nrow(te), length(usable),
                100 * anchor_true, 100 * anchor_blind))
  }
}

if (!length(res)) stop("[wsi1] nothing scored")
out <- dplyr::bind_rows(res)
front <- c("outcome", "held_out", "anchor", "tilt")
out <- out[, c(front, setdiff(names(out), front))]
readr::write_csv(out, file.path(TDIR, "loco_decomposition.csv"))

cat("\n=== the tilt sweep, pooled over cells ===\n")
s <- out |> group_by(anchor, tilt) |>
  summarise(cells = dplyr::n(),
            r = round(mean(r, na.rm = TRUE), 3),
            mae_pp = round(mean(mae_pp, na.rm = TRUE), 2),
            abs_bias_pp = round(mean(abs(bias_pp), na.rm = TRUE), 2),
            .groups = "drop") |> arrange(anchor, tilt)
print(as.data.frame(s), row.names = FALSE)

cat("\n=== against the un-anchored transported model (the `absolute` arm) ===\n")
a <- out |> distinct(outcome, held_out, absolute_r, absolute_mae_pp)
cat(sprintf("absolute:      mean r %.3f, MAE %.2f pp, %d cells\n",
            mean(a$absolute_r, na.rm = TRUE), mean(a$absolute_mae_pp, na.rm = TRUE), nrow(a)))
b0 <- out[out$anchor == "true"  & out$tilt == 0, ]
bb <- out[out$anchor == "blind" & out$tilt == 0, ]
cat(sprintf("anchor only (true):   MAE %.2f pp\n", mean(b0$mae_pp, na.rm = TRUE)))
cat(sprintf("anchor only (blind):  MAE %.2f pp\n", mean(bb$mae_pp, na.rm = TRUE)))

cat("\n=== best tilt per held-out country, true anchor ===\n")
bt <- out |> filter(anchor == "true") |> group_by(held_out, tilt) |>
  summarise(r = round(mean(r, na.rm = TRUE), 3),
            mae = round(mean(mae_pp, na.rm = TRUE), 2), .groups = "drop")
print(as.data.frame(bt |> group_by(held_out) |> slice_min(mae, n = 1)), row.names = FALSE)
readr::write_csv(s, file.path(TDIR, "loco_decomposition_summary.csv"))
cat(sprintf("\n-> %s\n-> %s\n", file.path("results","tables","loco_decomposition.csv"),
            file.path("results","tables","loco_decomposition_summary.csv")))
