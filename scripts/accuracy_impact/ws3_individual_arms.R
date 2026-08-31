# =============================================================================
# scripts/accuracy_impact/ws3_individual_arms.R
#
# WS3. Section 5 under one protocol, with the arms nested and the guard fixed.
#
# WHAT CHANGES FROM scripts/covariates/18_individual_anchor.R
# -----------------------------------------------------------
# 1. THE GUARD. WS7a found gw_wm_whbc (women's haemoglobin, g/dL) and gw_gchb
#    (child haemoglobin) passing is_biomarker_column(), so the published
#    questionnaire arm had a blood measurement in it. Both are now blocked.
#
# 2. THE ARMS ARE NESTED. The published questionnaire arm applied the guard to
#    Xvars_full, which also removed the MAP sickle-cell rasters and the DHS
#    mean-haemoglobin Admin-2 aggregates that the proxy arm keeps, because the
#    proxy arm uses Xvars unfiltered. The two arms therefore differed in more
#    than the questionnaire. allowed_under_arm() scopes the filter to the
#    concurrent survey so proxy is a strict subset of questionnaire.
#
# 3. A THIRD ARM (WS3e). "questionnaire + field haemoglobin" permits the
#    hb_field class: haemoglobin, anaemia status and anaemia category, measured
#    with a finger-prick and a portable photometer rather than a venous draw.
#    This is the deployable scenario for most countries, because every DHS
#    already collects it, and the published two-arm design could not answer what
#    it buys.
#
# 4. TWO PROTOCOLS (WS3a). Section 3's headline individual-level number and
#    Section 5's differ by FOLD CONSTRUCTION, which Stage 0 established by
#    reading the code: aggregate_admin2_sl() aggregates res$yhat_full, which
#    src/analysis/sl_helpers.R produces with origami::make_folds(cluster_ids),
#    a CLUSTER-blocked K-fold. Section 5 uses region-blocked leave-one-region-
#    out. Both are run here on the same cells so the gap is measured rather than
#    inferred.
#
#   PROFILE=smoke     Ghana child_iron, both protocols, three arms
#   WS3_CELLS=parity  the 4-cell protocol-parity subset only
#
#   Rscript scripts/accuracy_impact/ws3_individual_arms.R
# -> results/tables/individual_arms_2026-09.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))

PROFILE <- Sys.getenv("PROFILE", "full")
STORE <- here("_targets_full")
SUF <- switch(PROFILE, smoke = "_SMOKE", "")
K_SCREEN <- 40L; SEED <- 20260906L; MIN_N <- 5L; KFOLD <- 5L
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

# The four shared outcomes, matching Section 5's scope.
OUTCOMES <- c("child_iron", "child_vitA", "women_iron", "women_vitA")
# Protocol parity subset: spans all four countries and both nutrients.
PARITY <- list(c("Ghana","child_iron"), c("Gambia","women_iron"),
               c("Malawi","child_vitA"), c("SierraLeone","women_iron"))

fit_pred <- function(Xtr, ytr, Xte) {
  sel <- .awsl_screen(Xtr, ytr, K_SCREEN)
  s <- tryCatch(.awsl_stack(Xtr[, sel, drop = FALSE], ytr, rep(1, length(ytr))),
                error = function(e) NULL)
  if (is.null(s)) return(rep(NA_real_, nrow(Xte)))
  .awsl_predict(s, Xte[, sel, drop = FALSE])
}
prep <- function(d, vars) {
  vars <- intersect(vars, names(d))
  vars <- vars[vapply(vars, function(v) is.numeric(d[[v]]) || is.logical(d[[v]]) ||
                        inherits(d[[v]], "haven_labelled"), logical(1))]
  if (!length(vars)) return(matrix(numeric(0), nrow = nrow(d), ncol = 0))
  X <- vapply(vars, function(v) num(d[[v]]), numeric(nrow(d)))
  if (is.null(dim(X))) X <- matrix(X, nrow = nrow(d))
  colnames(X) <- vars
  for (j in seq_len(ncol(X))) {
    m <- stats::median(X[, j], na.rm = TRUE)
    X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0
  }
  X[, apply(X, 2, function(z) stats::sd(z) > 0), drop = FALSE]
}

cfgs <- get_country_configs()
sel_cells <- NULL
if (PROFILE == "smoke") sel_cells <- list(c("Ghana","child_iron"))
if (Sys.getenv("WS3_CELLS") == "parity") sel_cells <- PARITY
kk <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))

rows <- list()
set.seed(SEED)
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; lc <- tolower(cn)
  for (on in OUTCOMES) {
    if (!is.null(sel_cells) &&
        !any(vapply(sel_cells, function(z) kk(z[1]) == kk(cn) && z[2] == on, logical(1)))) next
    oc <- cc$outcomes[[on]]; if (is.null(oc)) next
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE),
                   error = function(e) NULL)
    if (is.null(od)) next
    d <- od$data
    if (!all(c("Admin2", oc$binary) %in% names(d))) next
    y <- num(d[[oc$binary]]); keep <- is.finite(y)
    d <- d[keep, , drop = FALSE]; y <- y[keep]
    if (length(unique(y)) < 2 || nrow(d) < 100) next
    clu <- if (cc$cluster_id %in% names(d)) as.character(d[[cc$cluster_id]]) else NA
    blk <- if ("Admin1" %in% names(d)) as.character(d$Admin1) else rep("a", nrow(d))
    if (dplyr::n_distinct(blk) < 3) next

    full <- od$Xvars_full %||% od$Xvars
    ARMS <- list(
      proxy            = od$Xvars,
      questionnaire    = full[allowed_under_arm(full, "questionnaire")],
      questionnaire_hb = full[allowed_under_arm(full, "questionnaire_hb")]
    )

    for (arm in names(ARMS)) {
      X <- prep(d, ARMS[[arm]]); if (ncol(X) < 5) next
      for (proto in c("region_loro", "cluster_kfold")) {
        # region_loro   leave-one-region-out, the strict protocol Section 5 uses
        # cluster_kfold cluster-blocked K-fold, which is what
        #               src/analysis/sl_helpers.R produces and therefore what
        #               Section 3's individual-level number rests on
        fold <- if (proto == "region_loro") blk else {
          cu <- unique(clu[!is.na(clu)])
          asg <- sample(rep_len(seq_len(KFOLD), length(cu))); names(asg) <- cu
          as.character(asg[clu])
        }
        if (dplyr::n_distinct(fold[!is.na(fold)]) < 3) next
        oof <- rep(NA_real_, nrow(X))
        for (f in unique(fold[!is.na(fold)])) {
          i <- which(fold == f); tr <- setdiff(seq_len(nrow(X)), i)
          if (!length(tr) || length(i) >= nrow(X) || length(unique(y[tr])) < 2) next
          oof[i] <- fit_pred(X[tr, , drop = FALSE], y[tr], X[i, , drop = FALSE])
        }
        ok <- is.finite(oof)
        if (sum(ok) < 50) next
        for (unit in c("district", "cluster")) {
          key <- if (unit == "district") paste(blk, d$Admin2) else clu
          if (all(is.na(key))) next
          agg <- data.frame(key = key[ok], obs = y[ok], pred = oof[ok]) |>
            group_by(key) |>
            summarise(n = dplyr::n(), obs = mean(obs), pred = mean(pred),
                      .groups = "drop") |>
            filter(n >= MIN_N)
          if (nrow(agg) < 8) next
          rows[[length(rows) + 1L]] <- data.frame(
            country = cc$country, outcome = on, arm = arm, protocol = proto,
            unit = unit, n_units = nrow(agg), n_pred = ncol(X),
            r = round(suppressWarnings(stats::cor(agg$obs, agg$pred)), 4),
            mae_pp = round(100 * mean(abs(agg$obs - agg$pred)), 2),
            stringsAsFactors = FALSE)
          cat(sprintf("  %-13s %-11s %-17s %-14s %-8s p=%4d units=%3d r=%+.3f\n",
                      cn, on, arm, proto, unit, ncol(X), nrow(agg),
                      utils::tail(rows, 1)[[1]]$r))
        }
      }
    }
  }
}
res <- dplyr::bind_rows(rows)
if (!nrow(res)) stop("No rows produced.")
readr::write_csv(res, here("results", "tables",
                           sprintf("individual_arms_2026-09%s.csv", SUF)))

D <- res[res$unit == "district", ]
cat("\n=== WS3a: protocol by arm, district ===\n")
print(as.data.frame(D |> group_by(protocol, arm) |>
  summarise(cells = dplyr::n(), mean_r = round(mean(r, na.rm = TRUE), 3),
            med_r = round(stats::median(r, na.rm = TRUE), 3),
            mae = round(mean(mae_pp, na.rm = TRUE), 2),
            mean_p = round(mean(n_pred)), .groups = "drop")), row.names = FALSE)

cat("\n=== WS3e: arm gains over proxy, strict protocol, district ===\n")
S <- D[D$protocol == "region_loro", ]
w <- S |> select(country, outcome, arm, r) |>
  tidyr::pivot_wider(names_from = arm, values_from = r)
if (all(c("proxy","questionnaire") %in% names(w))) {
  w$gain_quest <- round(w$questionnaire - w$proxy, 3)
  if ("questionnaire_hb" %in% names(w))
    w$gain_hb <- round(w$questionnaire_hb - w$proxy, 3)
  print(as.data.frame(w), row.names = FALSE)
  cat("\n--- aggregates ---\n")
  agg1 <- function(lab, z) cat(sprintf(
    "%-34s cells %2d | mean gain quest %+0.3f (better %d) | mean gain +hb %+0.3f (better %d)\n",
    lab, nrow(z), mean(z$gain_quest, na.rm = TRUE), sum(z$gain_quest > 0, na.rm = TRUE),
    if ("gain_hb" %in% names(z)) mean(z$gain_hb, na.rm = TRUE) else NA_real_,
    if ("gain_hb" %in% names(z)) sum(z$gain_hb > 0, na.rm = TRUE) else NA_integer_))
  agg1("all cells", w)
  agg1("WS3b: excluding Malawi", w[kk(w$country) != "malawi", ])
  agg1("WS3c: excluding Ghana women_iron",
       w[!(kk(w$country) == "ghana" & w$outcome == "women_iron"), ])
  agg1("WS3b+3c: excluding both",
       w[kk(w$country) != "malawi" &
         !(kk(w$country) == "ghana" & w$outcome == "women_iron"), ])
}

cat("\n=== WS3f: cluster against district, same protocol and arm ===\n")
u <- res[res$protocol == "region_loro", ] |> select(country, outcome, arm, unit, r) |>
  tidyr::pivot_wider(names_from = unit, values_from = r)
if (all(c("district","cluster") %in% names(u))) {
  u$gain <- round(u$cluster - u$district, 3)
  print(as.data.frame(u |> group_by(country) |>
    summarise(cells = dplyr::n(), mean_gain = round(mean(gain, na.rm = TRUE), 3),
              better = sprintf("%d/%d", sum(gain > 0, na.rm = TRUE), dplyr::n()),
              .groups = "drop")), row.names = FALSE)
}
