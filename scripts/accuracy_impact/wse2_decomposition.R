# =============================================================================
# scripts/accuracy_impact/wse2_decomposition.R
#
# WS-E1 and WS-E2.
#
# E1. Stage 0 established by reading ws4b_skill_curve.R that the D8 exclusion
#     was SYMMETRIC: the reduced 197-covariate set (all 97 dhs columns removed)
#     is used at both the DHS-target block and the micronutrient block. No rerun
#     is needed and none is performed; the check is recorded.
#
# E2. The stated interpretation of the skill curve is that these covariates carry
#     between-region socioeconomic signal and almost nothing within region. That
#     is an inference from the resolution sweep. Here it becomes a measured
#     column: each target's out-of-fold correlation is decomposed by residualising
#     both prediction and target on Admin-1.
#
#     r_between  correlation of the region means of prediction and target
#     r_within   correlation of the two after removing region means
#
#   Rscript scripts/accuracy_impact/wse2_decomposition.R
# -> results/tables/skill_curve_decomposition.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results","tables")
SEED <- 20260923L; K_SCREEN <- 20L; MIN_D <- 12L
set.seed(SEED)

H <- suppressMessages(readr::read_csv(
  here("data","covariates","harmonized","predictors_admin2_harmonized.csv"),
  show_col_types = FALSE)) |> as.data.frame()
COVS <- grep("^dhs", setdiff(names(H), c("country","Admin1","Admin2")),
             value = TRUE, invert = TRUE)
cat(sprintf("[E1] covariate set is the reduced one used for BOTH families: %d columns\n",
            length(COVS)))

decomp <- function(obs, pred, reg) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]; pred <- pred[ok]; reg <- as.factor(reg[ok])
  if (length(obs) < MIN_D || nlevels(reg) < 2) return(NULL)
  mo <- tapply(obs, reg, mean); mp <- tapply(pred, reg, mean)
  rb <- if (length(mo) >= 3 && stats::sd(mo) > 0 && stats::sd(mp) > 0)
    suppressWarnings(stats::cor(mo, mp)) else NA_real_
  ro <- obs - mo[reg]; rp <- pred - mp[reg]
  rw <- if (stats::sd(ro) > 0 && stats::sd(rp) > 0)
    suppressWarnings(stats::cor(ro, rp)) else NA_real_
  data.frame(r_overall = suppressWarnings(stats::cor(obs, pred)),
             r_between = rb, r_within = rw,
             var_share_between = round(stats::var(mo[reg]) / stats::var(obs), 4),
             n_units = length(obs), n_regions = nlevels(reg),
             stringsAsFactors = FALSE)
}
fit_oof <- function(m, ycol) {
  X <- as.matrix(m[, COVS, drop = FALSE]); X[!is.finite(X)] <- NA
  keep <- colMeans(is.finite(X)) == 1 & apply(X, 2, function(z) stats::sd(z, na.rm=TRUE) > 0)
  X <- X[, keep, drop = FALSE]; if (ncol(X) < 5) return(NULL)
  oof <- rep(NA_real_, nrow(m))
  for (rg in unique(m$Admin1)) {
    i <- which(m$Admin1 == rg); tr <- setdiff(seq_len(nrow(m)), i)
    if (length(tr) < 6) next
    f <- .ds_fit(X[tr,,drop=FALSE], m[[ycol]][tr], k_screen = K_SCREEN)
    p <- .ds_predict(f, X[i,,drop=FALSE]); if (!is.null(p)) oof[i] <- p
  }
  oof
}

rows <- list()
# ── DHS indicator family ────────────────────────────────────────────────────
FILES <- list(Gambia="Gambia_2019", Ghana="Ghana_2014", Malawi="Malawi_2015",
              `Sierra Leone`="Sierra Leone_2013")
for (cn in names(FILES)) {
  p <- here("data","DHS","clean", paste0(FILES[[cn]], "_dhs_admin2_direct.rds"))
  if (!file.exists(p)) next
  dd <- readRDS(p); hc <- H[H$country == gsub(" ","",cn) | H$country == cn, , drop=FALSE]
  if (!nrow(hc)) next
  for (ind in sort(unique(dd$indicator))) {
    z <- dd[dd$indicator == ind, , drop=FALSE]
    z <- z[is.finite(z$direct.est) & is.finite(z$direct.var) & z$direct.var > 1e-8, , drop=FALSE]
    if (nrow(z) < MIN_D) next
    k <- data.frame(Admin1 = trimws(as.character(z$admin1.name)),
                    Admin2 = trimws(as.character(z$admin2.name)),
                    y = z$direct.est, stringsAsFactors = FALSE)
    m <- dplyr::inner_join(k, hc, by = admin2_join_by(k, hc))
    if (nrow(m) < MIN_D || dplyr::n_distinct(m$Admin1) < 3) next
    oof <- fit_oof(m, "y"); if (is.null(oof)) next
    r <- decomp(m$y, oof, m$Admin1); if (is.null(r)) next
    r$country <- cn; r$target <- ind; r$family <- "DHS indicator"
    rows[[length(rows)+1L]] <- r
  }
}
# ── micronutrient family ────────────────────────────────────────────────────
cfgs <- get_country_configs()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]; hc <- H[H$country == cn, , drop=FALSE]; if (!nrow(hc)) next
  for (ocn in names(cc$outcomes)) {
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_", tolower(cn), "_", ocn),
                                         store=STORE), error=function(e) NULL)
    if (is.null(sv) || nrow(sv) < MIN_D) next
    m <- dplyr::inner_join(sv[, intersect(c("Admin1","Admin2","svy_prev"), names(sv)), drop=FALSE],
                           hc, by = admin2_join_by(sv, hc))
    m <- m[is.finite(m$svy_prev), , drop=FALSE]
    if (nrow(m) < MIN_D || dplyr::n_distinct(m$Admin1) < 3) next
    oof <- fit_oof(m, "svy_prev"); if (is.null(oof)) next
    r <- decomp(m$svy_prev, oof, m$Admin1); if (is.null(r)) next
    r$country <- cc$country; r$target <- ocn; r$family <- "micronutrient"
    rows[[length(rows)+1L]] <- r
  }
}
res <- dplyr::bind_rows(rows)
front <- c("country","target","family")
res <- res[, c(front, setdiff(names(res), front))]
readr::write_csv(res, file.path(TDIR, "skill_curve_decomposition.csv"))

cat("\n=== WS-E2: between- and within-region decomposition ===\n")
print(as.data.frame(res |> group_by(family) |> summarise(
  targets = dplyr::n(),
  med_r_overall = round(stats::median(r_overall, na.rm=TRUE),3),
  med_r_between = round(stats::median(r_between, na.rm=TRUE),3),
  med_r_within  = round(stats::median(r_within,  na.rm=TRUE),3),
  pct_within_positive = round(100*mean(r_within > 0, na.rm=TRUE),1),
  med_var_share_between = round(stats::median(var_share_between, na.rm=TRUE),3),
  .groups="drop")), row.names = FALSE)
cat(sprintf("\ntargets where r_between exceeds r_within: %d of %d\n",
            sum(res$r_between > res$r_within, na.rm=TRUE),
            sum(is.finite(res$r_between) & is.finite(res$r_within))))
