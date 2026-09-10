# =============================================================================
# scripts/policy_deck/04_civ_climate_soil_prediction.R
#
# COTE D'IVOIRE ON THE PROTOCOL-V2 CLIMATE + SOIL INDEX
#
# Cote d'Ivoire has a full proxy database and no biomarker survey. This script
# fits the two-domain (Climate and weather + Soil characteristics) index on the
# four surveyed countries and ranks CIV's 33 districts with it - the exact use
# case the transport results describe.
#
# The harness is copied from scripts/protocol_v2/02b_merge_and_loco.R estimand C
# so the numbers are comparable: rank-normalise within country, pool on the
# COMMON columns, build domain PCs with the orientation learned from training
# countries only, then the zero-tuning index.
#
# Two guards, because a new country silently joining a vocabulary is exactly how
# this pipeline has broken before:
#   GUARD 1  reproduce the published four-country transport without CIV
#   GUARD 2  repeat it with CIV in the column intersection, and report the cost
#
# Inputs : results/transportability/civ_canonical_admin2_full.rds (script pair
#          run earlier: extract_gee_admin2 -> cov_harmonize_country -> fill)
# Outputs: results/tables/policy_deck/civ_climate_soil_ranking.csv
#          results/tables/policy_deck/civ_transport_guards.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(sf)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
set.seed(20260909L)

P2   <- "results/tables/protocol_v2"
OUTT <- "results/tables/policy_deck"; dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)

TG  <- read.csv(file.path(P2, "targets_v2.csv"), stringsAsFactors = FALSE)
S   <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD  <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
CIV <- readRDS("results/transportability/civ_canonical_admin2_full.rds")

CS_DOMAINS <- c("Climate and weather", "Soil characteristics")
PREDS <- intersect(MD$column[MD$domain %in% CS_DOMAINS], names(S))
domain_of <- stats::setNames(MD$domain, MD$column)
COUNTRIES <- c("Gambia", "Ghana", "Malawi", "SierraLeone")
cat(sprintf("climate + soil vocabulary: %d columns; CIV holds %d\n",
            length(PREDS), length(intersect(PREDS, names(CIV)))))

# ── cells ────────────────────────────────────────────────────────────────────
build_cell <- function(cn, on, target = "level") {
  ycol <- if (target == "prev") "y_prev" else "y_level"
  wcol <- if (target == "prev") "n_eff"  else "n_eff_cont"
  t <- TG[TG$country == cn & TG$outcome == on, ]
  t <- t[is.finite(t[[ycol]]) & is.finite(t[[wcol]]), ]
  if (nrow(t) < 12) return(NULL)
  m <- inner_join(t, S[S$country == cn, c("Admin1", "Admin2", PREDS)],
                  by = c("Admin1", "Admin2"))
  if (nrow(m) < 12 || dplyr::n_distinct(m$Admin1) < 3) return(NULL)
  Xr <- prep_predictors_v2(as.matrix(m[, PREDS, drop = FALSE]))
  if (ncol(Xr) < 10) return(NULL)
  list(country = cn, n = nrow(m), y_nat = m[[ycol]],
       y_mod = if (target == "prev") .v2_logit(m[[ycol]]) else m[[ycol]],
       X = Xr, w = m[[wcol]], Admin1 = m$Admin1, Admin2 = m$Admin2)
}

civ_cell <- function() {
  m <- CIV[, c("Admin1", "Admin2", intersect(PREDS, names(CIV)))]
  Xr <- prep_predictors_v2(as.matrix(m[, intersect(PREDS, names(CIV)), drop = FALSE]))
  list(country = "CoteDIvoire", n = nrow(m), y_nat = rep(NA_real_, nrow(m)),
       y_mod = rep(NA_real_, nrow(m)), X = Xr, w = rep(1, nrow(m)),
       Admin1 = m$Admin1, Admin2 = m$Admin2)
}

# ── one LOCO pass; `extra` is an unlabelled country carried in the pool ───────
loco_pass <- function(on, extra = NULL, target = "level") {
  cl <- list()
  for (cn in COUNTRIES) {
    z <- tryCatch(build_cell(cn, on, target), error = function(e) NULL)
    if (!is.null(z)) cl[[cn]] <- z
  }
  if (length(cl) < 3) return(NULL)
  pool <- cl
  if (!is.null(extra)) pool[[extra$country]] <- extra
  common <- Reduce(intersect, lapply(pool, function(z) colnames(z$X)))
  if (length(common) < 10) return(NULL)

  Xm   <- do.call(rbind, lapply(pool, function(z) z$X[, common, drop = FALSE]))
  ctry <- rep(names(pool), vapply(pool, function(z) z$n, 0L))
  Y    <- unlist(lapply(pool, function(z) {
            if (all(is.na(z$y_mod))) rep(NA_real_, z$n) else as.numeric(scale(z$y_mod)) }))
  ynat <- unlist(lapply(pool, function(z) z$y_nat))
  wv   <- unlist(lapply(pool, function(z) z$w))
  labelled <- ctry %in% names(cl)

  out <- list(n_common = length(common), scores = NULL, civ = NULL)

  # held-out training countries: the reproduction guard
  sc <- list()
  for (cn in names(cl)) {
    te <- which(ctry == cn); tr <- which(labelled & ctry != cn)
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    p  <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, Dm, NULL),
                   error = function(e) rep(NA_real_, length(te)))
    r  <- suppressWarnings(cor(ynat[te], p, method = "spearman", use = "complete.obs"))
    sc[[cn]] <- data.frame(outcome = on, heldout = cn, n_units = length(te),
                           spearman = r, n_common = length(common))
  }
  out$scores <- bind_rows(sc)

  # the unlabelled country: trained on every labelled row
  if (!is.null(extra)) {
    te <- which(ctry == extra$country); tr <- which(labelled)
    Dm <- domain_representation_v2(Xm, domain_of, sign_rows = tr)
    p  <- tryCatch(arm_domain_index_v2(tr, te, Y, Xm, Dm, NULL),
                   error = function(e) rep(NA_real_, length(te)))
    out$civ <- data.frame(outcome = on, Admin1 = extra$Admin1, Admin2 = extra$Admin2,
                          index = p)
  }
  out
}

OUTCOMES <- c("child_iron", "child_vitA", "women_iron", "women_vitA",
              "women_folate", "women_b12")
CE <- civ_cell()
cat(sprintf("CIV cell: %d districts, %d usable columns after coverage screen\n",
            CE$n, ncol(CE$X)))

g1 <- list(); g2 <- list(); civ_pred <- list()
for (on in OUTCOMES) {
  a <- loco_pass(on, extra = NULL)          # GUARD 1: four countries only
  b <- loco_pass(on, extra = CE)            # GUARD 2: CIV in the intersection
  if (!is.null(a)) g1[[on]] <- transform(a$scores, arm = "4 countries only")
  if (!is.null(b)) {
    g2[[on]] <- transform(b$scores, arm = "with CIV in the vocabulary")
    if (!is.null(b$civ)) civ_pred[[on]] <- b$civ
  }
}
G <- bind_rows(bind_rows(g1), bind_rows(g2))
write.csv(G, file.path(OUTT, "civ_transport_guards.csv"), row.names = FALSE)

summ <- G |> group_by(arm) |>
  summarise(cells = n(), mean_spearman = mean(spearman, na.rm = TRUE),
            positive = sum(spearman > 0, na.rm = TRUE),
            n_common = median(n_common), .groups = "drop")
cat("\n=== GUARDS: leave-one-country-out on the climate + soil index (level) ===\n")
print(as.data.frame(summ), row.names = FALSE)
cat("\npublished reference (nested_domain_selection.csv, arm fixed_cs, level):",
    "mean 0.369 over 22 cells, 22 positive\n")

CP <- bind_rows(civ_pred)
CP <- CP |> group_by(outcome) |>
  mutate(rank = rank(-index, ties.method = "average"),
         pct  = 100 * (rank - 0.5) / n()) |> ungroup()
write.csv(CP, file.path(OUTT, "civ_climate_soil_ranking.csv"), row.names = FALSE)
cat(sprintf("\nCIV rankings written for %d outcomes x %d districts\n",
            dplyr::n_distinct(CP$outcome), dplyr::n_distinct(CP$Admin2)))

ci <- CP[CP$outcome == "child_iron", ]
cat("\nCote d'Ivoire, child iron - ten worst-ranked districts:\n")
print(head(ci[order(ci$rank), c("Admin1", "Admin2", "rank")], 10), row.names = FALSE)
