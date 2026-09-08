# =============================================================================
# scripts/protocol_v2/45_lsms_flunet_block.R   [AD-lsms]
#
# LSMS AND FLUNET, BACK IN AS AN ADD-ON BLOCK
#
# The January 2026 deck counted 275 LSMS and 3 FluNet predictors. The
# harmonisation notes (scripts/covariates/harmonize_extra_domains.R) record
# why neither reached the vocabulary: LSMS exists for Ghana only (GLSS7,
# data/LSMS/Ghana_LSMS_clean.RDS, a svyby object of region-level means and
# standard errors, 10 old regions x 276 columns), and FluNet is national
# weekly surveillance with no within-country variation and, for our four
# countries, data for Ghana and Sierra Leone only. Neither is a
# harmonisation failure; both are limits of the sources.
#
# This script builds what CAN be built and hands it to the add-on harness:
#   lsms_*     GLSS7 region means broadcast to Ghana's districts through the
#              16-to-10 region crosswalk (standard errors and ids dropped)
#   flunet_*   survey-year national influenza indicators (specimens
#              processed per week, share positive, influenza A share) for
#              Ghana (2017) and Sierra Leone (2013); NA elsewhere
# Under in-fill the LSMS block can only help Ghana; FluNet is constant within
# a country and can help nothing. That is the honest test.
#
#   Rscript scripts/protocol_v2/45_lsms_flunet_block.R
# -> data/covariates/harmonized/predictors_admin2_lsms_flunet.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)[, c("country", "Admin1", "Admin2")]
# ── LSMS (Ghana, GLSS7) ──────────────────────────────────────────────────────
L <- readRDS("data/LSMS/Ghana_LSMS_clean.RDS"); L <- as.data.frame(L, stringsAsFactors = FALSE)
keep <- names(L)[!grepl("^lsms_hh_se[.]|_se$|^lsms_hh_nh$|^lsms_hh_pid$|region", names(L))]
keep <- keep[vapply(keep, function(k) is.numeric(L[[k]]) && sum(is.finite(L[[k]])) >= 8 && stats::sd(L[[k]], na.rm = TRUE) > 0, TRUE)]
reg10 <- as.character(if ("lsms_admin1" %in% names(L)) L$lsms_admin1 else L$lsms_region)
cat("LSMS: regions", length(reg10), "| numeric mean columns kept", length(keep), "of", ncol(L), "\n"); cat("  regions:", paste(reg10, collapse = ", "), "\n")
XW <- tryCatch(read.csv("data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv", stringsAsFactors = FALSE), error = function(e) NULL)
gh <- S[S$country == "Ghana", ]
norm <- function(x) tolower(gsub("[^a-z]", "", tolower(x)))
if (!is.null(XW)) { cat("  crosswalk columns:", paste(names(XW), collapse = ", "), "\n")
  c16 <- XW[[grep("16|new", names(XW), ignore.case = TRUE)[1]]]; c10 <- XW[[grep("10|old", names(XW), ignore.case = TRUE)[1]]]
  gh$reg10 <- c10[match(norm(gh$Admin1), norm(c16))] } else gh$reg10 <- NA_character_
gh$reg10[is.na(gh$reg10)] <- gh$Admin1[is.na(gh$reg10)]          # regions unchanged between the two systems
j <- match(norm(gh$reg10), norm(reg10)); cat("  Ghana districts mapped to an LSMS region:", sum(!is.na(j)), "of", nrow(gh), "\n")
if (any(is.na(j))) cat("  unmapped Admin1 names:", paste(unique(gh$Admin1[is.na(j)]), collapse = ", "), "\n")
LS <- gh[, c("country", "Admin1", "Admin2")]; for (k in keep) LS[[k]] <- L[[k]][j]
# ── FluNet (national, survey year) ──────────────────────────────────────────
F <- read.csv("data/FluNet/VIW_FNT.csv", stringsAsFactors = FALSE)
source("R/survey_years.R"); yr <- survey_years(); iso <- c(Gambia = "GMB", Ghana = "GHA", Malawi = "MWI", SierraLeone = "SLE")
num <- function(x) suppressWarnings(as.numeric(x))
FL <- do.call(rbind, lapply(names(yr), function(cn) { f <- F[F$COUNTRY_CODE == iso[[cn]] & F$ISO_YEAR == yr[[cn]], ]
  if (!nrow(f)) return(data.frame(country = cn, flunet_specimens_per_week = NA_real_, flunet_share_positive = NA_real_, flunet_influenza_a_share = NA_real_))
  sp <- sum(num(f$SPEC_PROCESSED_NB), na.rm = TRUE); pa <- sum(num(f$INF_A), na.rm = TRUE); pall <- sum(num(f$INF_ALL), na.rm = TRUE)
  data.frame(country = cn, flunet_specimens_per_week = sp / nrow(f), flunet_share_positive = if (sp > 0) pall / sp else NA_real_, flunet_influenza_a_share = if (pall > 0) pa / pall else NA_real_) }))
print(FL, row.names = FALSE)
OUT <- left_join(S, LS, by = c("country", "Admin1", "Admin2")) |> left_join(FL, by = "country")
write.csv(OUT, file.path(HDIR, "predictors_admin2_lsms_flunet.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %d rows, %d LSMS columns (Ghana only), 3 FluNet columns (Ghana, Sierra Leone)\nDONE\n", nrow(OUT), length(keep)))
