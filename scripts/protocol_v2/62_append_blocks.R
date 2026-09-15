# =============================================================================
# scripts/protocol_v2/62_append_blocks.R   [HC-01 / RT-01 / NV-01, 2026-09-15]
#
# APPEND THE STAND-ALONE PREDICTOR BLOCKS TO THE SHARED ADMIN-2 SET
#
# Each block is a pair of files on the shared spine (554 units, NA where a
# country has no data) written by its own builder:
#   hces         data/covariates/harmonized/predictors_admin2_hces.csv        scripts/covariates/build_hces_diet_block.R
#   rtfp         data/covariates/harmonized/predictors_admin2_rtfp.csv        scripts/covariates/build_rtfp_price_block.R
#   ndvi_modis   data/covariates/harmonized/predictors_admin2_ndvi_modis.csv  scripts/covariates/extract_gee_ndvi_modis.py
#   mics         data/covariates/harmonized/predictors_admin2_mics.csv        scripts/covariates/build_mics_admin2_block.R
# plus <block>_metadata.csv with column, domain, source, subnational,
# assumption, n_countries, countries, completeness. A block whose files are
# absent is skipped with a message, so the step is safe on a partial rebuild.
# Nothing is filtered here: the audit and the exclusion policy (script 53) see
# every column with its completeness.
#
# Audit-only columns stay in the block files and are NOT appended:
#   hces_n_hh, hces_level (sample size, own-households vs parent-mean flag)
#   rtfp_n_markets        (markets inside the polygon)
#   mics_n_hh, mics_level (households behind the estimate, own-unit vs parent)
#
# Rebuild order (docs/findings/SANDBOX_LOG_2026-09.md): build_shared_predictor_set
# -> 07 -> 08 -> 59 -> [60] -> 62 -> 53 -> stamp_predictor_metadata -> audit.
# Re-running is idempotent: earlier columns of every listed block are replaced,
# and columns of a block that no longer exists are removed.
#
#   Rscript -e "source('scripts/protocol_v2/62_append_blocks.R')"
# -> data/covariates/harmonized/predictors_admin2_shared.csv (+ _metadata.csv)
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"
BLOCKS <- c(hces = "^hces_", rtfp = "^rtfp_", ndvi_modis = "^ndvi_modis_", mics = "^mics_(?!heat_)")   # block -> column prefix it owns (perl regex; mics_heat_ belongs to script 59)
AUDIT_ONLY <- c("hces_n_hh", "hces_level", "rtfp_n_markets", "mics_n_hh", "mics_level")
META_COLS <- c("column", "domain", "source", "n_countries", "countries", "completeness", "subnational")

SH  <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE, stringsAsFactors = FALSE)
SHM <- read.csv(file.path(HDIR, "predictors_admin2_shared_metadata.csv"), stringsAsFactors = FALSE)
SHM$subnational <- as.logical(toupper(as.character(SHM$subnational)))
stopifnot(ncol(SH) - 3L == nrow(SHM), setequal(setdiff(names(SH), c("country", "Admin1", "Admin2")), SHM$column))
for (f in c("predictors_admin2_shared.csv", "predictors_admin2_shared_metadata.csv"))
  if (!file.exists(file.path(HDIR, paste0(f, ".pre_blocks"))))
    file.copy(file.path(HDIR, f), file.path(HDIR, paste0(f, ".pre_blocks")))

key <- function(d) paste(d$country, d$Admin1, d$Admin2)
blocks <- list(); meta <- list()
for (b in names(BLOCKS)) {
  fd <- file.path(HDIR, sprintf("predictors_admin2_%s.csv", b)); fm <- file.path(HDIR, sprintf("predictors_admin2_%s_metadata.csv", b))
  if (!file.exists(fd) || !file.exists(fm)) { cat(sprintf("%s: block files absent, skipped\n", b)); next }
  d <- read.csv(fd, check.names = FALSE, stringsAsFactors = FALSE)
  m <- read.csv(fm, stringsAsFactors = FALSE)
  keep <- setdiff(intersect(m$column, names(d)), AUDIT_ONLY)
  stopifnot(length(keep) > 0, all(grepl(BLOCKS[[b]], keep, perl = TRUE)), nrow(d) == nrow(SH), setequal(key(d), key(SH)), !anyDuplicated(key(d)))
  d <- d[match(key(SH), key(d)), ]
  blocks[[b]] <- d[, c("country", "Admin1", "Admin2", keep)]
  meta[[b]] <- m |> filter(column %in% keep) |> transmute(column, domain, source, n_countries, countries, completeness, subnational = as.logical(toupper(as.character(subnational))))
  cat(sprintf("%s: %d columns (%d audit-only left in the block file)\n", b, length(keep), length(intersect(names(d), AUDIT_ONLY))))
}
stopifnot(length(blocks) > 0)
AB <- Reduce(function(x, y) left_join(x, y, by = c("country", "Admin1", "Admin2")), blocks)
MD <- bind_rows(meta)[, META_COLS]
newcols <- MD$column
old <- names(SH)[Reduce(`|`, lapply(BLOCKS, function(rx) grepl(rx, names(SH), perl = TRUE)))]   # every column any listed block owns
if (length(old)) cat("replacing earlier block columns:", length(old), "\n")
SH2  <- SH |> select(-any_of(c(newcols, old))) |> left_join(AB, by = c("country", "Admin1", "Admin2"))
SHM2 <- bind_rows(SHM |> filter(!column %in% c(newcols, old)), MD)
stopifnot(nrow(SH2) == nrow(SH), ncol(SH2) - 3L == nrow(SHM2), !anyDuplicated(names(SH2)), setequal(SHM2$column, setdiff(names(SH2), c("country", "Admin1", "Admin2"))))
write.csv(SH2, file.path(HDIR, "predictors_admin2_shared.csv"), row.names = FALSE)
write.csv(SHM2, file.path(HDIR, "predictors_admin2_shared_metadata.csv"), row.names = FALSE)
cat(sprintf("\nshared set: %d units x %d predictors (was %d)\n", nrow(SH2), ncol(SH2) - 3L, ncol(SH) - 3L))
print(as.data.frame(MD[, c("column", "domain", "n_countries", "completeness")]), row.names = FALSE)
cat("\nNext: 53_apply_exclusion_policy.R, then scripts/covariates/stamp_predictor_metadata.R, then audit_predictor_set.R\nDONE\n")
