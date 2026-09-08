# =============================================================================
# scripts/protocol_v2/49_espen_admin2_block.R   [HB-01]
#
# HELMINTH BURDEN AND CONTROL, FROM ESPEN, AT ADMIN-2
#
# Hookworm and schistosomiasis cause chronic blood loss and are the most
# mechanistically direct missing predictor for the iron outcomes; the domain
# was scaffolded in July but never populated because the keyed ESPEN API had
# closed. The portal's own export endpoint needs no key:
#   https://espen.afro.who.int/api/download-data/{ISO2}/{sth|sch}/iu/{from}/{to}
# (found 2026-09-07 in the site's route table; ISO-2 country codes, lowercase
# disease). The eight files (four countries x STH / SCH, implementation-unit
# level, 2014-2025) sit in data/ESPEN/raw/. The IU is ESPEN's ADM2 in every
# country: Ghana 292 IUs (2014 districts), Gambia 44, Sierra Leone 16, Malawi
# 29 (its districts, this project's Admin-1, broadcast to Traditional
# Authorities as the IHME block is).
#
# The export carries no continuous prevalence, only the programme's endemicity
# CLASS per IU-year, plus mass drug administration (MDA) delivery and
# coverage. Features per disease (dz = sth, sch):
#   espen_<dz>_prev_mid       class midpoint (%) at the year nearest the survey
#   espen_<dz>_prev_mid_base  class midpoint at the earliest reported year
#   espen_<dz>_mda_share      share of 2014-2018 with MDA delivered in the IU
#   espen_<dz>_cov_mean       mean reported epidemiological coverage (%) over
#                             the years MDA was delivered, 2014-2018
# Class midpoints: non-endemic 0; surveillance / < 2% -> 1; 1-9% -> 5;
# 2-9% -> 5.5; < 10% -> 5; < 20% -> 10; 10-19% -> 15; 10-49% -> 30;
# 20-49% -> 35; >= 50% -> 65; "Not reported" / "Unknown" -> NA.
# Sierra Leone's survey (2013) predates ESPEN's first year, so its
# survey-year value is the 2014 class.
#
#   Rscript scripts/protocol_v2/49_espen_admin2_block.R
# -> data/covariates/harmonized/predictors_admin2_espen.csv
#    data/covariates/cluster/predictors_cluster_espen.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
HDIR <- "data/covariates/harmonized"; CDIR <- "data/covariates/cluster"
source("R/survey_years.R"); SURVEY_YEAR <- survey_years()   # single source: metadata/survey_years.csv (Gambia 2018, Ghana 2017, Malawi 2016, Sierra Leone 2013)
ISO2 <- c(Gambia = "GM", Ghana = "GH", Malawi = "MW", SierraLeone = "SL")
JOIN_LEVEL <- c(Gambia = "Admin2", Ghana = "Admin2", Malawi = "Admin1", SierraLeone = "Admin2")
S <- read.csv(file.path(HDIR, "predictors_admin2_shared.csv"), check.names = FALSE)
spine <- distinct(S, country, Admin1, Admin2)
kk <- function(x) { x <- iconv(as.character(x), to = "ASCII//TRANSLIT", sub = ""); tolower(gsub("[^A-Za-z0-9]", "", x)) }
match_names <- function(src, tgt, label) {
  ks <- kk(src); kt <- kk(tgt); m <- match(ks, kt); n_exact <- sum(!is.na(m))
  if (any(is.na(m)) && requireNamespace("stringdist", quietly = TRUE)) {
    idx <- which(is.na(m)); dm <- stringdist::stringdistmatrix(ks[idx], kt, method = "jw", p = 0.1)
    best <- apply(dm, 1, which.min); ok <- dm[cbind(seq_along(idx), best)] <= 0.15; m[idx[ok]] <- best[ok] }
  cat(sprintf("    %-24s %3d exact, %3d after fuzzy, %3d unmatched of %3d\n", label, n_exact, sum(!is.na(m)), sum(is.na(m)), length(src)))
  if (any(is.na(m))) cat("      unmatched:", paste(head(src[is.na(m)], 12), collapse = " | "), "\n")
  tgt[m] }
mid_of <- function(x) { x <- tolower(x); dplyr::case_when(
  grepl("non-endemic", x) ~ 0, grepl("surveillance", x) | grepl("less than 2", x) ~ 1,
  grepl("2%-9%", x) ~ 5.5, grepl("1%-9%", x) ~ 5, grepl("less than 10", x) ~ 5, grepl("less than 20", x) ~ 10,
  grepl("10%-19%", x) ~ 15, grepl("20%-49%", x) ~ 35, grepl("10%-49%", x) ~ 30, grepl("50% and above", x) ~ 65,
  TRUE ~ NA_real_) }
num <- function(x) suppressWarnings(as.numeric(x))
ALIAS <- list(
  Ghana = c("Accra Metropolitan Area" = "Accra", "Wa" = "Wa Municipal", "Ho" = "Ho Municipal", "East Akim" = "Abuakwa South",
            "Lower Denkyira" = "Twifo-Hemang-Lower Denkyira", "Axim Municipal" = "Nzema East", "Guan" = "Krachi East"),
  Gambia = c("Basse" = "Fulladu East", "Jimara" = "Fulladu East", "Jumara" = "Fulladu East", "Tumana" = "Fulladu East",
             "Upper Fulladu West" = "Fulladu West", "Lower Fulladu West" = "Fulladu West", "Farafenni Town" = "Upper Baddibu"),
  SierraLeone = c("Karene" = "Bombali", "Falaba" = "Koinadugu"),
  Malawi = c("Likoma Islands" = "Likoma"))

blocks <- list()
for (cn in names(SURVEY_YEAR)) { cat(sprintf("\n%s (survey %d, join at %s)\n", cn, SURVEY_YEAR[[cn]], JOIN_LEVEL[[cn]]))
  per_dz <- list()
  for (dz in c("sth", "sch")) {
    f <- sprintf("data/ESPEN/raw/espen_%s_%s_iu_2014_2025.csv", ISO2[[cn]], dz); if (!file.exists(f)) { cat("  missing", f, "\n"); next }
    d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE); d$year <- num(d$year); d$mid <- mid_of(d$endemicity)
    d$mda <- !is.na(d$mdaScheme) & d$mdaScheme != "" & d$mdaScheme != "Not delivered"
    d$covv <- num(d$epiCov); d$covv[!is.finite(d$covv) | d$covv <= 0] <- num(d$cov)[!is.finite(d$covv) | d$covv <= 0]; d$covv[!is.finite(d$covv) | d$covv <= 0] <- NA
    yr <- SURVEY_YEAR[[cn]]
    g <- d |> group_by(iuCode, iu_name = admin2) |> summarise(
      prev_mid = { ok <- is.finite(mid); if (!any(ok)) NA_real_ else mid[ok][order(abs(year[ok] - yr), year[ok])][1] },
      prev_mid_base = { ok <- is.finite(mid); if (!any(ok)) NA_real_ else mid[ok][order(year[ok])][1] },
      mda_share = mean(mda[year >= 2014 & year <= 2018]),
      cov_mean = { s <- covv[mda & year >= 2014 & year <= 2018 & is.finite(covv)]; if (length(s)) mean(s) else 0 },
      .groups = "drop")
    names(g)[3:6] <- paste0("espen_", dz, "_", names(g)[3:6])
    cat(sprintf("  %s: %d IUs | class at survey year known for %d | MDA in >=1 year 2014-18 for %d\n", toupper(dz), nrow(g), sum(is.finite(g[[3]])), sum(g[[5]] > 0)))
    per_dz[[dz]] <- g }
  if (!length(per_dz)) next
  g <- Reduce(function(a, b) full_join(a, b, by = c("iuCode", "iu_name")), per_dz)
  lvl <- JOIN_LEVEL[[cn]]; tgt <- unique(spine[[lvl]][spine$country == cn])
  # ESPEN carries post-split or programme names the spine does not; map them to the spine unit that contains them
  # (each alias is applied only if its target exists in the spine)
  al <- ALIAS[[cn]]; if (!is.null(al)) { al <- al[al %in% tgt]; hit <- g$iu_name %in% names(al); g$iu_name[hit] <- al[g$iu_name[hit]]; if (any(hit)) cat(sprintf("    %d IU names re-pointed by alias\n", sum(hit))) }
  g$.unit <- match_names(g$iu_name, tgt, paste(cn, lvl))
  g <- g[!is.na(g$.unit), ]
  vals <- setdiff(names(g), c("iuCode", "iu_name", ".unit"))
  w <- g |> group_by(.unit) |> summarise(across(all_of(vals), ~ mean(.x, na.rm = TRUE)), n_iu = n(), .groups = "drop")
  w[vals] <- lapply(w[vals], function(v) { v[!is.finite(v)] <- NA; v })
  if (any(w$n_iu > 1)) cat(sprintf("    %d spine units receive the mean of %d-%d IUs\n", sum(w$n_iu > 1), min(w$n_iu[w$n_iu > 1]), max(w$n_iu[w$n_iu > 1])))
  names(w)[names(w) == ".unit"] <- lvl; w$n_iu <- NULL
  blocks[[cn]] <- spine[spine$country == cn, ] |> left_join(w, by = lvl)
  cov <- blocks[[cn]] |> summarise(across(all_of(vals), ~ mean(is.finite(.x)))); cat("    coverage of spine:", paste(sprintf("%s %.2f", sub("espen_", "", names(cov)), unlist(cov)), collapse = ", "), "\n") }
OUT <- bind_rows(blocks)
write.csv(OUT, file.path(HDIR, "predictors_admin2_espen.csv"), row.names = FALSE)
TC <- read.csv("results/tables/cluster_level/targets_cluster.csv") |> distinct(country, cluster, Admin1, Admin2)
CL <- TC |> left_join(OUT, by = c("country", "Admin1", "Admin2"))
write.csv(CL, file.path(CDIR, "predictors_cluster_espen.csv"), row.names = FALSE)
cat(sprintf("\nwritten: %d Admin-2 rows x %d espen columns; %d cluster rows (%.2f matched)\n", nrow(OUT), ncol(OUT) - 3, nrow(CL), mean(is.finite(CL$espen_sth_prev_mid))))
cat("\n-- district medians by country --\n")
print(as.data.frame(OUT |> group_by(country) |> summarise(across(starts_with("espen_"), ~ round(median(.x, na.rm = TRUE), 1)))), row.names = FALSE)
cat("\nDONE\n")
