# =============================================================================
# scripts/accuracy_impact/wsb3_ship.R
#
# WS-B3. Produce the deliverable objects from the tournament winner.
#
# The estimator is the empirical Bayes blend (R/eb_district_estimator.R). Every
# cell passes through the WS-D calibration gate and the WS1a reliability
# suppression rule before it reaches a deliverable.
#
#   Rscript scripts/accuracy_impact/wsb3_ship.R
# -> results/deliverables/district_estimates.rds  (+ .csv)
# -> results/deliverables/national_regional.csv
# -> results/tables/shipped_estimates_summary.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
STORE <- here("_targets_full"); TDIR <- here("results","tables")
DDIR <- here("results","deliverables")
dir.create(DDIR, showWarnings = FALSE, recursive = TRUE)
SEED <- 20260924L

REL <- read.csv(file.path(TDIR,"reliability_empirical.csv"), stringsAsFactors=FALSE)
REL <- REL[REL$scheme=="within",]
GATE <- tryCatch(read.csv(file.path(TDIR,"calibration_gate_report.csv"),
                          stringsAsFactors=FALSE), error=function(e) NULL)
THR <- read.csv(here("config","who_thresholds.csv"), stringsAsFactors=FALSE)
kk <- function(x) tolower(gsub("[^a-z]","",tolower(x)))
lam_of <- function(cn,ocn){ i <- which(kk(REL$country)==kk(cn) & REL$outcome==ocn)
  if (length(i)) REL$lambda_emp[i[1]] else NA_real_ }
rmax_of <- function(cn,ocn){ i <- which(kk(REL$country)==kk(cn) & REL$outcome==ocn)
  if (length(i)) REL$r_max_emp[i[1]] else NA_real_ }

# FITNESS_FOR_USE section 3: suppress the Admin-2 layer below this reliability.
SUPPRESS_RMAX <- 0.30

cfgs <- get_country_configs(); out <- list(); nat_rows <- list()
for (cn in names(cfgs)) {
  cc <- cfgs[[cn]]
  for (ocn in names(cc$outcomes)) {
    oc <- cc$outcomes[[ocn]]
    sv <- tryCatch(targets::tar_read_raw(paste0("svy_admin2_",tolower(cn),"_",ocn),
                                         store=STORE), error=function(e) NULL)
    od <- tryCatch(targets::tar_read_raw(paste0("outcome_data_",tolower(cn),"_",ocn),
                                         store=STORE), error=function(e) NULL)
    if (is.null(sv) || !nrow(sv)) next
    lam <- lam_of(cc$country, ocn); rmx <- rmax_of(cc$country, ocn)
    eb <- tryCatch(eb_district_estimate(sv, lambda_emp = lam), error=function(e) NULL)
    if (is.null(eb)) next
    eb <- eb_rank_intervals(eb, B = 2000L, seed = SEED)
    eb$country <- cc$country; eb$outcome <- ocn
    eb$r_max_emp <- rmx
    # Suppression: below the reliability cut the Admin-2 layer is withheld and
    # the reason travels with the row rather than the row disappearing.
    eb$admin2_released <- is.finite(rmx) && rmx >= SUPPRESS_RMAX
    eb$suppression_reason <- ifelse(eb$admin2_released, NA_character_,
      sprintf("empirical reliability %.3f below the %.2f cut", rmx, SUPPRESS_RMAX))
    # Calibration gate verdict for this cell, if the report covers it.
    g <- if (!is.null(GATE)) GATE[kk(GATE$country)==kk(cc$country) &
                                  GATE$outcome==ocn, , drop=FALSE] else NULL
    eb$calibration_status <- if (!is.null(g) && nrow(g)) g$status[1] else "not_assessed"
    out[[paste(cn,ocn)]] <- eb

    # National and regional, from the survey directly. FITNESS section 1 says
    # these are the reportable quantities; the district layer is a ranking.
    if (!is.null(od)) {
      nat <- tryCatch(national_design_based(od, cc, oc), error=function(e) NULL)
      a1  <- tryCatch(admin1_design_based(od, cc, oc), error=function(e) NULL)
      th <- THR[THR$outcome == ocn, , drop=FALSE]
      band <- function(p) if (!nrow(th) || !is.finite(p)) NA_character_ else
        if (p < th$band_none_upper[1]) "low" else
        if (p < th$band_mild_upper[1]) "mild" else
        if (p < th$band_moderate_upper[1]) "moderate" else "severe"
      if (!is.null(nat))
        nat_rows[[length(nat_rows)+1L]] <- data.frame(
          country=cc$country, outcome=ocn, level="national", unit="(national)",
          prev=round(nat$prev,6), n=nat$n,
          se=round(sqrt(nat$prev*(1-nat$prev)/max(nat$n,1)),6),
          who_band=band(nat$prev), band_source_verified=if(nrow(th)) th$source_verified[1] else NA,
          stringsAsFactors=FALSE)
      if (!is.null(a1)) nat_rows[[length(nat_rows)+1L]] <- data.frame(
          country=cc$country, outcome=ocn, level="admin1", unit=a1$Admin1,
          prev=round(a1$prev,6), n=a1$n,
          se=round(sqrt(a1$prev*(1-a1$prev)/pmax(a1$n,1)),6),
          who_band=vapply(a1$prev, band, character(1)),
          band_source_verified=if(nrow(th)) th$source_verified[1] else NA,
          stringsAsFactors=FALSE)
    }
  }
}
D <- dplyr::bind_rows(out)
front <- c("country","outcome","Admin1","Admin2")
D <- D[, c(front, setdiff(names(D), front))]
saveRDS(D, file.path(DDIR,"district_estimates.rds"))
readr::write_csv(D, file.path(DDIR,"district_estimates.csv"))
NR <- dplyr::bind_rows(nat_rows)
readr::write_csv(NR, file.path(DDIR,"national_regional.csv"))

summ <- D |> group_by(country, outcome) |> summarise(
  districts = dplyr::n(),
  median_lambda = round(stats::median(lambda, na.rm=TRUE),3),
  median_shift_pp = round(stats::median(abs(100*(eb_prev - svy_prev)), na.rm=TRUE),2),
  median_rank_width = round(stats::median(rank_hi - rank_lo, na.rm=TRUE),1),
  r_max_emp = round(r_max_emp[1],3),
  admin2_released = admin2_released[1],
  calibration = calibration_status[1],
  tau2_source = tau2_source[1], .groups="drop")
readr::write_csv(summ, file.path(TDIR,"shipped_estimates_summary.csv"))

cat("=== WS-B3: shipped district estimates ===\n")
print(as.data.frame(summ), row.names=FALSE)
cat(sprintf("\ncells: %d | Admin-2 layer released: %d | suppressed: %d\n",
            nrow(summ), sum(summ$admin2_released), sum(!summ$admin2_released)))
cat(sprintf("cells failing the calibration gate: %d\n",
            sum(summ$calibration == "calibration_failed")))
cat(sprintf("tau2 from the split-half reliability in %d of %d cells\n",
            sum(summ$tau2_source == "split_half_reliability"), nrow(summ)))
cat(sprintf("\nnational and regional rows: %d\n", nrow(NR)))
