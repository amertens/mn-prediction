# =============================================================================
# scripts/protocol_v2/70_il01_column_audit.R   [IL-02d]
#
# WHICH SURVEY COLUMNS IL-01 ACTUALLY USED
#
# Reproduces the survey-column selection of 46_individual_level_models.R (the
# LEAK regex, 70 percent coverage, the 80-column coverage cap) cell by cell,
# lists the columns that entered, and ranks them by univariate AUC against the
# outcome (orientation-free: max(AUC, 1 - AUC)). Also reports, per country,
# which age / sex / supplementation / genotype columns the regex removes.
# Found 2026-09-27: cluster numbers, child id, region and cluster-level
# sampling-design columns in the Ghana and Gambia sets; haemoglobin (m228
# child, m432 women) in the Malawi set; women's age removed in Ghana.
#
#   Rscript scripts/protocol_v2/70_il01_column_audit.R
# -> results/tables/protocol_v2/il02_il01_survey_columns.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(targets)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
invisible(capture.output(targets::tar_source("R")))
STORE <- "_targets_full"; MAXCOL <- 80
cfg <- get_country_configs()
LEAK <- "RBP|rbp|VAD|vad|LogFer|logfer|Ferr|ferr|IDA|ida$|Brinda|BRINDA|Thurn|Folate|folate|Fol|B12|b12|Zinc|zinc|zn_|Hb|hb$|Hgb|hgb|Haem|Hem|anaem|anem|Anem|CRP|crp|AGP|agp|infl|Infl|UIC|uic|Iod|iod|Salt|salt|MUAC|muac|weight|Weight|Wt$|wt$|_id$|ID$|Id$|cnum|clust|Clust|Date|date|month|Month|year|Year|psu|PSU|strat|Strat|Admin|admin|lat|Lat|lon|Lon|GPS|gps|Team|team|line|Line|hhid|caseid|Retinol|retinol|Vit|vit|Anemia|Malaria|malaria|RDT|rdt|Plasmod|Sickle|G6PD|Transferrin|sTfR|stfr|ZPP|zpp|Ret$|Def$|def$|Adj|adj"
uauc <- function(y, x) { ok <- is.finite(x) & is.finite(y); y <- y[ok]; x <- x[ok]; n1 <- sum(y == 1); n0 <- sum(y == 0); if (n1 < 5 || n0 < 5 || sd(x) == 0) return(NA_real_); r <- rank(x); a <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0); max(a, 1 - a) }
allcols <- list()
for (cn in names(cfg)) { cc <- cfg[[cn]]; lc <- tolower(cn)
  outs <- intersect(names(cc$outcomes), c("child_vitA", "women_vitA", "child_iron", "women_iron", "women_folate", "women_b12"))
  for (on in outs) { oc <- cc$outcomes[[on]]
    if (cn == "Malawi") {
      MW <- readRDS("data/IPD/Malawi/clean_malawi_mn_data.RDS")
      defcol <- c(child_vitA = "vitA_def", women_vitA = "vitA_def", child_iron = "iron_def", women_iron = "iron_def")[on]
      if (is.na(defcol) || !defcol %in% names(MW)) next
      grp <- if (grepl("^child", on)) is.finite(MW$psc_agecat) else is.finite(MW$women_agecat)
      d <- MW[grp, ]; y <- .v2_num(d[[defcol]])
      gw <- names(d)[grepl("^m[0-9]+[a-g]?$|^mvisits$|^mlang|^mtype$|^hhs_|^oil_vita$|^sugar_vita$|^salt$|^fast$", names(d))]
    } else {
      od <- tryCatch(tar_read_raw(paste0("outcome_data_", lc, "_", on), store = STORE), error = function(e) NULL); if (is.null(od)) next
      d <- od$data; yb <- tryCatch(resolve_uniform_outcome(d, cc, oc), error = function(e) NULL)
      y <- if (!is.null(yb)) .v2_num(yb) else .v2_num(d[[oc$binary]])
      outcome_vars <- unique(unlist(lapply(cc$outcomes, function(o) c(o$binary, o$continuous))))
      LEAK2 <- paste0(LEAK, "|^sf_|ferritin|_nmol|^fol|^zn|^map2_|^rbp|^vit|_def$")
      gw <- names(d)[grepl("^gw_", names(d)) & !grepl(LEAK2, names(d)) & !names(d) %in% outcome_vars]
    }
    gw <- gw[vapply(gw, function(k) { v <- .v2_num(d[[k]]); mean(is.finite(v)) >= 0.7 && stats::sd(v, na.rm = TRUE) > 0 && length(unique(v[is.finite(v)])) >= 2 }, TRUE)]
    cov <- vapply(gw, function(k) mean(is.finite(.v2_num(d[[k]]))), 0); gw80 <- gw[order(-cov)][seq_len(min(MAXCOL, length(gw)))]
    ua <- vapply(gw80, function(k) uauc(y, .v2_num(d[[k]])), 0)
    lab <- vapply(gw80, function(k) { l <- attr(d[[k]], "label"); if (is.null(l)) "" else substr(as.character(l)[1], 1, 60) }, "")
    allcols[[paste(cn, on)]] <- data.frame(country = cn, outcome = on, col = gw80, coverage = round(cov[gw80], 2), uni_auc = round(ua, 3), label = lab)
    cat(sprintf("\n== %s %s: eligible %d, kept %d (cap %d). Dropped-by-cap examples: %s\n", cn, on, length(gw), length(gw80), MAXCOL, paste(head(setdiff(gw, gw80), 8), collapse = ", ")))
    top <- allcols[[paste(cn, on)]] |> arrange(desc(uni_auc)) |> head(8)
    print(top, row.names = FALSE)
  }
  # which age / sex / VAS columns exist and did they survive the regex?
  if (cn != "Malawi") { od <- tar_read_raw(paste0("outcome_data_", lc, "_", outs[1]), store = STORE); nm <- names(od$data)
    cand <- nm[grepl("^gw_", nm) & grepl("age|Age|sex|Sex|VAS|VitASupp|Supp|Preg|Lact|genotype|Thal|Sickle|sickle", nm)]
    cat(sprintf("\n-- %s age/sex/supplement/genotype columns and whether IL-01's regex removes them:\n", cn))
    print(data.frame(col = cand, removed_by_regex = grepl(paste0(LEAK, "|^sf_|ferritin|_nmol|^fol|^zn|^map2_|^rbp|^vit|_def$"), cand)), row.names = FALSE)
  }
}
A <- bind_rows(allcols); write.csv(A, "results/tables/protocol_v2/il02_il01_survey_columns.csv", row.names = FALSE)
cat("\nDONE\n")
