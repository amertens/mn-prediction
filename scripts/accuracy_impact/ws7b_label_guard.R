# =============================================================================
# scripts/accuracy_impact/ws7b_label_guard.R
#
# WS7b. Derive the assay guard from Stata variable LABELS and compare it with
# the name-based guard, column by column.
#
# WHY THIS IS THE CORRECT FIX AND THE OTHER ONE IS A WORKAROUND
# -------------------------------------------------------------
# Section 7.5 states the residual risk plainly: the merged data carries no
# variable labels, so is_biomarker_column() is a regex over opaque names with no
# ground truth, and its failure mode is silent. It recommends obtaining the
# original Stata files and re-deriving the guard from labels.
#
# WS7a demonstrated the failure mode is not hypothetical. gw_wm_whbc is women's
# haemoglobin in g/dL and passed the guard for as long as the guard has existed,
# because the token is spelled in lower case. A label reading "haemoglobin
# (g/dL)" would have caught it on the first pass.
#
# WHAT IS AVAILABLE, AND WHAT IS NOT
# ----------------------------------
# Stage 0 found labelled .dta files for THREE of the four surveys:
#   Gambia         data/IPD/Gambia/*.dta
#   Ghana          data/IPD/Ghana/*.dta
#   Sierra Leone   data/IPD/Sierra Leone/*.dta
# Malawi has no .dta. Its source is an RDS whose 242 questionnaire columns are
# coded m01, m115a, m220h and carry ZERO labels, and the only local
# documentation is a 349-character README pointing at immpact@cdc.gov. Malawi is
# therefore reported as not coverable rather than silently omitted.
#
#   Rscript scripts/accuracy_impact/ws7b_label_guard.R
# -> results/tables/assay_guard_label_comparison.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(here)})
Sys.setenv(COVARIATE_VOCAB = "harmonized")
targets::tar_source(here("R"))
TDIR <- here("results", "tables")

DTA <- list(
  Gambia = c("data/IPD/Gambia/Gambia_Child data_Davis-Berkeley.dta",
             "data/IPD/Gambia/Gambia_Woman data_Davis-Berkeley.dta"),
  Ghana  = c("data/IPD/Ghana/GMS_ChildrenAllData_DavisBerkeley.dta",
             "data/IPD/Ghana/GMS_WomenAllData_DavisBerkeley.dta"),
  `Sierra Leone` = c("data/IPD/Sierra Leone/Sierra Leone_Child data_Davis-Berkeley.dta",
                     "data/IPD/Sierra Leone/Sierra Leone_Woman data_Davis-Berkeley.dta")
)

# A guard derived from what a label SAYS, not from what a name looks like.
# Deliberately generous on the blood-draw side: this predicate exists to find
# columns the name-based guard misses, so a false positive here is a finding to
# adjudicate and a false negative is the failure being hunted.
LABEL_DRAW <- paste(
  "ferritin", "transferrin", "\\bstfr\\b", "soluble transferrin",
  "retinol", "\\brbp\\b", "retinol.?binding",
  "c.?reactive", "\\bcrp\\b", "alpha.?1.?acid", "\\bagp\\b",
  "serum zinc", "plasma zinc", "serum folate", "red cell folate",
  "vitamin b.?12", "cobalamin", "\\bmma\\b", "methylmalonic",
  "blood sample", "venipuncture", "venepuncture", "phlebotom",
  "specimen", "assay", "\\bserum\\b", "\\bplasma\\b",
  "iron deficien", "vitamin a deficien", "zinc deficien", "folate deficien",
  "inflammation.?adjust", "brinda", "thurnham",
  sep = "|")
LABEL_HB <- paste("h(a)?emoglobin", "\\bhb\\b", "an(a)?emi", sep = "|")
# Labels that MENTION a blood word but describe knowledge, belief or an
# exposure. Section 7.4 verified these individually on the name side.
LABEL_KEEP <- paste("heard", "believe", "think", "opinion", "knowledge",
                    "fortif", "supplement", "tablet", "syrup", "consent",
                    "willing", "refus", sep = "|")

classify_label <- function(lab) {
  l <- tolower(trimws(as.character(lab)))
  keep <- grepl(LABEL_KEEP, l)
  draw <- grepl(LABEL_DRAW, l) & !keep
  hb   <- grepl(LABEL_HB, l) & !keep & !draw
  out <- rep("questionnaire", length(l))
  out[hb] <- "hb_field"; out[draw] <- "blood_draw"
  out[!nzchar(l)] <- "unlabelled"
  out
}

rows <- list()
for (cn in names(DTA)) {
  for (p in DTA[[cn]]) {
    fp <- here(p)
    if (!file.exists(fp)) { cat("[missing]", p, "\n"); next }
    d <- tryCatch(haven::read_dta(fp, n_max = 1L), error = function(e) NULL)
    if (is.null(d)) { cat("[unreadable]", p, "\n"); next }
    labs <- vapply(names(d), function(v) {
      a <- attr(d[[v]], "label"); if (is.null(a)) "" else as.character(a)[1]
    }, character(1))
    rows[[length(rows) + 1L]] <- data.frame(
      country = cn, file = basename(p), raw_name = names(d),
      label = unname(labs), stringsAsFactors = FALSE)
    cat(sprintf("  %-13s %-46s %4d cols, %4d labelled\n", cn, basename(p),
                ncol(d), sum(nzchar(labs))))
  }
}
if (!length(rows)) stop("No .dta files readable.")
L <- dplyr::bind_rows(rows)
L <- L[!duplicated(paste(L$country, L$raw_name)), ]

# The pipeline prefixes survey columns with gw_. Compare on that name.
L$pipeline_name <- paste0("gw_", L$raw_name)
L$label_class <- classify_label(L$label)
L$name_class  <- biomarker_column_class(L$pipeline_name)
L$name_blocked <- is_biomarker_column(L$pipeline_name)

# Restrict to columns the pipeline actually carries, so the comparison is about
# predictors that can reach a model rather than about the raw file.
cfgs <- get_country_configs(); pipe <- character(0)
for (cn2 in names(cfgs)) for (ocn in names(cfgs[[cn2]]$outcomes)) {
  od <- tryCatch(targets::tar_read_raw(
    paste0("outcome_data_", tolower(cn2), "_", ocn), store = here("_targets_full")),
    error = function(e) NULL)
  if (!is.null(od)) pipe <- c(pipe, od$Xvars_full)
}
pipe <- unique(pipe)
L$in_pipeline <- L$pipeline_name %in% pipe

L$agreement <- dplyr::case_when(
  L$label_class == "unlabelled" ~ "no label",
  L$label_class == L$name_class ~ "agree",
  L$label_class != "questionnaire" & L$name_class == "questionnaire" ~
    "LABEL BLOCKS, NAME ALLOWS",
  L$label_class == "questionnaire" & L$name_class != "questionnaire" ~
    "NAME BLOCKS, LABEL ALLOWS",
  TRUE ~ "different class")

readr::write_csv(L[, c("country","file","raw_name","pipeline_name","label",
                       "label_class","name_class","name_blocked","in_pipeline",
                       "agreement")],
                 file.path(TDIR, "assay_guard_label_comparison.csv"))

cat("\n=== coverage ===\n")
print(as.data.frame(L |> group_by(country) |>
  summarise(columns = dplyr::n(), labelled = sum(label_class != "unlabelled"),
            in_pipeline = sum(in_pipeline), .groups = "drop")), row.names = FALSE)
cat("\nMalawi: NOT COVERABLE. No .dta exists; its 242 questionnaire columns are\n",
    "coded m01/m115a/m220h with zero labels. Recorded as an open request.\n")

P <- L[L$in_pipeline & L$label_class != "unlabelled", ]
cat(sprintf("\n=== agreement on %d labelled columns the pipeline carries ===\n", nrow(P)))
print(as.data.frame(table(P$agreement)), row.names = FALSE)

crit <- P[P$agreement == "LABEL BLOCKS, NAME ALLOWS", ]
cat(sprintf("\n=== THE FINDING THAT MATTERS: %d column(s) a label blocks and the name does not ===\n",
            nrow(crit)))
if (nrow(crit)) print(as.data.frame(crit[, c("country","pipeline_name","label","label_class")]),
                      row.names = FALSE)
over <- P[P$agreement == "NAME BLOCKS, LABEL ALLOWS", ]
cat(sprintf("\n=== possible over-blocking: %d column(s) the name blocks and the label does not ===\n",
            nrow(over)))
if (nrow(over)) print(as.data.frame(utils::head(over[, c("country","pipeline_name","label")], 40)),
                      row.names = FALSE)
