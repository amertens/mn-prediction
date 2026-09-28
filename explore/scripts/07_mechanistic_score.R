# =============================================================================
# explore/scripts/07_mechanistic_score.R   [probe MX-01, scoring step]
#
# QUESTION. Do ~10 nutrient-specific mechanistic columns, built from the crop
# basket and the food-composition table, beat 383 nutrient-agnostic ones?
#
# THE CONTROL THAT MAKES THIS READABLE. A mechanistic feature set could win for
# the boring reason that it is a generic agro-ecology proxy - cereal-vs-root
# country is also climate. So every cell is scored TWICE: once with the
# features matched to its own nutrient, and once with a MISMATCHED set (the
# zinc cell gets the vitamin A features, etc.). If mismatched does as well as
# matched, there is no nutrient specificity and the entry must say so.
#
# ARMS
#   mech_matched     the nutrient's own mechanistic features
#   mech_mismatched  another nutrient's features (the specificity control)
#   mech_all         all 16 mechanistic columns
#   index_plus_mech  the record's index with the matched features added
#   + domain_index, spatial, null on identical folds
#
#   Rscript explore/scripts/07_mechanistic_score.R
# -> explore/out/07_mechanistic_cells.csv, 07_mechanistic_loco.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()

MX <- read.csv(file.path(EXP_ROOT, "explore/out/06_mechanistic_features.csv"),
               stringsAsFactors = FALSE)
MXMETA <- read.csv(file.path(EXP_ROOT, "explore/out/06_mechanistic_features_metadata.csv"),
                   stringsAsFactors = FALSE)
MXCOLS <- MXMETA$column
GENERAL <- MXMETA$column[MXMETA$nutrient == "general"]

#' Which mechanistic columns belong to which outcome's nutrient
nutrient_of <- function(outcome) {
  if (grepl("zinc", outcome)) "zinc"
  else if (grepl("iron", outcome)) "iron"
  else if (grepl("vitA", outcome)) "vitA"
  else if (grepl("folate", outcome)) "folate"
  else if (grepl("b12", outcome)) "b12"
  else "general"
}
# B12 is animal-source-food only and has no crop-basket feature: it gets the
# general composition axes, and the entry says the mechanism is not built.
cols_for <- function(nutr) {
  own <- MXMETA$column[MXMETA$nutrient == nutr]
  unique(c(own, GENERAL))
}
MISMATCH <- c(zinc = "vitA", vitA = "zinc", iron = "folate",
              folate = "iron", b12 = "zinc", general = "vitA")

#' A ridge on a named subset of the mechanistic columns
make_mech_arm <- function(pick) {
  function(tr, te, y, X, D, aux) {
    cc <- intersect(pick, colnames(X))
    if (length(cc) < 2) return(rep(mean(y[tr]), length(te)))
    .v2_enet(X[tr, cc, drop = FALSE], y[tr], X[te, cc, drop = FALSE], alpha = 0)
  }
}
#' The index on the domain axes, with the mechanistic block added as its own axis
make_index_plus <- function(pick) {
  function(tr, te, y, X, D, aux) {
    cc <- intersect(pick, colnames(X))
    D2 <- if (length(cc)) cbind(D, X[, cc, drop = FALSE]) else D
    arm_domain_index_v2(tr, te, y, X, D2, aux)
  }
}

rows <- list(); loco <- list()
ix <- exp_cell_index(E)

for (i in seq_len(nrow(ix))) {
  cn <- ix$country[i]; on <- ix$outcome[i]
  nutr <- nutrient_of(on)
  pick_m <- cols_for(nutr)
  pick_x <- cols_for(MISMATCH[[nutr]])
  arms <- c(exp_baseline_arms(),
            list(mech_matched    = make_mech_arm(pick_m),
                 mech_mismatched = make_mech_arm(pick_x),
                 mech_all        = make_mech_arm(MXCOLS),
                 index_plus_mech = make_index_plus(pick_m)))
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, cn, on, tgt, extra = MX,
                              extra_domain = "Crop-basket mechanism"),
                     error = function(e) NULL)
    if (is.null(cell)) next
    r1 <- exp_infill(cell, arms, reps = REPS); r1$nutrient <- nutr
    r2 <- exp_region(cell, arms); if (!is.null(r2)) r2$nutrient <- nutr
    rows[[paste(i, tgt, "A")]] <- r1
    rows[[paste(i, tgt, "B")]] <- r2
  }
  message("  ", cn, " ", on, " (", nutr, ")")
}

# ── transport ───────────────────────────────────────────────────────────────
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    nutr <- nutrient_of(on)
    arms <- c(exp_baseline_arms()[c("null_train_mean", "domain_index")],
              list(mech_matched    = make_mech_arm(cols_for(nutr)),
                   mech_mismatched = make_mech_arm(cols_for(MISMATCH[[nutr]])),
                   mech_all        = make_mech_arm(MXCOLS),
                   index_plus_mech = make_index_plus(cols_for(nutr))))
    cl <- exp_all_cells(E, tgt, outcomes = on, extra = MX,
                        extra_domain = "Crop-basket mechanism")
    if (length(cl) < 3) next
    l <- exp_loco(cl, arms, domain_of = c(E$domain_of,
            stats::setNames(rep("Crop-basket mechanism", length(MXCOLS)), MXCOLS)))
    if (!is.null(l)) l$nutrient <- nutr
    loco[[paste(tgt, on)]] <- l
    message("  LOCO ", tgt, " ", on)
  }
}

raw <- dplyr::bind_rows(rows)
sm <- exp_summarise(raw)
sm$nutrient <- vapply(sm$outcome, nutrient_of, "")
exp_write(sm, "07_mechanistic_cells")
lc <- dplyr::bind_rows(loco); exp_write(lc, "07_mechanistic_loco")

cat("\n== in-fill, level target: median Spearman by arm ==\n")
a <- sm[sm$estimand == "infill" & sm$target == "level", ]
print(aggregate(spearman ~ arm, data = a,
                FUN = function(z) round(median(z, na.rm = TRUE), 3)), row.names = FALSE)

cat("\n== SPECIFICITY: matched vs mismatched, by nutrient (in-fill, level) ==\n")
sp <- a[a$arm %in% c("mech_matched", "mech_mismatched", "domain_index"), ]
w <- reshape(sp[, c("country", "outcome", "nutrient", "arm", "spearman")],
             idvar = c("country", "outcome", "nutrient"), timevar = "arm",
             direction = "wide")
names(w) <- sub("^spearman\\.", "", names(w))
w$matched_minus_mismatched <- round(w$mech_matched - w$mech_mismatched, 3)
print(w[order(-w$matched_minus_mismatched), ], row.names = FALSE)

cat("\n== by nutrient (mean over cells, in-fill level) ==\n")
print(aggregate(cbind(mech_matched, mech_mismatched, domain_index) ~ nutrient,
                data = w, FUN = function(z) round(mean(z, na.rm = TRUE), 3)),
      row.names = FALSE)

cat("\n== transport (LOCO) ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
print(b[order(b$target, -b$spearman), ], row.names = FALSE)
