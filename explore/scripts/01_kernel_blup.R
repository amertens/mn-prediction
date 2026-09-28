# =============================================================================
# explore/scripts/01_kernel_blup.R   [probe KB-01]
#
# QUESTION. The project has established that capacity is a liability at
# n = 14-87: every tuned learner loses to a zero-tuning domain index. Does the
# quantitative-genetics answer to the same regime - parameterise by the
# similarity between units, and estimate the shrinkage by REML rather than
# selecting it by CV - beat the index?
#
# ARMS
#   blup_all        one linear kernel on every predictor
#   blup_cs         one kernel on climate + soil (the pre-registered domains)
#   blup_5k         five broad blocks (climate, soil, embedding, agriculture,
#                   space) as separate kernels, multi-kernel REML
#   blup_spatial    spatial kernel alone (geography with no covariates)
#   blup_cs_spatial climate+soil and space as separate kernels
#   + the record's comparators on identical folds
#
# Also writes the REML variance components per domain: a model-based
# decomposition of where district-level variation lives, which is what the
# ablation work has been approaching by deletion.
#
#   Rscript explore/scripts/01_kernel_blup.R
#   EXP_REPS=3 for a quick pass
# -> explore/out/01_kernel_blup_cells.csv, 01_kernel_blup_loco.csv,
#    01_kernel_varcomp.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
DOM <- E$domain_of

CLIM <- "Climate and weather"
SOIL <- "Soil characteristics"
EMB  <- "Satellite embedding"
AGRI <- "Agricultural production, land use"

cols_in <- function(X, doms) colnames(X)[which(DOM[colnames(X)] %in% doms)]

ARMS <- c(
  exp_baseline_arms(),
  list(
    blup_all = make_blup_arm(function(X, D, aux) list(all = k_linear(X))),
    blup_cs  = make_blup_arm(function(X, D, aux) {
      cc <- cols_in(X, c(CLIM, SOIL))
      if (length(cc) < 5) return(list())
      list(cs = k_linear(X[, cc, drop = FALSE]))
    }),
    blup_spatial = make_blup_arm(function(X, D, aux)
      list(space = k_spatial(aux$lon, aux$lat))),
    blup_cs_spatial = make_blup_arm(function(X, D, aux) {
      cc <- cols_in(X, c(CLIM, SOIL))
      k <- list(space = k_spatial(aux$lon, aux$lat))
      if (length(cc) >= 5) k$cs <- k_linear(X[, cc, drop = FALSE])
      k
    }),
    # Five broad, mechanistically distinct blocks rather than all 21 domain
    # kernels. Two reasons, one scientific and one practical: 21 domain kernels
    # are far too collinear for the split between them to be attributable, and
    # the 21-kernel arm cost ~170s per cell against ~1s for the rest while
    # never being a contender (0.29-0.37 in the first pass).
    blup_5k = make_blup_arm(function(X, D, aux) {
      k <- list(clim  = k_linear(X[, cols_in(X, CLIM), drop = FALSE]),
                soil  = k_linear(X[, cols_in(X, SOIL), drop = FALSE]),
                emb   = k_linear(X[, cols_in(X, EMB),  drop = FALSE]),
                agri  = k_linear(X[, cols_in(X, AGRI), drop = FALSE]),
                space = k_spatial(aux$lon, aux$lat))
      k[!vapply(k, is.null, TRUE)]
    })
  ))

# ── estimands A and B ───────────────────────────────────────────────────────
rows <- list(); vc <- list()
ix <- exp_cell_index(E)
for (i in seq_len(nrow(ix))) {
  cn <- ix$country[i]; on <- ix$outcome[i]
  for (tgt in c("prev", "level")) {
    cell <- tryCatch(exp_cell(E, cn, on, tgt), error = function(e) NULL)
    if (is.null(cell)) next
    t0 <- Sys.time()
    rows[[paste(i, tgt, "A")]] <- exp_infill(cell, ARMS, reps = REPS)
    rows[[paste(i, tgt, "B")]] <- exp_region(cell, ARMS)

    # Variance components (reporting, not prediction). A REDUCED, deliberately
    # small kernel set: with all 21 domain kernels the components are still
    # estimable but the kernels are highly collinear, so the split between them
    # is not uniquely attributable and moves between near-equivalent solutions.
    # Five broad, mechanistically distinct blocks is the most that can be read.
    Kd <- list(clim  = k_linear(cell$X[, cols_in(cell$X, CLIM), drop = FALSE]),
               soil  = k_linear(cell$X[, cols_in(cell$X, SOIL), drop = FALSE]),
               emb   = k_linear(cell$X[, cols_in(cell$X, EMB),  drop = FALSE]),
               agri  = k_linear(cell$X[, cols_in(cell$X, AGRI), drop = FALSE]),
               space = k_spatial(cell$aux$lon, cell$aux$lat))
    Kd <- Kd[!vapply(Kd, is.null, TRUE)]
    if (length(Kd)) {
      v <- tryCatch(blup_varcomp(cell$y_mod, Kd, seq_len(cell$n)),
                    error = function(e) NULL)
      if (!is.null(v)) vc[[paste(cn, on, tgt)]] <-
        cbind(data.frame(country = cn, outcome = on, target = tgt,
                         n_areas = cell$n, stringsAsFactors = FALSE), v)
    }
    message(sprintf("  %-12s %-13s %-5s  %.0fs", cn, on, tgt,
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

# ── estimand C ──────────────────────────────────────────────────────────────
loco <- list()
LOCO_ARMS <- ARMS[c("null_train_mean", "domain_index", "blup_all", "blup_cs",
                    "blup_spatial", "blup_cs_spatial", "blup_5k")]
for (tgt in c("prev", "level")) {
  for (on in unique(ix$outcome)) {
    cl <- exp_all_cells(E, tgt, outcomes = on)
    if (length(cl) < 3) next
    loco[[paste(tgt, on)]] <- exp_loco(cl, LOCO_ARMS, domain_of = DOM)
    message("  LOCO ", tgt, " ", on)
  }
}

# ── write and summarise ─────────────────────────────────────────────────────
raw <- dplyr::bind_rows(rows)
exp_write(exp_summarise(raw), "01_kernel_blup_cells")
lc <- dplyr::bind_rows(loco)
exp_write(lc, "01_kernel_blup_loco")
exp_write(dplyr::bind_rows(vc), "01_kernel_varcomp")

sm <- exp_summarise(raw)
cat("\n== in-country: median Spearman over cells ==\n")
a <- aggregate(spearman ~ estimand + target + arm, data = sm,
               FUN = function(z) round(median(z, na.rm = TRUE), 3))
print(a[order(a$estimand, a$target, -a$spearman), ], row.names = FALSE)

cat("\n== transport (LOCO): mean Spearman over held-out cells ==\n")
b <- aggregate(spearman ~ target + arm, data = lc,
               FUN = function(z) round(mean(z, na.rm = TRUE), 3))
b$positive <- aggregate(spearman ~ target + arm, data = lc,
                        FUN = function(z) sum(z > 0, na.rm = TRUE))$spearman
b$cells <- aggregate(spearman ~ target + arm, data = lc,
                     FUN = function(z) sum(is.finite(z)))$spearman
print(b[order(b$target, -b$spearman), ], row.names = FALSE)

VC <- dplyr::bind_rows(vc)
if (nrow(VC)) {
  cat("\n== REML variance share by domain (median over cells, level target) ==\n")
  v <- VC[VC$target == "level", ]
  s <- aggregate(share ~ kernel, data = v, FUN = function(z) round(median(z), 3))
  print(s[order(-s$share), ], row.names = FALSE)
}
