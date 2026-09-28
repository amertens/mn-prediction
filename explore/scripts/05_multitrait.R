# =============================================================================
# explore/scripts/05_multitrait.R   [probe MT-01]
#
# QUESTION. XO-01 found that borrowing ACROSS nutrients is null under transport,
# but that the SAME nutrient in the other population (women <-> children) beats
# the index (0.454 vs 0.403). That is a low-rank structure across the 24 cells,
# and the project currently fits every cell on its own.
#
# Quantitative genetics fits correlated traits jointly (multi-trait BLUP), and
# gains most exactly where single-trait estimation is noisiest - which at
# n = 14-87 is everywhere.
#
# WHAT BORROWING CAN AND CANNOT DO HERE. A held-out district has NO biomarker
# measured, so the model cannot condition on its other nutrients. The gain must
# come from estimating the SHARED spatial component on more data: all traits of
# the training districts, not one.
#
#   general factor g   = loadings-weighted mean of the country's standardised
#                        traits on TRAINING districts (loadings = first
#                        eigenvector of the training trait correlation matrix)
#   mt_blup            = BLUP of g on the kernel, scaled back by trait loading,
#                        plus a trait-specific BLUP of the residual
#   st_blup            = the same kernel, single trait (the controlled contrast)
#
# The held-out districts are held out for EVERY trait simultaneously, which the
# script asserts and prints - otherwise the general factor reads the held-out
# district's own survey through another biomarker.
#
#   Rscript explore/scripts/05_multitrait.R
# -> explore/out/05_multitrait_cells.csv, 05_multitrait_traitcor.csv
# =============================================================================
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/harness.R")
source("C:/Users/andre/OneDrive/Documents/mn-prediction/explore/R/methods_kernel.R")

REPS <- as.integer(Sys.getenv("EXP_REPS", "10"))
E <- exp_load()
DOM <- E$domain_of
CS <- c("Climate and weather", "Soil characteristics")
ix <- exp_cell_index(E)

#' All outcomes of one country on the districts they share
country_traits <- function(cn, tgt) {
  ons <- unique(ix$outcome[ix$country == cn])
  cs <- list()
  for (on in ons) {
    cc <- tryCatch(exp_cell(E, cn, on, tgt), error = function(e) NULL)
    if (!is.null(cc)) {
      cc$key <- paste(cc$Admin1, cc$Admin2, sep = "|")
      cs[[on]] <- cc
    }
  }
  if (length(cs) < 2) return(NULL)
  keys <- Reduce(intersect, lapply(cs, function(z) z$key))
  if (length(keys) < 15) return(NULL)
  Y <- do.call(cbind, lapply(cs, function(z) {
    v <- z$y_mod[match(keys, z$key)]
    as.numeric(scale(v))
  }))
  colnames(Y) <- names(cs)
  ref <- cs[[1]]
  i0 <- match(keys, ref$key)
  list(Y = Y, keys = keys, X = ref$X[i0, , drop = FALSE],
       lon = ref$aux$lon[i0], lat = ref$aux$lat[i0],
       Admin1 = ref$Admin1[i0], n = length(keys),
       cells = lapply(cs, function(z) {
         i <- match(keys, z$key)
         list(y_nat = z$y_nat[i], y_mod = z$y_mod[i], w = z$w[i])
       }))
}

rows <- list(); tcor <- list()
for (cn in unique(ix$country)) {
  for (tgt in c("prev", "level")) {
    ct <- country_traits(cn, tgt)
    if (is.null(ct)) { message("  skip ", cn, " ", tgt); next }
    n <- ct$n; TT <- ncol(ct$Y)
    cc <- colnames(ct$X)[which(DOM[colnames(ct$X)] %in% CS)]
    Kl <- list(cs = k_linear(ct$X[, cc, drop = FALSE]),
               space = k_spatial(ct$lon, ct$lat))
    Kl <- Kl[!vapply(Kl, is.null, TRUE)]
    if (!length(Kl)) next

    if (tgt == "level")
      tcor[[cn]] <- cbind(country = cn, as.data.frame(as.table(stats::cor(ct$Y))))

    for (r in seq_len(REPS)) {
      folds <- make_folds_v2("kfold_district", n, k = 5, rep_id = r)
      pred_st <- pred_mt <- matrix(NA_real_, n, TT, dimnames = list(NULL, colnames(ct$Y)))
      for (f in unique(folds)) {
        te <- which(folds == f); tr <- which(folds != f)
        if (length(tr) < 15 || !length(te)) next
        # ASSERTION: the held-out rows are held out for every trait at once
        stopifnot(length(intersect(tr, te)) == 0)

        Ytr <- ct$Y[tr, , drop = FALSE]
        R <- suppressWarnings(stats::cor(Ytr, use = "pairwise.complete.obs"))
        R[!is.finite(R)] <- 0
        eg <- eigen(R, symmetric = TRUE)
        l <- eg$vectors[, 1]
        if (mean(l > 0) < 0.5) l <- -l           # majority-positive orientation
        g_tr <- as.numeric(Ytr %*% l) / sum(l^2)

        gfull <- rep(NA_real_, n); gfull[tr] <- g_tr
        gh <- reml_blup(ifelse(is.na(gfull), 0, gfull), Kl, tr, te)$pred

        for (k in seq_len(TT)) {
          yk <- ct$Y[, k]
          pred_st[te, k] <- reml_blup(yk, Kl, tr, te)$pred
          resid <- rep(NA_real_, n); resid[tr] <- yk[tr] - l[k] * g_tr
          rh <- reml_blup(ifelse(is.na(resid), 0, resid), Kl, tr, te)$pred
          pred_mt[te, k] <- l[k] * gh + rh
        }
      }
      for (k in seq_len(TT)) {
        on <- colnames(ct$Y)[k]
        cel <- ct$cells[[on]]
        sc <- if (tgt == "prev") "prev" else "level"
        # predictions are on the within-country standardised scale; map back
        back <- function(p) p * stats::sd(cel$y_mod, na.rm = TRUE) +
          mean(cel$y_mod, na.rm = TRUE)
        for (a in c("st_blup", "mt_blup")) {
          p <- back(if (a == "st_blup") pred_st[, k] else pred_mt[, k])
          if (tgt == "prev") p <- .v2_expit(p)
          s <- score_v2(cel$y_nat, p, cel$w, scale = sc)
          rows[[paste(cn, tgt, on, a, r)]] <- cbind(
            data.frame(country = cn, outcome = on, target = tgt,
                       estimand = "infill", arm = a, rep = r, n_areas = n,
                       n_traits = TT, stringsAsFactors = FALSE), s)
        }
      }
    }
    message(sprintf("  %-12s %-5s  %d districts x %d traits", cn, tgt, n, TT))
  }
}

raw <- dplyr::bind_rows(rows)
sm <- exp_summarise(raw); exp_write(sm, "05_multitrait_cells")
TC <- dplyr::bind_rows(tcor)
if (nrow(TC)) { names(TC) <- c("country", "trait1", "trait2", "cor")
  exp_write(TC[TC$trait1 != TC$trait2, ], "05_multitrait_traitcor") }

cat("\n== single-trait vs multi-trait BLUP, paired by cell ==\n")
w <- reshape(sm[, c("country", "outcome", "target", "arm", "spearman")],
             idvar = c("country", "outcome", "target"), timevar = "arm",
             direction = "wide")
names(w) <- sub("^spearman\\.", "", names(w))
w$gain <- round(w$mt_blup - w$st_blup, 3)
print(w[order(-w$gain), ], row.names = FALSE, digits = 3)
for (tgt in c("level", "prev")) {
  s <- w[w$target == tgt, ]
  cat(sprintf("\n%-5s  mean gain %+.4f   improved %d of %d cells\n", tgt,
              mean(s$gain, na.rm = TRUE), sum(s$gain > 0, na.rm = TRUE),
              sum(is.finite(s$gain))))
}

if (nrow(TC)) {
  cat("\n== strongest trait correlations within country (level) ==\n")
  T2 <- TC[TC$trait1 != TC$trait2, ]
  T2 <- T2[!duplicated(t(apply(T2[, c("trait1", "trait2")], 1, sort))), ]
  print(head(T2[order(-abs(T2$cor)), ], 15), row.names = FALSE, digits = 3)
}
