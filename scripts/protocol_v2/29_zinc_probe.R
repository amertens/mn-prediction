# =============================================================================
# scripts/protocol_v2/29_zinc_probe.R   [ZN-01]
#
# ZINC: A RELIABLE TARGET THE PROXIES DO NOT TOUCH. WHY?
#
# HR-01 found the zinc outcomes are the extreme case in the headroom map: the
# empirical split-half ceiling is ~0.8 (the between-unit variation is REAL,
# not sampling noise), yet the in-fill model scores below zero. Every other
# outcome has some proxy signal; zinc has none. Zinc is measured in Malawi
# only (87 Traditional Authorities), so nothing here can be replicated across
# countries -- this is a diagnosis, not a claim.
#
# Four checks:
#   1. Is the zinc target coherent? child_zinc vs women_zinc rank correlation
#      across the same units. Two independent samples of the same places
#      agreeing means the geography is real.
#   2. Is zinc orthogonal to the agro-ecology axis the model learns? Rank
#      correlation of zinc prevalence with the other Malawi outcomes.
#   3. Does ANY indicator correlate with zinc within Malawi? All predictors,
#      Spearman, with a max-|r| permutation null for family-wise error, and
#      the same scan for child iron as a positive control.
#   4. Soil zinc specifically: sign and strength of every soil-zinc column.
#
#   Rscript scripts/protocol_v2/29_zinc_probe.R
# -> results/tables/protocol_v2/zinc_probe.csv
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
setwd("C:/Users/andre/OneDrive/Documents/mn-prediction")
source("R/protocol_v2.R")
OUTDIR <- "results/tables/protocol_v2"; NPERM <- 999L; set.seed(20260903L)
TG <- read.csv(file.path(OUTDIR, "targets_v2.csv"), stringsAsFactors = FALSE)
S  <- read.csv("data/covariates/harmonized/predictors_admin2_shared.csv", check.names = FALSE)
MD <- read.csv("data/covariates/harmonized/predictors_admin2_shared_metadata.csv")
PREDS <- drop_near_outcome_v2(intersect(MD$column, names(S)), MD)
domain_of <- stats::setNames(MD$domain, MD$column)
E <- tryCatch(read.csv("results/tables/reliability_empirical.csv"), error = function(e) NULL)

mw <- TG[TG$country == "Malawi", ]
wide <- mw |> select(Admin2, outcome, y_prev) |> tidyr::pivot_wider(names_from = outcome, values_from = y_prev)
cat("Malawi units:", nrow(wide), "| outcomes:", paste(setdiff(names(wide), "Admin2"), collapse = ", "), "\n")

cat("\n===== 1. ceiling and coherence of the zinc target =====\n")
if (!is.null(E)) { e <- E[E$country == "Malawi" & grepl("zinc", E$outcome) & E$scheme == "within", c("outcome", "r_max_emp", "r_max_emp_lo", "r_max_emp_hi")]
  if (nrow(e)) print(e, row.names = FALSE) }
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 10) return(NA_real_); suppressWarnings(cor(a[ok], b[ok], method = "spearman")) }
cat(sprintf("child_zinc ~ women_zinc across units: Spearman %+.3f\n", sp(wide$child_zinc, wide$women_zinc)))

cat("\n===== 2. is zinc on the same axis as the other outcomes? (Spearman across units) =====\n")
oth <- setdiff(names(wide), c("Admin2", "child_zinc", "women_zinc"))
print(data.frame(outcome = oth, child_zinc = round(vapply(oth, function(o) sp(wide$child_zinc, wide[[o]]), numeric(1)), 3),
                 women_zinc = round(vapply(oth, function(o) sp(wide$women_zinc, wide[[o]]), numeric(1)), 3)), row.names = FALSE)

cat("\n===== 3. bivariate scan within Malawi, family-wise permutation null (max |r|) =====\n")
scan <- function(on) {
  t <- mw[mw$outcome == on, c("Admin1", "Admin2", "y_prev")]
  m <- inner_join(t, S[S$country == "Malawi", c("Admin1", "Admin2", PREDS)], by = c("Admin1", "Admin2"))
  y <- m$y_prev; X <- as.matrix(m[, PREDS]); keep <- apply(X, 2, function(z) sum(is.finite(z)) > 0.8 * length(z) && sd(z, na.rm = TRUE) > 0)
  X <- X[, keep, drop = FALSE]; yr <- rank(y); Xr <- apply(X, 2, function(z) { z[!is.finite(z)] <- median(z, na.rm = TRUE); rank(z) })
  r <- suppressWarnings(cor(Xr, yr)); r <- as.numeric(r)
  mx <- replicate(NPERM, max(abs(suppressWarnings(cor(Xr, sample(yr)))), na.rm = TRUE))
  fwer <- vapply(abs(r), function(a) (1 + sum(mx >= a)) / (NPERM + 1), numeric(1))
  data.frame(outcome = on, predictor = colnames(X), domain = unname(domain_of[colnames(X)]), r = r, p_fwer = fwer, stringsAsFactors = FALSE)
}
res <- bind_rows(lapply(intersect(c("child_zinc", "women_zinc", "child_iron"), unique(mw$outcome)), scan))
write.csv(res, file.path(OUTDIR, "zinc_probe.csv"), row.names = FALSE)
for (on in unique(res$outcome)) {
  d <- res[res$outcome == on, ]; d <- d[order(-abs(d$r)), ]
  cat(sprintf("\n-- %s : %d predictors | family-wise significant (p<0.05): %d | top 8 --\n", on, nrow(d), sum(d$p_fwer < 0.05)))
  print(head(data.frame(predictor = d$predictor, domain = substr(d$domain, 1, 28), r = round(d$r, 3), p_fwer = round(d$p_fwer, 3)), 8), row.names = FALSE)
}

cat("\n===== 4. soil zinc columns vs zinc outcomes =====\n")
sz <- res[grepl("zinc", res$predictor, ignore.case = TRUE) & grepl("zinc", res$outcome), ]
if (nrow(sz)) print(data.frame(outcome = sz$outcome, predictor = sz$predictor, r = round(sz$r, 3), p_fwer = round(sz$p_fwer, 3)), row.names = FALSE) else cat("  no soil-zinc predictor columns present\n")
cat("\nDONE\n")
