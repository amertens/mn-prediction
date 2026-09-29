# =============================================================================
# explore/scripts/22_zinc_why.R   [probe ZW-01]
# Why is zinc unpredicted despite high split-half reliability?
# Hypothesis: its reproducible between-"district" variation is CLUSTER-level.
# With 85% of Malawi districts being a single cluster, a district-mean split-half
# counts the cluster effect as district signal, so r_max is inflated.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
EXP_ROOT <- "C:/Users/andre/OneDrive/Documents/mn-prediction"; setwd(EXP_ROOT)
source("explore/R/harness.R"); source("R/config.R"); source("R/data_prep.R")
num <- function(x) suppressWarnings(as.numeric(haven::zap_labels(x)))

#' nested variance components by moment matching: district, cluster-within-
#' district, residual. Balanced-ish design, so a two-stage ANOVA is adequate.
vc <- function(y, dis, cl) {
  d <- data.frame(y = y, dis = factor(dis), cl = factor(paste(dis, cl)))
  d <- d[is.finite(d$y), ]
  if (nrow(d) < 50) return(NULL)
  fit <- tryCatch(stats::aov(y ~ dis + cl, d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  a <- summary(fit)[[1]]
  ms <- a[["Mean Sq"]]; df <- a[["Df"]]
  nm <- trimws(rownames(a))
  ms_d <- ms[nm == "dis"]; ms_c <- ms[nm == "cl"]; ms_e <- ms[nm == "Residuals"]
  n_per_cl <- nrow(d) / dplyr::n_distinct(d$cl)
  n_per_di <- nrow(d) / dplyr::n_distinct(d$dis)
  s2e <- ms_e
  s2c <- max((ms_c - ms_e) / n_per_cl, 0)
  s2d <- max((ms_d - ms_c) / n_per_di, 0)
  tot <- s2d + s2c + s2e
  data.frame(var_district = s2d, var_cluster = s2c, var_resid = s2e,
             share_district = s2d/tot, share_cluster = s2c/tot,
             n = nrow(d), clusters = dplyr::n_distinct(d$cl),
             districts = dplyr::n_distinct(d$dis),
             pct_single_cluster = 100*mean(table(d$dis[!duplicated(paste(d$dis,d$cl))]) < 2))
}

CFG <- get_country_configs()
rows <- list()
for (cn in names(CFG)) {
  cc <- CFG[[cn]]
  dat <- tryCatch(load_merged_data(cc$data_path), error = function(e) NULL)
  if (is.null(dat)) next
  a2 <- cc$admin2_col; psu <- cc$psu_col
  if (!all(c(a2, psu) %in% names(dat))) next
  for (on in names(cc$outcomes)) {
    oc <- cc$outcomes[[on]]
    if (is.null(oc$continuous) || !oc$continuous %in% names(dat)) next
    keep <- outcome_population_mask(dat, cc, oc, label = "[ZW-01]")
    d <- dat[keep, ]
    y <- num(d[[oc$continuous]])
    if (!grepl("log", oc$continuous, ignore.case = TRUE) && all(y[is.finite(y)] > 0)) y <- log(y)
    v <- vc(y, as.character(d[[a2]]), as.character(d[[psu]]))
    if (is.null(v)) next
    rows[[paste(cc$country,on)]] <- cbind(
      data.frame(country = cc$country, outcome = on, stringsAsFactors = FALSE), v)
  }
}
V <- dplyr::bind_rows(rows)
V$cluster_to_district <- round(V$share_cluster / pmax(V$share_district, 1e-6), 1)
exp_write(V, "22_zinc_variance_components")

cat("\n== nested variance components of the continuous biomarker ==\n")
cat("   share_district = real Admin-2 geography; share_cluster = within-district\n")
cat("   cluster effect, which Admin-2 covariates cannot predict\n\n")
p <- V[order(-V$share_district), c("country","outcome","districts","clusters",
       "pct_single_cluster","share_district","share_cluster","cluster_to_district")]
print(p, row.names = FALSE, digits = 3)

cat("\n== does the district share explain what the model achieves? ==\n")
BM <- read.csv("results/tables/protocol_v2/benchmarks_v2_cells.csv", stringsAsFactors=FALSE)
b <- BM[BM$estimand=="infill" & BM$target=="level" & BM$arm=="domain_index",
        c("country","outcome","spearman")]; names(b)[3] <- "achieved"
V2 <- V; V2$country[V2$country=="Sierra Leone"] <- "SierraLeone"
m <- merge(V2, b, by=c("country","outcome"))
R <- read.csv(file.path(EXP_OUT,"16_reliability_all.csv"), stringsAsFactors=FALSE)
R$country[R$country=="Sierra Leone"] <- "SierraLeone"
m <- merge(m, R[,c("country","outcome","r_full_respondent")], by=c("country","outcome"))
m$r_max_splithalf <- sqrt(pmax(m$r_full_respondent,0))
m$r_max_district  <- sqrt(pmax(m$share_district/(m$share_district+m$share_cluster),0))
for (v in c("share_district","r_full_respondent")) {
  ct <- suppressWarnings(cor.test(m[[v]], m$achieved, method="spearman"))
  cat(sprintf("  spearman(%-18s, achieved) = %+.3f  p = %.4f  (n=%d)\n", v, ct$estimate, ct$p.value, nrow(m)))
}
cat("\n== zinc against the well-predicted cells ==\n")
print(m[order(-m$share_district), c("country","outcome","share_district","share_cluster",
      "r_full_respondent","achieved")], row.names=FALSE, digits=3)
