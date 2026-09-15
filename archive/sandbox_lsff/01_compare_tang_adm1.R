# =====================================================================
# Compare this project's Admin-2 predictions against the modelled dietary
# inadequacy of Tang et al. (Nature Food 2026, doi:10.1038/s43016-026-01412-2).
#
# Two independent data sources for the same geography:
#   ours   = biomarker deficiency prevalence, predicted from proxy covariates
#   theirs = dietary inadequacy, modelled from HCES apparent consumption
# Levels are not comparable; RANK agreement is the question.
#
# Ghana        n = 10 pre-2018 regions (project Admin2 -> 16 -> 10, pop-weighted)
# Cote d'Ivoire n = 33 regions (paper ADM1 == project Admin2, direct join)
# =====================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

root <- ".."
ref  <- read.csv(file.path(root, "data/WFP_LSFF_2026/tangS1_adm1_vulnerability.csv"),
                 fileEncoding = "UTF-8")
xwalk <- read.csv(file.path(root, "data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv"),
                  fileEncoding = "UTF-8")

# normalise names for joining across accent/hyphen conventions
norm <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("autonome d'?", "", x)
  x <- gsub("[^a-z]", "", x)
  x
}

# which paper column corresponds to which project outcome
nutrient_of <- c(child_vitA = "vitA_pct", women_vitA = "vitA_pct",
                 child_iron = "iron_pct", women_iron = "iron_pct",
                 women_folate = "folate_pct", women_b12 = "b12_pct")

# Spearman + permutation p (exact-ish, robust at n=10)
spear <- function(a, b, nperm = 20000) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]
  n <- length(a)
  if (n < 4) return(list(n = n, rho = NA_real_, p = NA_real_))
  rho <- suppressWarnings(cor(a, b, method = "spearman"))
  set.seed(42)
  null <- replicate(nperm, cor(a, sample(b), method = "spearman"))
  list(n = n, rho = rho, p = (1 + sum(abs(null) >= abs(rho))) / (nperm + 1))
}

# ---------------------------------------------------------------- Ghana
pop  <- readRDS(file.path(root, "dashboard/data/admin2_population.rds"))
pred <- readRDS(file.path(root, "dashboard/data/admin2_area_predictions.rds"))

gha <- pred %>%
  filter(country == "Ghana") %>%
  left_join(pop %>% filter(country == "Ghana") %>% select(Admin2, pop_child, pop_women),
            by = "Admin2") %>%
  left_join(xwalk, by = c("Admin1" = "admin1_16"))

stopifnot(!any(is.na(gha$adm1_paper_10)))

gha_adm1 <- gha %>%
  mutate(w = ifelse(grepl("^child", outcome), pop_child, pop_women)) %>%
  filter(is.finite(w)) %>%
  group_by(outcome, adm1_paper_10) %>%
  summarise(
    n_admin2      = dplyr::n(),
    n_admin2_obs  = sum(is.finite(obs_prev)),
    pred_prev_agg = weighted.mean(pred_prev, w, na.rm = TRUE),
    obs_prev_agg  = if (any(is.finite(obs_prev))) {
                      weighted.mean(obs_prev[is.finite(obs_prev)], w[is.finite(obs_prev)])
                    } else NA_real_,
    .groups = "drop"
  ) %>%
  rename(pred_prev = pred_prev_agg, obs_prev = obs_prev_agg) %>%
  left_join(ref %>% filter(country == "Ghana"),
            by = c("adm1_paper_10" = "adm1_paper"))

# ------------------------------------------------------- Cote d'Ivoire
civ <- readRDS(file.path(root, "dashboard/data/oos_cote_divoire.rds"))$predictions
ref_civ <- ref %>% filter(country == "Cote dIvoire") %>% mutate(key = norm(adm1_paper))
civ <- civ %>% mutate(key = norm(Admin2))

unmatched <- setdiff(civ$key, ref_civ$key)
if (length(unmatched)) message("CIV unmatched: ", paste(unmatched, collapse = ", "))
civ_j <- left_join(civ, ref_civ, by = "key")
stopifnot(sum(is.na(civ_j$mpi)) == 0)

# ------------------------------------------------------------- results
rows <- list()
for (oc in unique(gha_adm1$outcome)) {
  d <- gha_adm1 %>% filter(outcome == oc)
  nut <- nutrient_of[[oc]]
  for (cmp in c("mpi", nut)) {
    s <- spear(d$pred_prev, d[[cmp]])
    rows[[length(rows) + 1]] <- data.frame(
      country = "Ghana", unit = "10 pre-2018 regions", outcome = oc,
      ours = "model prediction", theirs = cmp, n = s$n, rho = s$rho, p = s$p)
    if (any(is.finite(d$obs_prev))) {
      s2 <- spear(d$obs_prev, d[[cmp]])
      rows[[length(rows) + 1]] <- data.frame(
        country = "Ghana", unit = "10 pre-2018 regions", outcome = oc,
        ours = "observed survey", theirs = cmp, n = s2$n, rho = s2$rho, p = s2$p)
    }
  }
}
for (oc in unique(civ_j$outcome)) {
  d <- civ_j %>% filter(outcome == oc)
  nut <- nutrient_of[[oc]]
  for (cmp in c("mpi", nut)) {
    s <- spear(d$pred_prev, d[[cmp]])
    rows[[length(rows) + 1]] <- data.frame(
      country = "Cote d'Ivoire", unit = "33 regions (= our Admin2)", outcome = oc,
      ours = "model prediction (out-of-sample)", theirs = cmp,
      n = s$n, rho = s$rho, p = s$p)
  }
}
res <- bind_rows(rows)
# 28 tests over correlated outcomes -- report FDR, not raw p
res$q <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
res$rho <- round(res$rho, 3); res$p <- round(res$p, 4); res$q <- round(res$q, 3)

dir.create(file.path(root, "results/tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(root, "results/tables/tang_lsff_adm1_agreement.csv"), row.names = FALSE)
write.csv(gha_adm1, "gha_adm1_joined.csv", row.names = FALSE)
write.csv(civ_j %>% select(-key), "civ_adm2_joined.csv", row.names = FALSE)

cat("\n===== RANK AGREEMENT: ours vs Tang et al. dietary inadequacy =====\n\n")
print(res, row.names = FALSE)
cat("\nGhana Admin-2 units with observed biomarker data, by region:\n")
print(gha_adm1 %>% filter(outcome == "child_vitA") %>%
        select(adm1_paper_10, n_admin2, n_admin2_obs, pred_prev, obs_prev, mpi, vitA_pct) %>%
        as.data.frame(), row.names = FALSE)
