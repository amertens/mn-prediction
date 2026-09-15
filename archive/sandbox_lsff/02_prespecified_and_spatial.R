# Restrict to the PRE-SPECIFIABLE family (our outcome vs the SAME nutrient's
# % vulnerable), then ask the question this project has already learned to ask:
# is any agreement just a shared north-south gradient?
suppressPackageStartupMessages({library(dplyr); library(sf)})
root <- ".."

gha <- read.csv("gha_adm1_joined.csv", fileEncoding = "UTF-8")
civ <- read.csv("civ_adm2_joined.csv", fileEncoding = "UTF-8")
nut <- c(child_vitA="vitA_pct", women_vitA="vitA_pct", child_iron="iron_pct",
         women_iron="iron_pct", women_folate="folate_pct", women_b12="b12_pct")

# ---- centroid latitude ------------------------------------------------
civ_b <- readRDS(file.path(root,"dashboard/data/oos_cote_divoire.rds"))$boundaries
suppressWarnings({ cc <- st_coordinates(st_centroid(st_geometry(civ_b))) })
civ_lat <- data.frame(Admin2 = civ_b$Admin2, lat = cc[,2])
civ <- left_join(civ, civ_lat, by = "Admin2")

a1 <- readRDS(file.path(root,"dashboard/data/admin1_boundaries.rds"))
g1 <- a1[["ghana"]]
suppressWarnings({ gc_ <- st_coordinates(st_centroid(st_geometry(g1))) })
xw <- read.csv(file.path(root,"data/WFP_LSFF_2026/ghana_region_crosswalk_16_to_10.csv"),
               fileEncoding = "UTF-8")
gha_lat <- data.frame(admin1_16 = g1$Admin1, lat = gc_[,2]) %>%
  left_join(xw, by = "admin1_16") %>%
  group_by(adm1_paper_10) %>% summarise(lat = mean(lat), .groups="drop")
gha <- left_join(gha, gha_lat, by = "adm1_paper_10")

# ---- Spearman, and Spearman partialling out latitude ------------------
perm_p <- function(stat_fun, a, b, nperm = 20000) {
  set.seed(42); s <- stat_fun(a, b)
  null <- replicate(nperm, stat_fun(a, sample(b)))
  c(stat = s, p = (1 + sum(abs(null) >= abs(s))) / (nperm + 1))
}
rho_fun <- function(a,b) cor(a, b, method="spearman")

out <- list()
for (df_name in c("gha","civ")) {
  d0 <- get(df_name)
  for (oc in unique(d0$outcome)) {
    d <- d0 %>% filter(outcome == oc)
    y <- d$pred_prev; x <- d[[ nut[[oc]] ]]; lat <- d$lat
    ok <- is.finite(y) & is.finite(x) & is.finite(lat)
    y <- y[ok]; x <- x[ok]; lat <- lat[ok]
    r  <- perm_p(rho_fun, y, x)
    # partial: rank-residualise both on latitude
    ry <- resid(lm(rank(y) ~ rank(lat))); rx <- resid(lm(rank(x) ~ rank(lat)))
    rp <- perm_p(function(a,b) cor(a,b,method="pearson"), ry, rx)
    out[[length(out)+1]] <- data.frame(
      country = ifelse(df_name=="gha","Ghana","Cote d'Ivoire"),
      n = length(y), outcome = oc, theirs = nut[[oc]],
      rho = r[["stat"]], p = r[["p"]],
      rho_partial_lat = rp[["stat"]], p_partial = rp[["p"]],
      rho_lat_ours = cor(y, lat, method="spearman"),
      rho_lat_theirs = cor(x, lat, method="spearman"))
  }
}
res <- bind_rows(out)
res$q <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], round, 3)
write.csv(res, file.path(root,"results/tables/tang_lsff_prespecified.csv"), row.names=FALSE)
cat("\n== PRE-SPECIFIED nutrient-matched family (our prediction vs same nutrient) ==\n")
cat("   rho_lat_* : correlation of each source with centroid latitude\n\n")
print(res, row.names = FALSE)
