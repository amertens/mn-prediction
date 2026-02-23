

pred_mat = readRDS(here("data","IPD","Ghana","GMS_cluster_prevalence_bootdraws.rds"))


# -------------------- 5) Build a 5-km prediction grid over Ghana --------------------
ghana_ll  <- rnaturalearth::ne_countries(scale = "medium", country = "Ghana", returnclass = "sf") %>% st_transform(4326)
ghana_utm <- st_transform(ghana_ll, 32630)

grid_centers <- st_make_grid(ghana_utm, cellsize = 5000, what = "centers") %>% st_sf()
grid_sf      <- st_intersection(grid_centers, ghana_utm)  # keep inside Ghana

# predictors for the grid (lon/lat + any area-level covariates you later add)
grid_ll        <- st_transform(grid_sf, 4326)
grid_coords_ll <- st_coordinates(grid_ll)
X_grid <- data.frame(
  lon = grid_coords_ll[,1],
  lat = grid_coords_ll[,2]
)


q_lo <- apply(pred_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
q_hi <- apply(pred_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
mu   <- colMeans(pred_mat, na.rm = TRUE)
sdv  <- apply(pred_mat, 2, sd, na.rm = TRUE)
thr  <- 0.20
exc  <- colMeans(pred_mat > thr, na.rm = TRUE)

pred_unc_sf <- sf::st_transform(grid_sf, 4326) |>
  dplyr::mutate(
    pred_mean = mu,
    pred_sd   = sdv,
    pred_lo   = q_lo,
    pred_hi   = q_hi,
    exc_20    = exc
  )



# Example: quick maps
library(ggplot2)
p_mean <- ggplot() +
  geom_sf(data = pred_unc_sf, aes(color = pred_mean), size = 2) +
  scale_color_viridis_c("Mean prevalence", limits = c(0, 0.7)) +
  coord_sf(expand = FALSE) + theme_minimal() + ggtitle("Predicted Child VitA Deficiency")

p_wid <- ggplot() +
  geom_sf(data = pred_unc_sf, aes(color = pred_hi - pred_lo), size = 1.5) +
  scale_color_viridis_c("95% CI width") +
  coord_sf(expand = FALSE) + theme_minimal() + ggtitle("Spatial Uncertainty in Child VitA Deficiency Predictions")

p_exc <- ggplot() +
  geom_sf(data = pred_unc_sf, aes(color = exc_20), size = 2) +
  scale_color_viridis_c("Pr(VAD > 20%)", limits = c(0, 1)) +
  coord_sf(expand = FALSE) + theme_minimal()

p_mean; p_wid; p_exc


# # --- (Optional) aggregate to ADM2 with uncertainty ----------
# poly.adm <- geodata::gadm(country="GH", level=2, path=tempdir())
# poly.adm <- sf::st_as_sf(poly.adm) %>% select(NAME_1, NAME_2) %>% rename(Admin1 = NAME_1, Admin2 = NAME_2)
# adm2 <- st_transform(poly.adm, crs = 4326)
# # For each bootstrap draw, area- or population-weight grid points, then summarize:
# # (Here simple area-weight: each grid point equal weight within polygons)
# M <- nrow(pred_mat);
# G <- nrow(pred_unc_sf)
# adm2_ids <- sf::st_within(pred_unc_sf, adm2) %>% lengths() # sanity
# idx_list <- sf::st_within(pred_unc_sf, adm2)
# agg_draws <- lapply(idx_list, function(ids) ids) # placeholder
# #Better: build a sparse mapping matrix and do matrix multiplies for speed.
