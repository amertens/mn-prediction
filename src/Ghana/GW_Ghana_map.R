################################################################################
# Packages ---------------------------------------------------------------------
################################################################################
libs <- c("tidyverse",      # data wrangling + ggplot2
          "here",
          "sf",             # simple‑features objects
          "rnaturalearth",  # country borders
          "rnaturalearthdata",
          "ggspatial")      # nice scale bars / north arrows
invisible(lapply(libs, require, character.only = TRUE))

################################################################################
# 1.  Get a Ghana boundary polygon ---------------------------------------------
################################################################################
ghana <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf",
                                     country = "Ghana") |>
  st_transform(4326)              # WGS‑84 (lon/lat) – same as your GPS

################################################################################
# 2.  Convert the GPS data to an sf object -------------------------------------
################################################################################

gps <- read.csv(here("data", "IPD", "Ghana", "Ghana_GMS_GPS_cleaned.csv"))

clusters_sf <- gps |>
  distinct(EACode, longitude, latitude, urban..rural, Households) |>   # one point per EA
  st_as_sf(coords = c( "longitude","latitude"), crs = 4326)

################################################################################
# 3.  Quick static map with ggplot2 --------------------------------------------
################################################################################
p <- ggplot(ghana) +
  geom_sf(fill = "grey95", colour = "grey40", linewidth = .3) +
  geom_sf(data = clusters_sf,
          aes(colour = urban..rural, size=Households),      # colour by Urban / Rural
          alpha = .8) +
  scale_colour_manual(values = c("urban" = "#377eb8",
                                 "rural" = "#4daf4a"),
                      name = "Cluster type") +
  annotation_scale(location = "bl", width_hint = 0.4, text_cex = .6) +
  annotation_north_arrow(location = "tl",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(expand = FALSE) +
  labs(title = "Ghana: 2017 survey cluster locations",
       caption = "") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(colour = "transparent"))
p

saveRDS(p, file=here("figures", "GW_Ghana_map.rds"))
