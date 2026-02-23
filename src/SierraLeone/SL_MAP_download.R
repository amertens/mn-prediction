
rm(list=ls())
library(here)
library(terra)
library(malariaAtlas)

#all rasters to get
map_rasters <- listRaster()

#getting a country specific shapefile
MDG_shp <- getShp(ISO = "SLE", admin_level = "admin0")
autoplot(MDG_shp)

for(i in 1:nrow(map_rasters)){
  rast=NULL
  rast = try(getRaster(dataset_id = map_rasters$dataset_id[i], shp = MDG_shp, year = 2012))
  if(!is.null(rast) & class(rast)!="try-error"){
    writeRaster(rast, here(paste0("data/Malaria Atlas/SierraLeone/",map_rasters$dataset_id[i],".tif")), overwrite=TRUE)
  }
}
