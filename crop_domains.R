
library(readr)
library(sp)
library(raster)
library(gstat)
library(rgdal)
library(RNetCDF)
library(ncdf4)
library(stringr)

library(raster)
library(leaflet)
library(htmlwidgets)

# setwd("D:/Dust_Event_UAE_2015/WRF_trial_runs")
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRF_trial_runs")
setwd("D:/Dust_Event_UAE_2015/WRF_trial_runs/DUST_AOD_FK/extinction55/12km")

# get extent from WRFChem output-----------------------------------------------
# use the raster stack with all the WRF_chem output

# get only one image
# WRF_STACK_image <- raster("DUST_WRFChem_02April2015_stack_6_DAYS_LARGE.tif", band = 2)
# WRF_STACK_image_d01_12km <- raster("AOD_12km_WRFChem_DUST300.tif", band =10)
WRF_STACK_image_d01_12km <- raster("AOD_WRFChem_02April2015_aod_dust_opt3_chem_opt300.tif", band =10)
plot(WRF_STACK_image_d01_12km)

# get the extent of the raster
e <- extent(WRF_STACK_image_d01_12km)
plot(e)
# make a spatial polygon from the extent
p <- as(e, "SpatialPolygons")
plot(p)
crs(p) <- "proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# save shp file for the rectangular DOMAIN from WRFChem
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRFChem_domain")
shapefile(p, "domain_d01_MIKE.shp", overwrite=TRUE)


#########################################################################
#########################################################################

#### importing the NEW ADMINISTRATIVE WORLDWIDE shapefile to use as a masking 
dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS"
### shapefile for WW
shp_WW <- readOGR(dsn = dir, layer = "new_ADMIN_00")
# ----- Transform to EPSG 4326 - WGS84 (required)
shp_WW <- spTransform(shp_WW, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
shp_WW@data$name <- 1:nrow(shp_WW)
# plot(shp_WW)

d01_shp_WW <- crop(shp_WW, e)
plot(d01_shp_WW)
# save shp file for the rectangular DOMAIN from WRFChem
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRFChem_domain")
shapefile(d01_shp_WW, "ADMIN_domain_MIKE.shp.shp", overwrite=TRUE)



########################################################################
########################################################################
########################################################################

##### use smaller WRF chem domain (for the 4km resolution data) ######################################

# get only one image
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRF_trial_runs")
setwd("D:/Dust_Event_UAE_2015/WRF_trial_runs/DUST_AOD_FK/extinction55/4km")
# WRF_STACK_image <- raster("AOD_WRFChem_02April2015_stack_5_DAYS.tif", band = 5)
WRF_STACK_image_d02_4km <- raster("AOD_4km_WRFChem_DUST300.tif", band =10)
plot(WRF_STACK_image_d02_4km)

# get the extent of the raster
e <- extent(WRF_STACK_image_d02_4km)
plot(e)
# make a spatial polygon from the extent
p <- as(e, "SpatialPolygons")
plot(p)
crs(p) <- "proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# save shp file for the rectangular DOMAIN from WRFChem
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRFChem_domain")
shapefile(p, "domain_d02_4km_WRFChem_small.shp", overwrite=TRUE)


#########################################################################
#########################################################################

#### importing the NEW ADMINISTRATIVE WORLDWIDE shapefile to use as a masking 
dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS"
### shapefile for WW
shp_WW <- readOGR(dsn = dir, layer = "new_ADMIN_00")
# ----- Transform to EPSG 4326 - WGS84 (required)
shp_WW <- spTransform(shp_WW, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
shp_WW@data$name <- 1:nrow(shp_WW)
# plot(shp_WW)

d01_shp_WW <- crop(shp_WW, e)
plot(d01_shp_WW)
# save shp file for the rectangular DOMAIN from WRFChem
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRFChem_domain")
shapefile(d01_shp_WW, "ADMIN_domain_d01_4km_WRFChem.shp", overwrite=TRUE)




