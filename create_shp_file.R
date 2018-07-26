

library(readr)
library(sp)
library(raster)
library(gstat)
library(rgdal)
library(RNetCDF)
library(ncdf4)
library(stringr)
library(rgeos)
library(leaflet)
library(htmlwidgets)

e <- extent(50,60,20, 28)  # UAE extent
plot(e)

# make a spatial polygon from the extent
p <- as(e, "SpatialPolygons")
plot(p)
proj4string(p) <- CRS("+proj=longlat +datum=WGS84")
# crs(p) <- "proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# save shp file for the rectangular DOMAIN
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS/prova_shapes")
shapefile(p, "rectangular_domain.shp", overwrite=TRUE)

# reload and plot domain

dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS/prova_shapes"
shp_rect <- readOGR(dsn = dir, layer = "rectangular_domain")
# ----- Transform to EPSG 4326 - WGS84 (required)
shp_rect <- spTransform(shp_rect, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
shp_rect@data$name <- 1:nrow(shp_rect)
plot(shp_rect)



dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/HISTORICAL_dust/UAE_boundary"
### shapefile for UAE
shp_UAE <- readOGR(dsn = dir, layer = "uae_emirates")

# ----- Transform to EPSG 4326 - WGS84 (required)
shp_UAE <- spTransform(shp_UAE, CRS("+init=epsg:4326"))
# names(shp)

shp_UAE@data$name <- 1:nrow(shp_UAE)
plot(shp_rect)
plot(shp_UAE, add=TRUE, lwd=1)


# d01_shp_WW <- crop(shp_WW, e)
# plot(d01_shp_WW)
# setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/Dust_Event_UAE_2015/WRFChem_domain")
# shapefile(d01_shp_WW, "ADMIN_domain_MIKE.shp.shp", overwrite=TRUE)



require(sf)
coordinates_ABU_DHABI <- read.table(text="
                    longitude    latitude
                    54.646    24.43195",
                    header=TRUE)

coord_ABU_DHABI_point <- st_as_sf(x = coordinates_ABU_DHABI, 
                        coords = c("longitude", "latitude"),
                        crs = "+proj=longlat +datum=WGS84")
# simple plot
plot(coord_ABU_DHABI_point)

# convert to sp object if needed
coord_ABU_DHABI_point <- as(coord_ABU_DHABI_point, "Spatial")

# shp_buff <- gBuffer(shp_UAE, width=40, byid=TRUE, quadsegs=10)
shp_buff <- gBuffer(coord_ABU_DHABI_point, width=0.5)
shp_buff <- spTransform(shp_buff, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(shp_buff)
plot(shp_UAE, add=TRUE, lwd=1)

# save shp file for the crcular buffer
setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS/prova_shapes")
shapefile(shp_buff, "circular_buffer.shp", overwrite=TRUE)
dir <- "Z:/_SHARED_FOLDERS/Air Quality/Phase 2/admin_GIS/prova_shapes"
# reload and plot domain
shp_buff <- readOGR(dsn = dir, layer = "circular_buffer")



# load .tif file for the DUST MASK

setwd("Z:/_SHARED_FOLDERS/Air Quality/Phase 2/DUST SEVIRI/masks")
patt <- ".tif"
filenames <- list.files(pattern = patt)

filenames <- filenames[2]
mask <- raster(filenames)

# get points from the raster (lat, lon, points)
values <- rasterToPoints(mask)
colnames(values) <- c("x", "y", "values")
values <- as.data.frame(values) 
head(values)

crs <- projection(shp_buff) ### get projections from shp file
values <- SpatialPointsDataFrame(values[,1:2], values, 
                                       proj4string=CRS(crs))

plot(values)



pts.poly <- over(values, shp_buff[ ,"ID"])
NO2_2013_OMI$id <- pts.poly$GADMID

# Aggregate by zone/polygon
###  Make a dataframe ###
data_points <- NO2_2013_OMI@data
head(data_points)


#### Sum of emission in UK polygon [molecules/cm2]
data_points <- data_points %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(NO2_sum = sum(NO2_2013))


NO2_sum_UK <- (data_points$NO2_sum[1]) ### molecules/cm2
Surface_UK <- (shp@data$SQKM)*1e+10 ### from km2 to cm2
NO2_MASS <- 46.055 ### g/mol
N_Avogadro <- 6.022140857e+23
Total_mass_grams <- (NO2_sum_UK * Surface_UK * NO2_MASS)/ N_Avogadro
Total_mass_ktonnes <- Total_mass_grams/1000000000


