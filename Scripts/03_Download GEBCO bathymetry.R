
### Download GEBCO (2019) bathymetry data for each region ###

library(tidyverse)
library(terra)
library(sf)
library(cmocean)

source('Scripts/helper functions.R')


## Load processed tracks
turts <- read.csv("Processed_data/Processed_Cm_Tags_SSM.csv")




### Define extent for Gulf of Mexico tracks (including Bimini turts) ###

# define turtle track bbox
turts.gom <- turts %>%
  filter(lon < 0 & lat > 0) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
gom.bbox <- st_bbox(turts.gom)

# define {terra} SpatRaster based on bbox
bbox.gom <- terra::rast(xmin = gom.bbox[1], xmax = gom.bbox[3], ymin = gom.bbox[2], ymax = gom.bbox[4],
                        resolution = 15/3600, crs = "epsg:4326")
bbox.gom <- extend(bbox.gom, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)


## Download GEBCO 2019 bathymetry for GoM extent
gom.rast <- get_gebco(bbox.gom)


## Plot bathymetry data and tracks
gom.rast.df <- as.data.frame(gom.rast, xy = TRUE)
names(gom.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = gom.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("", name = "deep", direction = -1,
                     limits = c(min(gom.rast.df$depth), 0)) +
  geom_sf(data = turts.gom, aes(color = id)) +
  theme_bw() +
  coord_sf()






### Define extent for Brazil tracks ###

# define turtle track bbox
turts.brazil <- turts %>%
  filter(lon < 0 & lat < 0) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
brazil.bbox <- st_bbox(turts.brazil)

# define {terra} SpatRaster based on bbox
bbox.brazil <- terra::rast(xmin = brazil.bbox[1], xmax = brazil.bbox[3], ymin = brazil.bbox[2],
                           ymax = brazil.bbox[4], resolution = 15/3600, crs = "epsg:4326")
bbox.brazil <- extend(bbox.brazil, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)


## Download GEBCO 2019 bathymetry for Brazil extent
brazil.rast <- get_gebco(bbox.brazil)


## Plot bathymetry data and tracks
brazil.rast.df <- as.data.frame(brazil.rast, xy = TRUE)
names(brazil.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = brazil.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("", name = "deep", direction = -1,
                     limits = c(min(brazil.rast.df$depth), 0)) +
  geom_sf(data = turts.brazil, aes(color = id)) +
  theme_bw() +
  coord_sf()






### Define extent for Qatar tracks ###

# define turtle track bbox
turts.qatar <- turts %>%
  filter(lon > 0 & lat > 0) %>%
  st_as_sf(., coords = c("lon","lat"), crs = 4326)
qatar.bbox <- st_bbox(turts.qatar)

# define {terra} SpatRaster based on bbox
bbox.qatar <- terra::rast(xmin = qatar.bbox[1], xmax = qatar.bbox[3], ymin = qatar.bbox[2],
                           ymax = qatar.bbox[4], resolution = 15/3600, crs = "epsg:4326")
bbox.qatar <- extend(bbox.qatar, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)


## Download GEBCO 2019 bathymetry for Qatar extent
qatar.rast <- get_gebco(bbox.qatar)


## Plot bathymetry data and tracks
qatar.rast.df <- as.data.frame(qatar.rast, xy = TRUE)
names(qatar.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = qatar.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("", name = "deep", direction = -1,
                     limits = c(min(qatar.rast.df$depth), 0)) +
  geom_sf(data = turts.qatar, aes(color = id)) +
  theme_bw() +
  coord_sf()



### Export bathymetry rasters ###

terra::writeRaster(gom.rast, "Environ_data/GoM bathymetry.tif")
terra::writeRaster(brazil.rast, "Environ_data/Brazil bathymetry.tif")
terra::writeRaster(qatar.rast, "Environ_data/Qatar bathymetry.tif")
