

### Download dyanmic environmental variables as raster layers from ERDDAP ###
### Define extents on which to download GEBCO (2021) bathymetry data for each region from website ###
# https://download.gebco.net

library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(cmocean)
library(rerddapXtracto)
library(tictoc)

source("Scripts/helper functions.R")


## Load processed tracks
turts <- read.csv("Processed_data/Processed_Cm_Tracks_SSM_30min.csv")


### Define extent for Gulf of Mexico tracks (including Bimini turts) ###

# define turtle track bbox
turts.gom <- turts %>%
  # filter(lon < 0 & lat > 0) %>%
  st_as_sf(., coords = c("mu.x","mu.y"), crs = 3395, remove = FALSE) %>%
  st_transform(4326)
gom.bbox <- st_bbox(turts.gom)

# define {terra} SpatRaster based on bbox
bbox.gom <- terra::rast(xmin = gom.bbox[1], xmax = gom.bbox[3], ymin = gom.bbox[2], ymax = gom.bbox[4],
                        resolution = 15/3600, crs = "EPSG:4326")
bbox.gom <- extend(bbox.gom, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)
ext(bbox.gom)




##################
### Bathymetry ###
##################

## Load GEBCO bathymetry for GoM extent
gom.rast <- rast("Environ_data/GoM bathymetry.tif")


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
ext(bbox.brazil)

## Load GEBCO bathymetry for Brazil extent
brazil.rast <- rast("Environ_data/Brazil bathymetry.tif")


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
ext(bbox.qatar)

## Load GEBCO bathymetry for Qatar extent
qatar.rast <- rast("Environ_data/Qatar bathymetry.tif")


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





###############################
### Sea surface temperature ###
###############################

## Plot SST basemap based on defined bbox for particular month-year

xpos <- ext(bbox.gom)[1:2]
ypos <- ext(bbox.gom)[3:4]
tpos <- as_date(range(turts.gom$date)) %>% as.character()
sstInfo <- rerddap::info('jplMURSST41mday')

tic()
sst.bbox <- rxtracto_3D(sstInfo, parameter = 'sst', xcoord = xpos, ycoord = ypos, tcoord = tpos)
toc()
# takes 6.5 min to run

# plotBBox(sst.bbox, plotColor = 'thermal')



# Create {terra} SpatRaster for export and data.frame to plot raster in ggplot
sst.rast <- array2rast(lon = sst.bbox$longitude, lat = sst.bbox$latitude, var = sst.bbox$sst,
                       time = sst.bbox$time, extent = ext(bbox.gom))

sst.rast.df <- as.data.frame(sst.rast[[1:12]], xy = TRUE) %>%  #select 1st 12 layers as example for viz
  pivot_longer(cols = -c("x","y"), names_to = "date", values_to = "sst") %>%
  arrange(date)


# Plot tracks overlaid w/ SST
turts.gom.l <- turts.gom %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_raster(data = sst.rast.df, aes(x, y, fill = sst)) +
  scale_fill_cmocean() +
  geom_sf(data = turts.gom.l, aes(color = factor(ptt))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_sf() +
  facet_wrap(~date)


### Export SST rasters as CSV (w/ WGS84 CRS; EPSG:4326)

# writeRaster(sst.rast, "Environ_data/GoM SST example.tif")
# writeRaster(sst.rast, "Environ_data/GoM SST.tif")











################################
### Net Primary Productivity ###
################################

## Plot NPP basemap based on defined bbox for particular month-year

xpos <- ext(bbox.gom)[1:2]
ypos <- ext(bbox.gom)[3:4]
tpos <- as_date(range(turts.gom$date)) %>% as.character()
nppInfo <- rerddap::info('erdMH1ppmday')
nppInfo$alldata$NC_GLOBAL[38,]

tic()
npp.bbox <- rxtracto_3D(nppInfo, parameter = 'productivity', xcoord = xpos, ycoord = ypos, tcoord = tpos,
                        zcoord = c(0,0))
toc()
# takes 30 sec to run


plotBBox(npp.bbox, plotColor = 'algae')


# Create {terra} SpatRaster for export and data.frame to plot raster in ggplot
npp.rast <- array2rast(lon = npp.bbox$longitude, lat = npp.bbox$latitude, var = npp.bbox$productivity,
                       time = npp.bbox$time, extent = ext(bbox.gom))

npp.rast.df <- as.data.frame(npp.rast[[1:12]], xy = TRUE) %>%  #select 1st 12 layers as example for viz
  pivot_longer(cols = -c("x","y"), names_to = "date", values_to = "npp") %>%
  arrange(date)


# Plot tracks overlaid w/ NPP
turts.gom.l <- turts.gom %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_raster(data = npp.rast.df, aes(x, y, fill = npp)) +
  scale_fill_cmocean(name = "algae") +
  geom_sf(data = turts.gom.l, aes(color = factor(ptt))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_sf() +
  facet_wrap(~date)


### Export NPP rasters as CSV (w/ WGS84 CRS; EPSG:4326)

# writeRaster(npp.rast, "Environ_data/GoM NPP.tif")
