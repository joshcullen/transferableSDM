

### Download dynamic environmental variables as raster layers from ERDDAP ###
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
gom.turts <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_4hr_aniMotum.csv")
br.turts <- read_csv("Processed_data/Processed_Brazil_Cm_Tracks_SSM_4hr_aniMotum.csv")
qa.turts <- read_csv("Processed_data/Processed_Qatar_Cm_Tracks_SSM_4hr_aniMotum.csv")

glimpse(gom.turts)

all.turts <- list(gom = gom.turts,
                  br = br.turts,
                  qa = qa.turts)

# Remove any observations before 2012-01-02; first available data for VIIRS Kd(PAR)
# all.turts <- all.turts %>%
#   map(., ~{.x %>%
#       filter(date > "2012-01-01")}  #results only in removal of PTT 104833 from GoM
#       )


### Define spatial extent for each region ###

# define turtle track bbox
all.turts.sf <- all.turts %>%
  map(., ~{.x %>%
      drop_na(x, y) %>%
      st_as_sf(coords = c("x","y"), crs = 3395, remove = FALSE) %>%
      st_transform(4326)}
      )

bbox <- all.turts.sf %>%
  map(st_bbox)

# define {terra} SpatRaster based on bbox
bbox2 <- bbox %>%
  map(., ~{terra::rast(xmin = .x$xmin, xmax = .x$xmax, ymin = .x$ymin, ymax = .x$ymax,
                  resolution = 15/3600, crs = "EPSG:4326") %>%
      extend(50)})  #pad each side by 5 cells (w/ 15 arc sec resolution)
ext(bbox2$gom)  #check resulting extent




##################
### Bathymetry ###
##################

## Load GEBCO bathymetry for GoM extent
gom.rast <- rast("Environ_data/GoM bathymetry.tif")


## Plot bathymetry data and tracks

turts.gom.l <- all.turts.sf %>%
  pluck("gom") %>%
  group_by(id) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

gom.rast.df <- as.data.frame(gom.rast, xy = TRUE)
names(gom.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = gom.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("Depth (m)", name = "deep", direction = -1,
                     limits = c(min(gom.rast.df$depth), 0)) +
  geom_sf(data = turts.gom.l, aes(color = factor(id))) +
  guides(size = "legend", color = "none") +
  theme_bw() +
  coord_sf(expand = FALSE)






### Define extent for Brazil tracks ###

# define turtle track bbox
# turts.brazil <- turts %>%
#   filter(lon < 0 & lat < 0) %>%
#   st_as_sf(., coords = c("lon","lat"), crs = 4326)
# brazil.bbox <- st_bbox(turts.brazil)
turts.brazil.l <- all.turts.sf %>%
  pluck("br") %>%
  group_by(id) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

# define {terra} SpatRaster based on bbox
# bbox.brazil <- terra::rast(xmin = brazil.bbox[1], xmax = brazil.bbox[3], ymin = brazil.bbox[2],
#                            ymax = brazil.bbox[4], resolution = 15/3600, crs = "epsg:4326")
# bbox.brazil <- extend(bbox.brazil, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)
# ext(bbox.brazil)

## Load GEBCO bathymetry for Brazil extent
brazil.rast <- rast("Environ_data/Brazil bathymetry.tif")


## Plot bathymetry data and tracks
brazil.rast.df <- as.data.frame(brazil.rast, xy = TRUE)
names(brazil.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = brazil.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("", name = "deep", direction = -1,
                     limits = c(min(brazil.rast.df$depth), 0)) +
  geom_sf(data = turts.brazil.l, aes(color = factor(id))) +
  guides(size = "legend", color = "none") +
  theme_bw() +
  coord_sf(expand = FALSE)






### Define extent for Qatar tracks ###

# define turtle track bbox
# turts.qatar <- turts %>%
#   filter(lon > 0 & lat > 0) %>%
#   st_as_sf(., coords = c("lon","lat"), crs = 4326)
# qatar.bbox <- st_bbox(turts.qatar)
turts.qatar.l <- all.turts.sf %>%
  pluck("qa") %>%
  group_by(id) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

# define {terra} SpatRaster based on bbox
# bbox.qatar <- terra::rast(xmin = qatar.bbox[1], xmax = qatar.bbox[3], ymin = qatar.bbox[2],
#                            ymax = qatar.bbox[4], resolution = 15/3600, crs = "epsg:4326")
# bbox.qatar <- extend(bbox.qatar, 50)  #pad each side by 5 cells (w/ 15 arc sec resolution)
# ext(bbox.qatar)

## Load GEBCO bathymetry for Qatar extent
qatar.rast <- rast("Environ_data/Qatar bathymetry.tif")


## Plot bathymetry data and tracks
qatar.rast.df <- as.data.frame(qatar.rast, xy = TRUE)
names(qatar.rast.df)[3] <- "depth"

ggplot() +
  geom_raster(data = qatar.rast.df, aes(x, y, fill = depth)) +
  scale_fill_cmocean("", name = "deep", direction = -1,
                     limits = c(min(qatar.rast.df$depth), 0)) +
  geom_sf(data = turts.qatar.l, aes(color = factor(id))) +
  guides(size = "legend", color = "none") +
  theme_bw() +
  coord_sf(expand = FALSE)





###############################
### Sea surface temperature ###
###############################

## Plot SST basemap based on defined bbox for particular month-year

rerddap.params <- bbox2 %>%
  map2(.x = ., .y = all.turts,
       ~{data.frame(xpos = ext(.x)[1:2],
                    ypos = ext(.x)[3:4],
                    tpos = range(.y$date) %>%
                      as_date() %>%
                      as.character())}
       )

# xpos <- ext(bbox.gom)[1:2]
# ypos <- ext(bbox.gom)[3:4]
# tpos <- range(turts.gom$datetime) %>%
#   as_date() %>%
#   as.character()
sstInfo <- rerddap::info('jplMURSST41mday')

sst.bbox.list <- vector("list", length(rerddap.params))
tic()
sst.bbox.list <- rerddap.params %>%
  map(., ~rxtracto_3D(sstInfo, parameter = 'sst', xcoord = .x$xpos, ycoord = .x$ypos, tcoord = .x$tpos))
toc()
# takes 12.25 min to run

# plotBBox(sst.bbox, plotColor = 'thermal')



# Create {terra} SpatRaster for export and data.frame to plot raster in ggplot
sst.rast.list <- vector("list", length(rerddap.params))
sst.rast.list <- sst.bbox.list %>%
  map2(.x = ., .y = bbox2, ~array2rast(lon = .x$longitude, lat = .x$latitude, var = .x$sst,
                       time = .x$time, extent = ext(.y))
  )

#select 1st 3 layers per Region as example for viz
tic()
sst.rast.sub <- sst.rast.list %>%
  map(., ~{as.data.frame(.x[[1:3]], xy = TRUE, na.rm = FALSE) %>%
  pivot_longer(cols = -c("x","y"), names_to = "date", values_to = "sst") %>%
  arrange(date)}
  ) %>%
  bind_rows(.id = "Region")
toc()  #took 30 sec


# Plot tracks overlaid w/ SST
# turts.gom.l <- turts.gom %>%
#   group_by(ptt) %>%
#   summarize(do_union = FALSE) %>%
#   st_cast("MULTILINESTRING")

ggplot() +
  geom_raster(data = sst.rast.sub %>%
                filter(Region == "gom"), aes(x, y, fill = sst)) +
  scale_fill_cmocean() +
  geom_sf(data = turts.gom.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)

ggplot() +
  geom_raster(data = sst.rast.sub %>%
                filter(Region == "br"), aes(x, y, fill = sst)) +
  scale_fill_cmocean() +
  geom_sf(data = turts.brazil.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)

ggplot() +
  geom_raster(data = sst.rast.sub %>%
                filter(Region == "qa"), aes(x, y, fill = sst)) +
  scale_fill_cmocean() +
  geom_sf(data = turts.qatar.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)


### Export SST rasters as GeoTIFF (w/ WGS84 CRS; EPSG:4326)

# writeRaster(sst.rast, "Environ_data/GoM SST example.tif")
# writeRaster(sst.rast.list$gom, "Environ_data/GoM SST.tif")
# writeRaster(sst.rast.list$br, "Environ_data/Brazil SST.tif")
# writeRaster(sst.rast.list$qa, "Environ_data/Qatar SST.tif")











################################
### Net primary productivity ###
################################

## Plot NPP basemap based on defined bbox for particular month-year

nppInfo <- rerddap::info('erdMH1ppmday')
nppInfo$alldata$NC_GLOBAL[38,]  #actually 0.0415 deg

npp.bbox.list <- vector("list", length(rerddap.params))
tic()
npp.bbox.list <- rerddap.params %>%
  map(., ~rxtracto_3D(nppInfo, parameter = 'productivity', xcoord = .x$xpos, ycoord = .x$ypos, tcoord = .x$tpos, zcoord = 0))
toc()
# takes 1.5 min to run


# Check that full time-series is available for NPP by Region (ie compare against SST)
map(sst.bbox.list, ~{.x[["time"]] %>% length()})
map(npp.bbox.list, ~{.x[["time"]] %>% length()})
# missing 1 NPP month from GoM
# missing 3 NPP month from Brazil
# missing 0 NPP month from Qatar


# Which months are missing per Region?
gom.ind <- which(!sst.bbox.list$gom$time %in% npp.bbox.list$gom$time)
sst.bbox.list$gom$time[gom.ind]  # August 2020 missing from GoM

br.ind <- which(!sst.bbox.list$br$time %in% npp.bbox.list$br$time)
sst.bbox.list$br$time[br.ind]  # August 2020, March 2022, April 2022, and July 2022 missing from br
# Apparently no Aug 2022 SST data, but doesn't matter since last transmission is July 29, 2022



# Download 8-day composite for missing monthly data
nppInfo.8d <- rerddap::info('erdMH1pp8day')
nppInfo.8d$alldata$NC_GLOBAL[38,]  #actually 0.0415 deg

rerddap.params.8d <- bbox2[1:2] %>%
  map(.,
       ~{data.frame(xpos = ext(.x)[1:2],
                    ypos = ext(.x)[3:4])}
  )
rerddap.params.8d$gom$tpos <- c("2020-08-01",'2020-08-31')
rerddap.params.8d$br$tpos <- c("2020-08-01",'2020-08-31')
rerddap.params.8d$br$tpos2 <- c("2022-03-01",'2022-07-31')  #I'll need to only retain March, April, and July 2022 data

npp.bbox.8d <- vector("list", 2)  #for storing Aug 2020 data for GoM and Brazil
tic()
npp.bbox.8d <- rerddap.params.8d %>%
  map(., ~rxtracto_3D(nppInfo.8d, parameter = 'productivity', xcoord = .x$xpos, ycoord = .x$ypos, tcoord = .x$tpos, zcoord = 0))
toc()
# takes 23 sec to run

npp.bbox.8d23 <- vector("list", 1)  #for storing 2023 data for Brazil
tic()
npp.bbox.8d23 <- rerddap.params.8d[2] %>%
  map(., ~rxtracto_3D(nppInfo.8d, parameter = 'productivity', xcoord = .x$xpos, ycoord = .x$ypos, tcoord = .x$tpos2, zcoord = 0))
toc()
# takes 1.25 min to run



# Create {terra} SpatRaster for export and merge mean NPP layers for missing months w/ rest of dataset
npp.rast.list <- vector("list", length(rerddap.params))
npp.rast.list <- npp.bbox.list %>%
  map2(.x = ., .y = bbox2, ~array2rast(lon = .x$longitude, lat = .x$latitude, var = .x$productivity,
                                       time = .x$time, extent = ext(.y))
  )

# Aug 2020
npp.rast.8d <- vector("list", length(rerddap.params.8d))
npp.rast.8d <- npp.bbox.8d %>%
  map2(.x = ., .y = bbox2[1:2], ~array2rast(lon = .x$longitude, lat = .x$latitude, var = .x$productivity,
                                       time = .x$time, extent = ext(.y))
  ) %>%
  map(., mean, na.rm = TRUE)
names(npp.rast.8d$gom) <- "2020-08-16"
names(npp.rast.8d$br) <- "2020-08-16"

# Brazil 2023
npp.rast.8d23 <- vector("list", 1)
npp.rast.8d23 <- npp.bbox.8d23 %>%
  map2(.x = ., .y = bbox2[2], ~array2rast(lon = .x$longitude, lat = .x$latitude, var = .x$productivity,
                                            time = .x$time, extent = ext(.y))
  )
ind.23 <- as_date(names(npp.rast.8d23$br)) %>%
                  month()
month.val <- c(3,4,7)
npp.rast.8d23.sub <- vector("list", 3)
for (i in 1:length(month.val)) {
  npp.rast.8d23.sub[[i]] <- npp.rast.8d23$br[[which(ind.23 %in% month.val[i])]] %>%
    mean(na.rm = TRUE)
}
npp.rast.8d23.sub <- rast(npp.rast.8d23.sub)
names(npp.rast.8d23.sub) <- c("2022-03-16","2022-04-16","2022-07-16")




# Merge back in missing monthly data
npp.rast.list$gom <- c(npp.rast.list$gom, npp.rast.8d$gom)
npp.rast.list$br <- c(npp.rast.list$br, npp.rast.8d$br, npp.rast.8d23.sub)

# Modify so as in chronological order
npp.rast.list$gom <- npp.rast.list$gom[[c(1:106, 109, 107:108)]]
npp.rast.list$br <- npp.rast.list$br[[c(1:51, 73, 52:69, 74:75, 70:71, 76)]]




#select 1st 3 layers per Region as example for viz
tic()
npp.rast.sub <- npp.rast.list %>%
  map(., ~{as.data.frame(.x[[1:3]], xy = TRUE, na.rm = FALSE) %>%
      pivot_longer(cols = -c("x","y"), names_to = "date", values_to = "npp") %>%
      arrange(date)}
  ) %>%
  bind_rows(.id = "Region")
toc()  #took 1.5 sec



ggplot() +
  geom_raster(data = npp.rast.sub %>%
                filter(Region == "gom"), aes(x, y, fill = npp)) +
  scale_fill_cmocean(name = "algae") +
  geom_sf(data = turts.gom.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)

ggplot() +
  geom_raster(data = npp.rast.sub %>%
                filter(Region == "br"), aes(x, y, fill = npp)) +
  scale_fill_cmocean(name = "algae") +
  geom_sf(data = turts.brazil.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)

ggplot() +
  geom_raster(data = npp.rast.sub %>%
                filter(Region == "qa"), aes(x, y, fill = npp)) +
  scale_fill_cmocean(name = "algae") +
  geom_sf(data = turts.qatar.l, aes(color = factor(id))) +
  scale_color_viridis_d(guide = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ date)





### Export NPP rasters as GeoTIFF (w/ WGS84 CRS; EPSG:4326)

# writeRaster(npp.rast.list$gom, "Environ_data/GoM NPP.tif", overwrite = TRUE)
# writeRaster(npp.rast.list$br, "Environ_data/Brazil NPP.tif", overwrite = TRUE)
# writeRaster(npp.rast.list$qa, "Environ_data/Qatar NPP.tif", overwrite = TRUE)





#
#
#
#
#
#
# ###########################################################
# ### Diffuse Attenuation Coefficient at 490 nm [Kd(490)] ###
# ###########################################################
#
# ## Plot Kd490 basemap based on defined bbox for particular month-year
#
# xpos <- ext(bbox.gom)[1:2]
# ypos <- ext(bbox.gom)[3:4]
# tpos <- range(turts.gom$datetime) %>%
#   as_date() %>%
#   as.character()
# # tpos[1] <- "2014-06-01"  #change to June 1st so that June 2014 layer is also downloaded
# kd490Info <- rerddap::info('erdVH2018k490mday')
# kd490Info$alldata$NC_GLOBAL[41,]
#
# tic()
# kd490.bbox <- rxtracto_3D(kd490Info, parameter = 'k490', xcoord = xpos, ycoord = ypos, tcoord = tpos)
# toc()
# # takes 13 sec to run
#
#
# # plotBBox(kd490.bbox, plotColor = 'turbid')
#
#
# # Create {terra} SpatRaster for export and data.frame to plot raster in ggplot
# kd490.rast <- array2rast(lon = kd490.bbox$longitude, lat = kd490.bbox$latitude, var = kd490.bbox$k490,
#                        time = kd490.bbox$time, extent = ext(bbox.gom))
# # names(kd490.rast) <- gsub(names(kd490.rast), pattern = " 12:00:00", replacement = "")
#
# kd490.rast.df <- as.data.frame(kd490.rast[[1:12]], xy = TRUE, na.rm = FALSE) %>%  #select 1st 12 layers as example for viz
#   pivot_longer(cols = -c("x","y"), names_to = "date", values_to = "kd490") %>%
#   arrange(date)
#
#
# # Plot tracks overlaid w/ NPP
# turts.gom.l <- turts.gom %>%
#   group_by(ptt) %>%
#   summarize(do_union = FALSE) %>%
#   st_cast("MULTILINESTRING")
#
# ggplot() +
#   geom_raster(data = kd490.rast.df, aes(x, y, fill = kd490)) +
#   scale_fill_cmocean(name = "turbid") +
#   geom_sf(data = turts.gom.l, aes(color = factor(ptt))) +
#   scale_color_viridis_d(guide = "none") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   coord_sf() +
#   facet_wrap(~date)


### Export Kd(490) rasters as GeoTIFF (w/ WGS84 CRS; EPSG:4326)

# writeRaster(kd490.rast, "Environ_data/GoM Kd490.tif")

