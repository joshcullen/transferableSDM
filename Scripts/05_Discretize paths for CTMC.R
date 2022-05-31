
### Create gridded paths for turt tracks ###
# Original code modified from Supplement of Brennan et al 2018 (https://static-content.springer.com/esm/art%3A10.1007%2Fs10980-018-0642-z/MediaObjects/10980_2018_642_MOESM3_ESM.html)

library(tidyverse)
library(lubridate)
library(sf)
library(crawlUtils)
library(furrr)
library(future)
# library(terra)
# library(ctmcmove)

source("Scripts/make_hex_grid2.R")


## Load turtle tracks
dat <- read.csv("Processed_data/Imputed_Cm_Tracks_SSM_30min.csv") %>%
  mutate(across(datetime, as_datetime))

glimpse(dat)

dat.sf.p <- dat %>%
  # filter(ptt == 181800) %>%
  st_as_sf(., coords = c('mu.x', 'mu.y'), crs = 3395)

cellsize <- hex_size(area = 4000^2)

plan("multisession", workers = availableCores() - 2)
hex_grid <- dat.sf.p %>%
  split(.$ptt) %>%
  future_map(~hex_grid_sfsample(sf_tracks = ., cellsize = cellsize, buffer = 15000)) %>%
  set_names(unique(dat.sf.p$ptt))

dat.sf.l <- dat.sf.p %>%
  group_by(ptt, rep, bout) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

plotly::ggplotly(
  ggplot() +
    geom_sf(data = hex_grid$`181807`$poly, fill = "transparent", size = 0.25) +
    geom_sf(data = dat.sf.l %>%
              filter(ptt == 181807), aes(color = rep), size = 0.25, alpha = 0.6) +
    geom_sf(data = dat.sf.p %>%
              filter(ptt == 181807), aes(color = rep), size = 0.25, alpha = 0.6) +
    theme_bw()
)


## Load environmental layers
sst <- rast("Environ_data/GoM SST example.tif")
sst.proj <- sst %>%
  project('EPSG:3395') %>%
  stack()  #MUST be of class RasterStack in order for path2ctmc() to work



## Create discrete path
dat.list <- dat %>%
  filter(ptt %in% c(181800, 181807)) %>%
  split(.$ptt)

n.turts <- length(dat.list)

ctmc.list <- list()
for (i in 1:n.turts){
  ctmc.list[[i]] <- path2ctmc(xy = as.matrix(dat.list[[i]][, c('mu.x','mu.y')]), t = as.numeric(dat.list[[i]][, "datetime"]),
                              method="ShortestPath", rast = sst.rast, directions = 4, print.iter = TRUE)
}



## Turn CTMC path into format for Poisson GLM

int <- sst.rast[[1]]
examplerast <- int
glm.list <- list()

n.list <- length(ctmc.list)

for (i in 1:n.list){
  glm.list[[i]] <- ctmc2glm(ctmc.list[[i]], stack.static = sst.rast, stack.grad = NULL, grad.point.decreasing = FALSE,
                            crw = TRUE)
}

for(i in 1:n.list){
  glm.list[[i]]$id <- rep(as.numeric(names(dat.list)[i]), times = length(glm.list[[i]]$z)) # To add individual elk ids to each list item
}

turt.dat <- bind_rows(glm.list)
# elkdata = glm.list[[1]]
# for(i in 2:n.list){
#   elkdata <- rbind(elkdata, glm.list[[i]]) # To bind all the data together in one data frame
# }

idx <- which(turt.dat$tau<10^-5 | turt.dat$tau>100) # Identifies instantaneous transitions and overly long stays
turt.dat <- turt.dat[-idx,] # Removes any

head(elkdata) # We do not provide code here for extracting NDVI values closest in date to GPS locations, but for the next steps we uploaded an elkdata file that contains this NDVI information. Also note: the dem.1 and forest.1 fields are gradient values.
