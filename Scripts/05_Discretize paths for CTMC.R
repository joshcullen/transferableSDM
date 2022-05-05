
### Create gridded paths for turt tracks ###
# Original code modified from Supplement of Brennan et al 2018 (https://static-content.springer.com/esm/art%3A10.1007%2Fs10980-018-0642-z/MediaObjects/10980_2018_642_MOESM3_ESM.html)

library(tidyverse)
library(lubridate)
library(terra)
library(ctmcmove)


## Load turtle tracks
dat <- read.csv("Processed_data/Processed_Cm_Tracks_SSM.csv") %>%
  mutate(across(date, as_datetime))


## Load environmental layers
sst <- read.csv("Environ_data/GoM SST example.csv")
sst.rast <- sst %>%
  split(.$time) %>%
  map({. %>%
      dplyr::select(-time)}) %>%
  map(~rasterFromXYZ(., crs = "EPSG:4326", res = c(0.01, 0.01), digits = 2)) %>%
  stack()  #MUST be of class RasterStack in order for path2ctmc() to work



## Create discrete path
dat.list <- dat %>%
  filter(id %in% c(181800, 181807)) %>%
  split(.$id)

n.turts <- length(dat.list)

ctmc.list <- list()
for (i in 1:n.turts){
  ctmc.list[[i]] <- path2ctmc(xy = as.matrix(dat.list[[i]][, c('lon','lat')]), t = as.numeric(dat.list[[i]][, "date"]), method="ShortestPath",
                              rast = sst.rast, directions = 4, print.iter = TRUE)
}



## Turn CTMC path into format for Poisson GLM

int <- stack_rast[[1]]
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
