
### Fit RSF as GLMM ###

library(tidyverse)
library(lubridate)
library(bayesmove)
library(terra)
# library(INLA)
# library(future)
# library(furrr)
library(sf)
library(sfarrow)
library(tictoc)
library(amt)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_4hr_aniMotum.csv")
gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(dat)
summary(dat)


dat <- dat %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))



###########################################################################
### Calculate step lengths, turning angles, and dt; filter observations ###
###########################################################################

dat <- prep_data(dat = dat, coord.names = c('x','y'), id = "id")
dat$speed <- dat$step / dat$dt


# Check speed for any outliers
round(quantile(dat$speed, c(0.01, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE), 2)
any(dat$step == 0, na.rm = TRUE)  #no step lengths equal to 0


# Remove any observations w/ speeds faster than 3 m/s and w/ large time gaps (i.e., 3 days)
dat2 <- dat %>%
  filter(speed <= 3) %>%
  filter(!is.na(bout))




#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'k490', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('k490', 'npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
  # terra::time(cov_list[[var]]) <- as_date(names(cov_list[[var]]))
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}

## Center and scale (by 1 SD) all covar layers
# cov_list_s <- sapply(cov_list, scale_across_covar)


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')
# cov_list_s <- map(cov_list_s, terra::project, 'EPSG:3395')




#####################################################################
### Generate random (available) points and extract environ covars ###
#####################################################################

# Convert to class for use of {amt} functions
dat3 <- make_track(dat2, .x = x, .y = y, .t = date, id = id, month.year = month.year, crs = 3395) %>%
  nest(data = -id)

# Remove points that intersect land
### PICKBACK UP FROM HERE

# Create KDE_href for each ID and add as nested column
dat4 <- dat3 %>%
  mutate(ud = map(data, ~hr_akde(.x, levels = 0.99)))

# Check mean step length
summary(dat2$step)

dat4 <- dat4 %>%
  mutate(ud_buff = map(ud, ~{.x %>%
      hr_isopleths() %>%
      slice(2) %>%
      sf::st_buffer(dist = 2e4)}  #add 2 km buffer (mean step length)
      ))

# Clip UDs by land mask
dat5 <- dat4 %>%
  mutate(ud_buff = map(ud_buff, ~st_mask(.x, gom.sf)))


# Generate available points
dat5 <- dat5 %>%
  mutate(avail_pts = map2(.x = ud_buff, .y = data, ~random_points(.x, n = 1000, presence = .y)))
