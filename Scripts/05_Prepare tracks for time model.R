
### Extract environmental covariates from imputed green turtle tracks

library(bayesmove)
library(tidyverse)
library(sf)
library(terra)
library(lubridate)
library(furrr)
library(future)
library(vroom)
library(tictoc)

source('Scripts/helper functions.R')

####################################
### Import imputed turtle tracks ###
####################################

dat <- vroom('Processed_data/Imputed_Cm_Tracks_SSM_2hr.csv', delim = ",")
dat <- dat %>%
  mutate(month.year = as_date(datetime),
         .after = 'datetime') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  rename(x = mu.x, y = mu.y, date = datetime, id = ptt)



# Calculate step lengths, turning angles, NSD, and dt
tic()
dat<- prep_data(dat = dat, coord.names = c('x','y'), id = "rep")
toc()  # takes 5.3 min to run

dat$speed <- dat$step / dat$dt


# Check speed for any outliers
round(quantile(dat$speed, c(0.01, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE), 2)
any(dat$step == 0, na.rm = TRUE)


## There are some extremely fast outliers (9146 m/s) that need to be removed before subsequent analysis

# Remove observations w/ NA step length (i.e., last obs of each PTT)
dat.filt <- dat %>%
  filter(!is.na(step))


# Remove any observations w/ speeds faster than 3 m/s (approx 99th quantile and filter used in CTCRW model)
dat.filt <- dat.filt %>%
  filter(speed <= 3)





#######################################
### Import Environmental Covariates ###
#######################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'Chla', 'Kd490', 'SST')

# Change names for NPP and SST to match KdPAR (YYYY-MM-01)
for (var in c('Chla', 'Kd490', 'SST')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}


## Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "SST")) {
  cov_list[[var]] <- resample(cov_list[[var]], cov_list$Chla, method = "bilinear")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')


########################################################
### Extract values from raster layers for each track ###
########################################################

dat.filt <- dat.filt %>%
  rename(id = ptt) %>%
  filter(id != 104833)

#for running in parallel; define number of cores to use
plan(multisession, workers = availableCores() - 2)
path<- extract.covars(data = dat.filt, layers = cov_list, dyn_names = c('Chla','Kd490','SST'),
                      ind = "month.year", imputed = TRUE)
#takes 2.8 hrs to run on desktop


plan(sequential)


save(dat, dat.filt, cov_list, path, file = "Data_products/Extracted environ covars.RData")


## 75% of observed steps have complete data

