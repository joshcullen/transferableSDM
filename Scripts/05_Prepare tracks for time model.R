
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

############################
### Import turtle tracks ###
############################

dat <- read_csv('Processed_data/Processed_GoM_Cm_Tracks_SSM_12hr_aniMotum.csv')

dat <- dat %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))



# Calculate step lengths, turning angles, NSD, and dt
dat<- prep_data(dat = dat, coord.names = c('x','y'), id = "id")
dat$speed <- dat$step / dat$dt


# Check speed for any outliers
round(quantile(dat$speed, c(0.01, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE), 2)
any(dat$step == 0, na.rm = TRUE)
sum(dat$step == 0, na.rm = TRUE)  #none


# Create available steps per observed step
set.seed(2022)

plan(multisession, workers = availableCores() - 2)
dat2 <- add_avail_steps(dat, nsteps = 100)
plan(sequential)  #took 48 sec to run


# Remove observations w/ NA step length (i.e., last obs of each PTT)
dat.filt <- dat2 %>%
  filter(!is.na(step))


# Remove any observations w/ speeds faster than 3 m/s
dat.filt <- dat.filt %>%
  filter(speed <= 3)






#######################################
### Import Environmental Covariates ###
#######################################

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
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- resample(cov_list[[var]], cov_list$npp, method = "average")
}


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')


########################################################
### Extract values from raster layers for each track ###
########################################################

dat.filt <- dat.filt %>%
  filter(id != 104833)  #no data available for Kd490 and NPP in 2011

# Extract avg covar values per used step (for time model)
plan(multisession, workers = availableCores() - 2)
obs.path <- extract.covars(data = dat.filt %>%
                             filter(obs == 1), layers = cov_list, dyn_names = c('k490', 'npp', 'sst'),
                       along = TRUE, ind = "month.year", imputed = FALSE)
#takes 1.5 min to run on desktop (18 cores)
plan(sequential)


# add strata and obs columns
obs.path1 <- cbind(obs.path, dat.filt[dat.filt$obs == 1, c("strata","obs")])


# save(dat, dat.filt, cov_list, obs.path1, file = "Data_products/Extracted environ covars.RData")


nrow(drop_na(obs.path1, bathym, k490, npp, sst)) / nrow(obs.path1)
## 70.9% of observed and available steps have complete data




## Export data
write_csv(obs.path1, "Processed_data/Input for time model.csv")
write_csv(dat.filt, "Processed_data/Input for tSSF.csv")
# arrow::write_parquet(path1, "Processed_data/Input for time model.parquet")
