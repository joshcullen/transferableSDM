
### Explore tSSF results ###

library(tidyverse)
library(lubridate)
library(rstan)
library(MCMCvis)
library(bayesplot)
library(terra)
library(sf)
library(furrr)
library(future)
library(tictoc)
library(progressr)

source('Scripts/helper functions.R')



### Load model fit, model input, and tracks ###

mod <- readRDS('Data_products/tSSF_model_GLMM_stanfit.rds')
tSSF.input <- read_csv("Processed_data/Input for tSSF.csv")
dat <- read_csv('Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foieGras.csv')

# Change time step from secs to mins
tSSF.input$dt <- tSSF.input$dt/60

gom <- sfarrow::st_read_parquet("Environ_data/GoM_land.parquet") %>%
  st_transform(3395)



#######################################
### Import Environmental Covariates ###
#######################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
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




##############################
### Extract environ covars ###
##############################

# Extract covar values at end of each used/available step
plan(multisession, workers = availableCores() - 2)
tSSF.input2 <- extract.covars(data = tSSF.input, layers = cov_list, dyn_names = c('k490', 'npp', 'sst'),
                              along = FALSE, ind = "month.year", imputed = FALSE)
#takes 1.5 min to run on desktop (18 cores)
plan(sequential)


# Center and scale covariates
tSSF.input3 <- tSSF.input2 %>%
  drop_na(bathym, k490, npp, sst) %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         k490.s = scale(k490) %>%
           as.vector(),
         npp.s = scale(npp) %>%
           as.vector(),
         sst.s = scale(sst) %>%
           as.vector())



scale.pars <- tSSF.input3 %>%
  dplyr::select(bathym, k490, npp, sst) %>%
  pivot_longer(cols = everything(), names_to = "covar", values_to = "value") %>%
  group_by(covar) %>%
  summarize(mean = mean(value),
            sd = sd(value)) %>%
  ungroup() %>%
  data.frame() %>%
  split(.$covar)


cov_list.s <- cov_list %>%
  map2(.x = ., .y = scale.pars, ~{(.x - .y$mean) / .y$sd})
