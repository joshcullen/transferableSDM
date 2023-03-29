
### Fit HGPR RSF at different spatial scales ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(future)
library(furrr)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10_5km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_10_10km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_10km.csv")
rsf.pts_10_20km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_20km.csv")
rsf.pts_10_40km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_40km.csv")
rsf.pts_10_80km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_80km.csv")

rsf.list <- list(
  sc.5 = rsf.pts_10_5km,
  sc.10 = rsf.pts_10_10km,
  sc.20 = rsf.pts_10_20km,
  sc.40 = rsf.pts_10_40km,
  sc.80 = rsf.pts_10_80km
)

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10_5km)
summary(rsf.pts_10_5km)



#####################
### Fit HGPR RSFs ###
#####################

# Remove rows w/ incomplete observations; log-transform covars
rsf.list2 <- rsf.list %>%
  map(., ~{.x %>%
      drop_na(bathym, npp, sst) %>%
      mutate(log.bathym = log(abs(bathym)),
             log.npp = log(npp),
             log.sst = log(sst))}
      )


### Down-weighted Poisson regression

# Define pixel area for each spatial scale (in m^2)
Area.list <- list(sc.5 = 4759.836 ^ 2,
                  sc.10 = 9532.72 ^ 2,
                  sc.20 = 19065.44 ^ 2,
                  sc.40 = 38130.88 ^ 2,
                  sc.80 = 76261.76 ^ 2)
rsf.list2 <- rsf.list2 %>%
  map2(.x = ., .y = Area.list,
       ~{.x %>%
           mutate(wts = case_when(obs == 0 ~ .y / sum(obs == 0),
                                  obs == 1 ~ 1e-6))}
  )



# Add ID as integer for INLA
rsf.list2 <- rsf.list2 %>%
  map(.x = .,
       ~{.x %>%
           mutate(id1 = as.integer(factor(id)))}
  )

# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup <- n_distinct(rsf.list2$sc.5$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)

nbasis <- 5  #number of basis functions for approximating GP
degree <- 2  #degree for defining 1D mesh of GP
alpha <- 2  #for calculating Matern covariance matrix


plan(multisession, workers = availableCores() - 2)
hgpr.res <- future_map(rsf.list2, ~fit_hgpr(data = ., covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                                            nbasis = 5, degree = 2, alpha = 2),
                      .options = furrr_options(seed = TRUE))
plan(sequential)



# hgpr.mod_5km <- fit_hgpr(data = rsf.list2$sc.5, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
#                          nbasis = 5, degree = 2, alpha = 2)

summary(hgpr.mod_5km)
