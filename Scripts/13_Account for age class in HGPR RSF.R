
### Fit correlative HGPR RSF at 5 km scale and accounting for life stage ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(future)
library(furrr)
library(tidyterra)
library(patchwork)
library(ggh4x)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10_5km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv") %>%
  mutate(across(id, as.character))
dat <- read_csv("Processed_data/Prefiltered_Cm_Tracks.csv") %>%
  filter(Region == "GoM") %>%
  dplyr::select(Age, Ptt) %>%
  distinct()

# Join Age to RSF input
rsf.pts_10_5km <- rsf.pts_10_5km %>%
  left_join(dat, by = c("id" = "Ptt"))

glimpse(rsf.pts_10_5km)
summary(rsf.pts_10_5km)

# Load land spatial layer
gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")





################################
### Fit correlative HGPR RSF ###
################################

# Remove rows w/ incomplete observations; log-transform covars
rsf.pts_10_5kms <- rsf.pts_10_5km %>%
  drop_na(bathym, npp, sst) %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))



# Down-weighted Poisson regression
A <- 2345557  #area of study region (Gulf of Mexico) in km^2; from region used to generate 'available' pts
rsf.pts_10_5kms$wts <- ifelse(rsf.pts_10_5kms$obs == 0, A / sum(rsf.pts_10_5kms$obs == 0), 1e-6)




# Add Age and ID as integer for INLA
rsf.pts_10_5kms <- rsf.pts_10_5kms %>%
  mutate(Age1 = as.integer(factor(Age, levels = c("Juv","Adult"))),
         id1 = as.integer(factor(id)))



# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup.age <- n_distinct(rsf.pts_10_5kms$Age1)
ngroup.id <- n_distinct(rsf.pts_10_5kms$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)


hgpr.age <- fit_hgpr(data = rsf.pts_10_5kms, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                     nbasis = 5, degree = 2, alpha = 2, age.class = TRUE, int.strategy = 'auto',
                     method = "corr")
# took 1 min to run

summary(hgpr.age)





##############################################
### Age class-level marginal effects plots ###
##############################################

# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),
                                 degree = 2,
                                 boundary = 'free')
}


# Generate matrices for covariate raster data (for prediction)
A.me.age <- vector("list", length(covars))
newdat.list <- map(mesh.seq, function(x) {seq(from = x[1], to = x[2], length.out = 500)})
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.me.age[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.list[[covars[[i]]]])
}

# Replicate list elements for mapping over lists
A.me.age <- rep(A.me.age, each = max(rsf.pts_10_5kms$Age1))

# Store resulting GP coeffs per covar into a list
pred.coeffs.age <- hgpr.age$summary.random[1:3] %>%
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, age = rep(1:2, each = 6))) %>%
  map(., ~split(.x, .x$age)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10_5kms$Age1)), rep(1:max(rsf.pts_10_5kms$Age1), length(covars)), sep = "_"))

# Make predictions via linear algebra
pred.vals.age <- A.me.age %>%
  map2(.x = ., .y = pred.coeffs.age,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10_5kms$Age1)), rep(1:max(rsf.pts_10_5kms$Age1), length(covars)), sep = "_")) %>%
  bind_rows() %>%
  mutate(bathym = newdat.list$log.bathym,
         npp = newdat.list$log.npp,
         sst = newdat.list$log.sst) %>%
  mutate(across(everything(), exp)) %>%
  pivot_longer(cols = -c(bathym, npp, sst), names_to = "label", values_to = "mean") %>%
  separate(col = label, into = c('covar','age'), sep = "_")


# Wrangle data so that x and y values match up in long format
pred.vals.age <- pred.vals.age %>%
  mutate(x = case_when(str_detect(covar, "log.bathym") ~ bathym,
                       str_detect(covar, "log.npp") ~ npp,
                       str_detect(covar, "log.sst") ~ sst)) %>%
  dplyr::select(-c(bathym, npp, sst)) %>%
  mutate(covar = gsub(pattern = "log.", replacement = "", x = covar)) %>%
  mutate(across(age, factor))
levels(pred.vals.age$age) <- c("Juv","Adult")



# Facet of covars by age class
ggplot() +
  geom_line(data = pred.vals.age, aes(x = x, y = mean, color = age), linewidth = 1.5) +
  theme_bw() +
  labs(y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(age ~ covar, scales = "free")


# Depth
ggplot() +
  geom_line(data = pred.vals.age %>%
              filter(covar == 'bathym'), aes(x = x, y = mean, color = age), linewidth = 1) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ age, scales = "free")


# NPP
ggplot() +
  geom_line(data = pred.vals.age %>%
              filter(covar == 'npp'), aes(x = x / 1000, y = mean, color = age), linewidth = 1) +
  theme_bw() +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ age, scales = "free")


# SST
ggplot() +
  geom_line(data = pred.vals.age %>%
              filter(covar == 'sst'), aes(x = x, y = mean, color = age), linewidth = 1) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ age, scales = "free")






#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
names(cov_list) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list[["bathym"]][cov_list[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')








####################################
### Generate predictive surfaces ###
####################################

rast.sep.20 <- data.frame(log.bathym = log(abs(terra::values(cov_list$bathym))) %>%
                            as.vector(),
                          log.npp = log(terra::values(cov_list$npp$`2020-09-01`)) %>%
                            as.vector(),
                          log.sst = log(terra::values(cov_list$sst$`2020-09-01`)) %>%
                            as.vector()) %>%
  mutate(row_id = 1:nrow(.)) %>%
  drop_na(log.bathym, log.npp, log.sst)



# Generate matrices for covariate raster data (for prediction)
A.age.sep.20 <- vector("list", length(covars))
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.age.sep.20[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = rast.sep.20[[covars[[i]]]])
}


# Replicate list elements for mapping over lists
A.me.age <- rep(A.me.age, each = max(rsf.pts_10_5kms$Age1))

# Store resulting GP coeffs per covar into a list
pred.coeffs.age <- hgpr.age$summary.random[1:3] %>%
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, age = rep(1:2, each = 6))) %>%
  map(., ~split(.x, .x$age)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10_5kms$Age1)), rep(1:max(rsf.pts_10_5kms$Age1), length(covars)), sep = "_"))

# Make predictions via linear algebra
pred.vals.juv <- A.age.sep.20 %>%
  map2(.x = ., .y = pred.coeffs.age[grep("1", names(pred.coeffs.age))],
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  bind_cols() %>%
  rowSums()  #sum up all predictions across covars

pred.vals.adult <- A.age.sep.20 %>%
  map2(.x = ., .y = pred.coeffs.age[grep("2", names(pred.coeffs.age))],
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  bind_cols() %>%
  rowSums()  #sum up all predictions across covars


# Define dummy raster for storing predictions
rast.pred.juv <- rast.pred.adult <- cov_list$bathym
terra::values(rast.pred.juv) <- terra::values(rast.pred.adult) <- NA  # initially store all NAs for locs w/o predictions
terra::values(rast.pred.juv)[rast.sep.20$row_id] <- pred.vals.juv
terra::values(rast.pred.adult)[rast.sep.20$row_id] <- pred.vals.adult

# Define spatial extent
bbox <- ext(rast.pred.juv)


# Generate predictive surface for GP at pop-level
p.juv <- ggplot() +
  geom_spatraster(data = rast.pred.juv) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "Juv Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)

p.adult <- ggplot() +
  geom_spatraster(data = rast.pred.adult) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "Adult Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)


p.juv + p.adult






###########################
### Export model object ###
###########################

saveRDS(hgpr.age, "Data_products/HGPR_corr_model_fit_scale_age.rds")


