
### Compare model transferability by method ###

library(tidyverse)
library(INLA)
# library(mgcv)
# library(gratia)
# library(gbm)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(MetBrewer)
library(tidyterra)
library(patchwork)

source('Scripts/helper functions.R')


##########################
### Load fitted models ###
##########################

hglm_corr.fit <- readRDS("Data_products/HGLM_corr_model_fit.rds")
hglm_hybrid.fit <- readRDS("Data_products/HGLM_hybrid_model_fit.rds")
hgpr_corr.fit <- readRDS("Data_products/HGPR_corr_model_fit.rds")
hgpr_hybrid.fit <- readRDS("Data_products/HGPR_hybrid_model_fit.rds")




################################
### Load validation datasets ###
################################

dat.br <- read_csv("Processed_data/Brazil_Cm_Tracks_behav.csv")
dat.qa <- read_csv("Processed_data/Qatar_Cm_Tracks_behav.csv")

# Create indexing column "month.year" and only retain Resident locs
dat.br <- dat.br %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  filter(behav == 'Resident')

dat.qa <- dat.qa %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  filter(behav == 'Resident')


# Load spatial land layers
br.sf <- st_read_parquet("Environ_data/Brazil_land.parquet")
qa.sf <- st_read_parquet("Environ_data/Qatar_land.parquet")




########################################
### Load environmental raster layers ###
########################################

### Brazil ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Brazil", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_br <- sapply(files, rast)
names(cov_list_br) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_br[[var]]) <- gsub(names(cov_list_br[[var]]), pattern = "-..$", replacement = "-01")
}

# Set all positive bathymetric values (i.e., elevation) as NA
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > 0] <- NA

# Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list_br[[var]] <- terra::resample(cov_list_br[[var]], cov_list_br$npp, method = "average")
}

# Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > -1e-9] <- NA


# Transform CRS to match tracks
cov_list_br <- map(cov_list_br, terra::project, 'EPSG:3395')






### Qatar ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Qatar", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_qa <- sapply(files, rast)
names(cov_list_qa) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_qa[[var]]) <- gsub(names(cov_list_qa[[var]]), pattern = "-..$", replacement = "-01")
}

# Set all positive bathymetric values (i.e., elevation) as NA
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > 0] <- NA

# Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list_qa[[var]] <- terra::resample(cov_list_qa[[var]], cov_list_qa$npp, method = "average")
}

# Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > -1e-9] <- NA


# Transform CRS to match tracks
cov_list_qa <- map(cov_list_qa, terra::project, 'EPSG:3395')







#################################
### Validate correlative HGLM ###
#################################

### Brazil ###

my.ind.br <- names(cov_list_br$npp)
br.rast.tmp <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
br.rast.hglm <- list(corr = br.rast.tmp, hybrid = br.rast.tmp)  #store results for both corr and hybrid

names(br.rast.hglm$corr) <- my.ind.br
names(br.rast.hglm$hybrid) <- my.ind.br

# Define coeff values from HGLM
coeff_corr <- hglm_corr.fit$summary.fixed$mean
coeff_hybrid <- hglm_hybrid.fit$summary.fixed$mean

# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  vars <- data.frame(intercept = 1,
                     log.bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log() %>%
                       sapply(., function(x) x^2),
                     log.npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log(),
                     log.npp2 = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log() %>%
                       sapply(.,function(x) x^2),
                     log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log(),
                     log.sst2 = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log() %>%
                       sapply(., function(x) x^2))


  # Make predictions on intensity of use from models
  br.hglm_corr <- as.matrix(vars) %*% coeff_corr
  br.hglm_hybrid <- as.matrix(vars) %*% coeff_hybrid

  # Store results in raster stack
  terra::values(br.rast.hglm$corr[[i]]) <- br.hglm_corr  #keep on log-scale
  terra::values(br.rast.hglm$hybrid[[i]]) <- br.hglm_hybrid

}
skrrrahh('khaled2')
toc()  #took 45 sec


# Normalize predictions on 0-1 scale
br.rast.hglm2 <- br.rast.hglm |>
  map(normalize)


## Assess model performance via Continuous Boyce Index ##

## Correlative HGLM
boyce.br.full.hglm_corr <- vector("list", nlyr(br.rast.hglm2$corr))
boyce.br.sub.hglm_corr <- vector("list", nlyr(br.rast.hglm2$corr))
tic()
for (i in 1:nlyr(br.rast.hglm2$corr)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_main <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hglm_corr[[i]] <- boyce(fit = br.rast.hglm2$corr[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.hglm_corr[[i]] <- boyce(fit = br.rast.hglm2$corr[[i]],
                                  obs = obs_main,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 2.5 sec



# perc.use.br.full.hglm_corr <- boyce.br.full.hglm_corr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.full.hglm_corr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.6 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.full.hglm %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()
#
#
#
# perc.use.br.sub.hglm <- boyce.br.sub.hglm %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.sub.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #1.3 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.sub.hglm %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   ylim(0,1) +
#   theme_bw()


boyce.br.full.hglm_corr <- boyce.br.full.hglm_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGLM_corr")

boyce.br.sub.hglm_corr <- boyce.br.sub.hglm_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_main",
             Method = "HGLM_corr")





## Hybrid HGLM
boyce.br.full.hglm_hybrid <- vector("list", nlyr(br.rast.hglm2$hybrid))
boyce.br.sub.hglm_hybrid <- vector("list", nlyr(br.rast.hglm2$hybrid))
tic()
for (i in 1:nlyr(br.rast.hglm2$hybrid)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_main <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hglm_hybrid[[i]] <- boyce(fit = br.rast.hglm2$hybrid[[i]],
                                        obs = obs_full,
                                        nbins = 10,
                                        bin.method = "seq",
                                        PEplot = FALSE,
                                        rm.duplicate = FALSE,
                                        method = "spearman")

  boyce.br.sub.hglm_hybrid[[i]] <- boyce(fit = br.rast.hglm2$hybrid[[i]],
                                       obs = obs_main,
                                       nbins = 10,
                                       bin.method = "seq",
                                       PEplot = FALSE,
                                       rm.duplicate = FALSE,
                                       method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 2.5 sec



# perc.use.br.full.hglm_hybrid <- boyce.br.full.hglm_hybrid %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
# apply(perc.use.br.full.hglm_hybrid, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.6 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.full.hglm %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()
#
#
#
# perc.use.br.sub.hglm <- boyce.br.sub.hglm %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.sub.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #1.3 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.sub.hglm %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   ylim(0,1) +
#   theme_bw()


boyce.br.full.hglm_hybrid <- boyce.br.full.hglm_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGLM_hybrid")

boyce.br.sub.hglm_hybrid <- boyce.br.sub.hglm_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_main",
             Method = "HGLM_hybrid")





### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.tmp <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
qa.rast.hglm <- list(corr = qa.rast.tmp, hybrid = qa.rast.tmp)  #store results for both corr and hybrid

names(qa.rast.hglm$corr) <- my.ind.qa
names(qa.rast.hglm$hybrid) <- my.ind.qa

# Define coeff values from HGLM
coeff_corr <- hglm_corr.fit$summary.fixed$mean
coeff_hybrid <- hglm_hybrid.fit$summary.fixed$mean

# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  vars <- data.frame(intercept = 1,
                     log.bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log() %>%
                       sapply(., function(x) x^2),
                     log.npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log(),
                     log.npp2 = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log() %>%
                       sapply(.,function(x) x^2),
                     log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log(),
                     log.sst2 = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log() %>%
                       sapply(., function(x) x^2))


  # Make predictions on intensity of use from models
  qa.hglm_corr <- as.matrix(vars) %*% coeff_corr
  qa.hglm_hybrid <- as.matrix(vars) %*% coeff_hybrid

  # Store results in raster stack
  terra::values(qa.rast.hglm$corr[[i]]) <- qa.hglm_corr  #keep on log-scale
  terra::values(qa.rast.hglm$hybrid[[i]]) <- qa.hglm_hybrid

}
skrrrahh('khaled2')
toc()  #took 1 sec


# Normalize predictions on 0-1 scale
qa.rast.hglm2 <- qa.rast.hglm |>
  map(normalize)


## Assess model performance via Continuous Boyce Index ##

## Correlative HGLM
boyce.qa.hglm_corr <- vector("list", nlyr(qa.rast.hglm2$corr))
tic()
for (i in 1:nlyr(qa.rast.hglm2$corr)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hglm_corr[[i]] <- boyce(fit = qa.rast.hglm2$corr[[i]],
                              obs = obs,
                              nbins = 10,
                              bin.method = "seq",
                              PEplot = FALSE,
                              rm.duplicate = FALSE,
                              method = "spearman")
  }
skrrrahh("khaled3")
toc()  #took 1 sec


# perc.use.qa.hglm <- boyce.qa.hglm %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.qa.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.6 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.qa.hglm %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()



## Hybrid HGLM
boyce.qa.hglm_hybrid <- vector("list", nlyr(qa.rast.hglm2$hybrid))
tic()
for (i in 1:nlyr(qa.rast.hglm2$hybrid)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hglm_hybrid[[i]] <- boyce(fit = qa.rast.hglm2$hybrid[[i]],
                                   obs = obs,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec




boyce.qa.hglm_corr <- boyce.qa.hglm_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGLM_corr")

boyce.qa.hglm_hybrid <- boyce.qa.hglm_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGLM_hybrid")










#####################
### Validate HGPR ###
#####################


# Define vector of covar names
covars <- c("log.bathym","log.npp","log.sst")

# Define 1D meshes to be used for prediction across sites
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)

# nbasis <- 5  #number of basis functions for approximating GP
# degree <- 2  #degree for defining 1D mesh of GP
# alpha <- 2  #for calculating Matern covariance matrix


# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),
                                 degree = 2,
                                 boundary = 'free')
}




### Brazil ###

my.ind.br <- names(cov_list_br$npp)
br.rast.tmp <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
br.rast.hgpr <- list(corr = br.rast.tmp, hybrid = br.rast.tmp)  #store results for both corr and hybrid

names(br.rast.hgpr$corr) <- my.ind.br
names(br.rast.hgpr$hybrid) <- my.ind.br


# Define GP coeff values from HGPR (correlative and hybrid models)
coeff_corr <- hgpr_corr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))
coeff_hybrid <- hgpr_hybrid.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Define coeff values of fixed terms from HGPR (hybrid model)
coeff_hybrid_fixed <- hgpr_hybrid.fit$summary.fixed$mean


# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  gp_vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log()) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(log.bathym, log.npp, log.sst)

  fixed_vars <- data.frame(Intercept = 1,
                      log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                        log(),
                      log.sst2 = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                        log() %>%
                        . ^ 2) %>%
    mutate(row_id = 1:nrow(.)) %>%
    filter(row_id %in% gp_vars$row_id)



  # Generate matrices for covariate raster data (for prediction)
  A.mat <- vector("list", length(covars))
  for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = gp_vars[[covars[[j]]]])
  }


  # Make predictions on intensity of use from model for GP terms
  br.hgpr_hybrid_gp <- A.mat %>%
    map2(.x = ., .y = coeff_hybrid,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    rowSums()  #sum up all predictions across covars

  br.hgpr_corr <- A.mat %>%
    map2(.x = ., .y = coeff_corr,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    rowSums()  #sum up all predictions across covars

  # Make predictions using linear terms
  br.hgpr_hybrid_fixed <- as.matrix(fixed_vars[,1:3]) %*% coeff_hybrid_fixed

  ## Store results in raster stack
  # initially store all NAs for locs w/o predictions
  terra::values(br.rast.hgpr$hybrid[[i]]) <- terra::values(br.rast.hgpr$corr[[i]]) <- NA

  # store predictions
  terra::values(br.rast.hgpr$hybrid[[i]])[gp_vars$row_id] <- br.hgpr_hybrid_gp + br.hgpr_hybrid_fixed[,1]
  terra::values(br.rast.hgpr$corr[[i]])[gp_vars$row_id] <- br.hgpr_corr
}
skrrrahh('khaled2')
toc()  #took 30 sec


# Normalize predictions on 0-1 scale
br.rast.hgpr2 <- br.rast.hgpr |>
  map(normalize)


## Assess model performance via Continuous Boyce Index ##

## Correlative HGPR
boyce.br.full.hgpr_corr <- vector("list", nlyr(br.rast.hgpr2$corr))
boyce.br.sub.hgpr_corr <- vector("list", nlyr(br.rast.hgpr2$corr))
tic()
for (i in 1:nlyr(br.rast.hgpr2$corr)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_main <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hgpr_corr[[i]] <- boyce(fit = br.rast.hgpr2$corr[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.hgpr_corr[[i]] <- boyce(fit = br.rast.hgpr2$corr[[i]],
                                  obs = obs_main,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 2 sec



# perc.use.br.full.hgpr <- boyce.br.full.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.full.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #5 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.full.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()
#
#
#
# perc.use.br.sub.hgpr <- boyce.br.sub.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.sub.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.2 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.sub.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   ylim(0,1) +
#   theme_bw()


boyce.br.full.hgpr_corr <- boyce.br.full.hgpr_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGPR_corr")

boyce.br.sub.hgpr_corr <- boyce.br.sub.hgpr_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_main",
             Method = "HGPR_corr")




## Hybrid HGPR
boyce.br.full.hgpr_hybrid <- vector("list", nlyr(br.rast.hgpr2$hybrid))
boyce.br.sub.hgpr_hybrid <- vector("list", nlyr(br.rast.hgpr2$hybrid))
tic()
for (i in 1:nlyr(br.rast.hgpr2$hybrid)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_main <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hgpr_hybrid[[i]] <- boyce(fit = br.rast.hgpr2$hybrid[[i]],
                                        obs = obs_full,
                                        nbins = 10,
                                        bin.method = "seq",
                                        PEplot = FALSE,
                                        rm.duplicate = FALSE,
                                        method = "spearman")

  boyce.br.sub.hgpr_hybrid[[i]] <- boyce(fit = br.rast.hgpr2$hybrid[[i]],
                                       obs = obs_main,
                                       nbins = 10,
                                       bin.method = "seq",
                                       PEplot = FALSE,
                                       rm.duplicate = FALSE,
                                       method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 2 sec



# perc.use.br.full.hgpr <- boyce.br.full.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.full.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #5 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.full.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()
#
#
#
# perc.use.br.sub.hgpr <- boyce.br.sub.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.br.sub.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.2 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.br.sub.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   ylim(0,1) +
#   theme_bw()


boyce.br.full.hgpr_hybrid <- boyce.br.full.hgpr_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGPR_hybrid")

boyce.br.sub.hgpr_hybrid <- boyce.br.sub.hgpr_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_main",
             Method = "HGPR_hybrid")







### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.tmp <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
qa.rast.hgpr <- list(corr = qa.rast.tmp, hybrid = qa.rast.tmp)  #store results for both corr and hybrid

names(qa.rast.hgpr$corr) <- my.ind.qa
names(qa.rast.hgpr$hybrid) <- my.ind.qa


# Define GP coeff values from HGPR (correlative and hybrid models)
coeff_corr <- hgpr_corr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))
coeff_hybrid <- hgpr_hybrid.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Define coeff values of fixed terms from HGPR (hybrid model)
coeff_hybrid_fixed <- hgpr_hybrid.fit$summary.fixed$mean



# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  gp_vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log()) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(log.bathym, log.npp, log.sst)

  fixed_vars <- data.frame(Intercept = 1,
                      log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                        log(),
                      log.sst2 = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                        log() %>%
                        . ^ 2) %>%
    mutate(row_id = 1:nrow(.)) %>%
    filter(row_id %in% gp_vars$row_id)



  # Generate matrices for covariate raster data (for prediction)
  A.mat <- vector("list", length(covars))
  for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = gp_vars[[covars[[j]]]])
  }


  # Make predictions on intensity of use from model for GP terms
  qa.hgpr_hybrid_gp <- A.mat %>%
    map2(.x = ., .y = coeff_hybrid,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    rowSums()  #sum up all predictions across covars

  qa.hgpr_corr <- A.mat %>%
    map2(.x = ., .y = coeff_corr,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
    rowSums()  #sum up all predictions across covars

  # Make predictions using linear terms
  qa.hgpr_hybrid_fixed <- as.matrix(fixed_vars[,1:3]) %*% coeff_hybrid_fixed

  ## Store results in raster stack
  # initially store all NAs for locs w/o predictions
  terra::values(qa.rast.hgpr$hybrid[[i]]) <- terra::values(qa.rast.hgpr$corr[[i]]) <- NA

  # store predictions
  terra::values(qa.rast.hgpr$hybrid[[i]])[gp_vars$row_id] <- qa.hgpr_hybrid_gp + qa.hgpr_hybrid_fixed[,1]
  terra::values(qa.rast.hgpr$corr[[i]])[gp_vars$row_id] <- qa.hgpr_corr
}
skrrrahh('khaled2')
toc()  #took 1 sec


# Normalize predictions on 0-1 scale
qa.rast.hgpr2 <- qa.rast.hgpr |>
  map(normalize)


## Assess model performance via Continuous Boyce Index ##

## Correlative HGPR
boyce.qa.hgpr_corr <- vector("list", nlyr(qa.rast.hgpr2$corr))
tic()
for (i in 1:nlyr(qa.rast.hgpr2$corr)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hgpr_corr[[i]] <- boyce(fit = qa.rast.hgpr2$corr[[i]],
                              obs = obs,
                              nbins = 10,
                              bin.method = "seq",
                              PEplot = FALSE,
                              rm.duplicate = FALSE,
                              method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


# perc.use.qa.hgpr <- boyce.qa.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.qa.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.6 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.qa.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()



boyce.qa.hgpr_corr <- boyce.qa.hgpr_corr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGPR_corr")






## Hybrid HGPR
boyce.qa.hgpr_hybrid <- vector("list", nlyr(qa.rast.hgpr2$hybrid))
tic()
for (i in 1:nlyr(qa.rast.hgpr2$hybrid)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hgpr_hybrid[[i]] <- boyce(fit = qa.rast.hgpr2$hybrid[[i]],
                                   obs = obs,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


# perc.use.qa.hgpr <- boyce.qa.hgpr %>%
#   map(., pluck, "perc.use") %>%
#   set_names(1:length(.)) %>%
#   bind_rows() %>%
#   janitor::remove_empty(which = "cols") %>%
#   apply(., 2, function(x) cumsum(rev(x)))
#
# # check fewest bins that contain >=90% of all obs
# apply(perc.use.qa.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
#   mean()  #2.6 bins
#
# # Viz plot of cumulative percentage of obs per bin (highest to lowest)
# perc.use.qa.hgpr %>%
#   data.frame() %>%
#   mutate(bin = factor(10:1, levels = 10:1)) %>%
#   pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
#   ggplot(aes(bin, cum.perc)) +
#   geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
#   geom_line(aes(group = month.year, color = month.year)) +
#   theme_bw()



boyce.qa.hgpr_hybrid <- boyce.qa.hgpr_hybrid %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGPR_hybrid")







####################################
### Summarize Validation Results ###
####################################

boyce.fit <- rbind(boyce.br.full.hglm_corr, boyce.br.sub.hglm_corr, boyce.qa.hglm_corr,
                   boyce.br.full.hglm_hybrid, boyce.br.sub.hglm_hybrid, boyce.qa.hglm_hybrid,
                   boyce.br.full.hgpr_corr, boyce.br.sub.hgpr_corr, boyce.qa.hgpr_corr,
                   boyce.br.full.hgpr_hybrid, boyce.br.sub.hgpr_hybrid, boyce.qa.hgpr_hybrid)

boyce.mean <- boyce.fit %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(Method, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.mean, aes(x = Region, y = mean, group = Method),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_fill_met_d('Egypt', labels = c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)")) +
  # scale_fill_manual(values = c("#dd5129",
  #                              scales::muted("#dd5129", l=80),
  #                              "#0f7ba2",
  #                              scales::muted("#0f7ba2", l=80)),
  #                   labels = c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)")) +
  scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))




# Create example prediction maps per method

bbox <- ext(br.rast.hglm$corr)

# time-matched observations
tmp.br <- dat.br %>%
  filter(month.year == "2022-02-01")

# Break rasters into bins used for Boyce Index
br.rast.hglm_corr.d <- classify(br.rast.hglm2$corr, seq(0, 1, by = 0.1))
br.rast.hglm_corr.d <- br.rast.hglm_corr.d + 1
br.rast.hglm_corr.df <- as.data.frame(br.rast.hglm_corr.d, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1))) %>%
  mutate(Method = "HGLM (corr)")

br.rast.hglm_hybrid.d <- classify(br.rast.hglm2$hybrid, seq(0, 1, by = 0.1))
br.rast.hglm_hybrid.d <- br.rast.hglm_hybrid.d + 1
br.rast.hglm_hybrid.df <- as.data.frame(br.rast.hglm_hybrid.d, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1))) %>%
  mutate(Method = "HGLM (hybrid)")

br.rast.hgpr_corr.d <- classify(br.rast.hgpr2$corr, seq(0, 1, by = 0.1))
br.rast.hgpr_corr.d <- br.rast.hgpr_corr.d + 1
br.rast.hgpr_corr.df <- as.data.frame(br.rast.hgpr_corr.d, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1))) %>%
  mutate(Method = "HGPR (corr)")

br.rast.hgpr_hybrid.d <- classify(br.rast.hgpr2$hybrid, seq(0, 1, by = 0.1))
br.rast.hgpr_hybrid.d <- br.rast.hgpr_hybrid.d + 1
br.rast.hgpr_hybrid.df <- as.data.frame(br.rast.hgpr_hybrid.d, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1))) %>%
  mutate(Method = "HGPR (hybrid)")



# p.hglm.br <- ggplot() +
#   geom_raster(data = br.rast.hglm.df, aes(x, y, fill = `2022-02-01`)) +
#   scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
#   geom_sf(data = br.sf) +
#   # geom_point(data = tmp.br, aes(x, y), color = "blue", alpha = 0.7, size = 1) +
#   labs(x="",y="", title = "HGLM") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.hgam.br <- ggplot() +
#   geom_raster(data = br.rast.hgam.df, aes(x, y, fill = `2022-02-01`)) +
#   scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
#   geom_sf(data = br.sf) +
#   # geom_point(data = tmp.pts, aes(x, y), color = "blue", alpha = 0.7, size = 1) +
#   labs(x="",y="", title = "HGAM") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.brt.br <- ggplot() +
#   geom_raster(data = br.rast.brt.df, aes(x, y, fill = `2022-02-01`)) +
#   scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
#   geom_sf(data = br.sf) +
#   # geom_point(data = tmp.pts, aes(x, y), color = "blue", alpha = 0.7, size = 1) +
#   labs(x="",y="", title = "BRT") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.hgpr.br <- ggplot() +
#   geom_raster(data = br.rast.hgpr.df, aes(x, y, fill = `2022-02-01`)) +
#   scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
#   geom_sf(data = br.sf) +
#   # geom_point(data = tmp.pts, aes(x, y), color = "blue", alpha = 0.7, size = 1) +
#   labs(x="",y="", title = "HGPR") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
#
# # Make composite plot
# p.hglm.br + p.hgam.br + p.brt.br + p.hgpr.br +
#   plot_layout(ncol = 2, guides = "collect")
#
# ggsave("Tables_Figs/Figure 4.png", width = 4, height = 5, units = "in", dpi = 400)


spat.preds <- rbind(br.rast.hglm_corr.df,
                    br.rast.hglm_hybrid.df,
                    br.rast.hgpr_corr.df,
                    br.rast.hgpr_hybrid.df) %>%
  mutate(across(Method, \(x) factor(x, levels = c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)"))))

ggplot() +
  geom_raster(data = spat.preds, aes(x, y, fill = `2022-02-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = br.sf) +
  geom_text(data = data.frame(x = c(-4802714,0,0,0),
                              y = c(-1000000,0,0,0),
                              Method = factor(c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)"),
                                              levels = c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)")),
                              lab = c("Brazil","","","")), aes(x, y, label = lab),
            size = 5, fontface = "italic") +
  geom_point(data = tmp.br, aes(x, y), fill = "chartreuse", alpha = 0.6, size = 1,
             shape = 21, stroke = 0.25) +
  labs(x="",y="", title = "February 2022") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = 0)) +
  facet_wrap(~ Method)

# ggsave("Tables_Figs/Figure 2.png", width = 4, height = 5, units = "in", dpi = 400)





###############################################################
### Create composite plot of marginal effects across models ###
###############################################################

# Load analysis dataset
rsf.pts_10s <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv") %>%
  drop_na(bathym, npp, sst) %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))


### HGLM ###

# Store model coeffs
hglm.coeffs <- list(corr = hglm_corr.fit$summary.fixed[-1,],
                    hybrid = hglm_hybrid.fit$summary.fixed[-1,]) %>%
  map(~{.x %>%
      mutate(param = factor(rownames(.x), levels = rownames(.x)))}) %>%
  bind_rows(.id = "dataset") %>%
  mutate(coeff = rep(rownames(hglm_corr.fit$summary.fixed[-1,]), 2)) %>%
  select(dataset, mean, coeff) %>%
  pivot_wider(names_from = dataset, values_from = mean, id_cols = coeff) %>%
  select(-coeff) %>%
  as.matrix()

# Create seq of predictor values
newdat.hglm <- list(depth = data.frame(bathym = seq(min(rsf.pts_10s$log.bathym),
                                                     max(rsf.pts_10s$log.bathym),
                                                     length.out = 500),
                                        bathym.2 = seq(min(rsf.pts_10s$log.bathym),
                                                       max(rsf.pts_10s$log.bathym),
                                                       length.out = 500) ^ 2,
                                        npp = 0,
                                        npp.2 = 0,
                                        sst = 0,
                                        sst.2 = 0) %>%
                      as.matrix(),
                    npp = data.frame(bathym = 0,
                                     bathym.2 = 0,
                                     npp = seq(min(rsf.pts_10s$log.npp),
                                               max(rsf.pts_10s$log.npp),
                                               length.out = 500),
                                     npp.2 = seq(min(rsf.pts_10s$log.npp),
                                                 max(rsf.pts_10s$log.npp),
                                                 length.out = 500) ^ 2,
                                     sst = 0,
                                     sst.2 = 0) %>%
                      as.matrix(),
                    sst = data.frame(bathym = 0,
                                     bathym.2 = 0,
                                     npp = 0,
                                     npp.2 = 0,
                                     sst = seq(min(rsf.pts_10s$log.sst),
                                               max(rsf.pts_10s$log.sst),
                                               length.out = 500),
                                     sst.2 = seq(min(rsf.pts_10s$log.sst),
                                                 max(rsf.pts_10s$log.sst),
                                                 length.out = 500) ^ 2) %>%
                      as.matrix())


# Store predictions
marg.eff.hglm <- vector("list", length = 3) %>%
  set_names(names(newdat.hglm))

for (i in seq_along(newdat.hglm)) {

  # Find first column per covariate to add for x-axis
  ind.col <- which(newdat.hglm[[i]][1,] != 0) %>%
    min()

  marg.eff.hglm[[i]] <- newdat.hglm[[i]] %*% hglm.coeffs %>%
    data.frame() %>%
    mutate(x = exp(newdat.hglm[[i]][,ind.col])) %>%
    pivot_longer(cols = -x, names_to = "Method", values_to = "mean") |>
    mutate(Method = paste0("HGLM (", Method, ")"))
}

marg.eff.hglm2 <- bind_rows(marg.eff.hglm, .id = "covar") %>%
  filter(!(covar == "depth" & x > 300))  #to constrain plot to only show depth up to 300 m
marg.eff.hglm2$mean <- exp(marg.eff.hglm2$mean)  #change to natural response scale






### HGPR ###

# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)


# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),
                                 degree = 2,
                                 boundary = 'free')
}

# Generate matrices for sequences of covariate ranges
A.me.pop <- vector("list", length(covars))
newdat.list <- map(mesh.seq, function(x) {seq(from = x[1], to = x[2], length.out = 500)})
for (i in 1:length(covars)) { #make marginal predictions to viz response by covar
  A.me.pop[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.list[[covars[[i]]]])
}


# Store resulting GP coeffs per covar into a list
hgpr.coeffs_hybrid <- hgpr_hybrid.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))
hgpr.coeffs_corr <- hgpr_corr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Make predictions via linear algebra
marg.eff.hgpr_hybrid <- A.me.pop %>%
  map2(.x = ., .y = hgpr.coeffs_hybrid,
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
       }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list),
         Method = "HGPR (hybrid)") %>%
  mutate(across(mean:x, exp)) %>%
  mutate(covar = case_when(covar == "log.bathym" ~ "depth",
                           covar == "log.npp" ~ "npp",
                           covar == "log.sst" ~ "sst")) %>%
  filter(!(covar == "depth" & x > 300)) %>%  #to constrain plot to only show depth up to 300 m
  filter(!(covar == "npp" & x > max(rsf.pts_10s$npp))) %>%  #to match range of other models
  filter(!(covar == "sst" & x > max(rsf.pts_10s$sst) |
             covar == "sst" & x < min(rsf.pts_10s$sst)))  #to match range of other models

marg.eff.hgpr_corr <- A.me.pop %>%
  map2(.x = ., .y = hgpr.coeffs_corr,
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
       }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list),
         Method = "HGPR (corr)") %>%
  mutate(across(mean:x, exp)) %>%
  mutate(covar = case_when(covar == "log.bathym" ~ "depth",
                           covar == "log.npp" ~ "npp",
                           covar == "log.sst" ~ "sst")) %>%
  filter(!(covar == "depth" & x > 300)) %>%  #to constrain plot to only show depth up to 300 m
  filter(!(covar == "npp" & x > max(rsf.pts_10s$npp))) %>%  #to match range of other models
  filter(!(covar == "sst" & x > max(rsf.pts_10s$sst) |
             covar == "sst" & x < min(rsf.pts_10s$sst)))  #to match range of other models


# Bind both model predictions together
marg.eff.hgpr <- rbind(marg.eff.hgpr_corr,
                       marg.eff.hgpr_hybrid)


# Make predictions of linear SST terms
sst.newdata <- data.frame(sst = newdat.list$log.sst,
                          sst.2 = newdat.list$log.sst ^ 2) %>%
  as.matrix()
fixed.sst.coeff <- hgpr_hybrid.fit$summary.fixed$mean[-1]

fixed.sst.pred <- sst.newdata %*% fixed.sst.coeff %>%
  data.frame(mean = .) %>%
  mutate(sst = exp(sst.newdata[,1]),
         mean = exp(mean)) %>%
  filter(!(sst > max(rsf.pts_10s$sst) |
             sst < min(rsf.pts_10s$sst)))  #to match range of other models

# Add linear SST predictions to GP SST predictions
marg.eff.hgpr2 <- marg.eff.hgpr
marg.eff.hgpr2[marg.eff.hgpr2$covar == "sst" &
                 marg.eff.hgpr2$Method == "HGPR (hybrid)",]$mean <- marg.eff.hgpr2[marg.eff.hgpr2$covar == "sst" &
                                                                                     marg.eff.hgpr2$Method == "HGPR (hybrid)",]$mean +
  fixed.sst.pred$mean




### Combine all predictions into single DF and plot

marg.eff <- rbind(marg.eff.hglm2,
                 marg.eff.hgpr2) %>%
  mutate(method = factor(Method, levels = c("HGLM (corr)","HGLM (hybrid)","HGPR (corr)","HGPR (hybrid)")),
         x = case_when(covar == "npp" ~ x / 1000,
                       TRUE ~ x),
         covar = case_when(covar == "depth" ~ "Depth (m)",
                           covar == "npp" ~ "NPP (g C m^-2 d^-1)",
                           covar == "sst" ~ "SST (°C)"))



# Plot pop-level marginal effects
ggplot() +
  geom_line(data = marg.eff, aes(x = x, y = mean, color = Method), linewidth = 1) +
  scale_color_met_d(palette_name = 'Egypt', guide = "none") +
  theme_bw(base_size = 14) +
  labs(x = "", y = "Relative Intensity of Use") +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank()
        ) +
  # facet_wrap(Method ~ covar, ncol = 3, scales = "free", strip.position = "bottom")
  ggh4x::facet_grid2(Method ~ covar,
                     independent = "y",
                     scales = "free",
                     switch = "x")

# ggsave("Tables_Figs/Figure S5.png", width = 11, height = 9, units = "in", dpi = 400)



### Export Boyce Index results ###

# write_csv(boyce.fit, "Data_products/boyce_alg_results.csv")
