
### Compare model transferability by method ###

library(tidyverse)
library(INLA)
library(mgcv)
library(gbm)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)

source('Scripts/helper functions.R')


##########################
### Load fitted models ###
##########################

hglm.fit <- readRDS("Data_products/HGLM_model_fit.rds")
hgam.fit <- readRDS("Data_products/HGAM_model_fit.rds")
brt.fit <- readRDS("Data_products/BRT_model_fit.rds")
hgpr.fit <- readRDS("Data_products/HGPR_model_fit.rds")




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







#####################
### Validate HGLM ###
#####################

### Brazil ###

my.ind.br <- names(cov_list_br$npp)
br.rast.hglm <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
names(br.rast.hglm) <- my.ind.br

# Define coeff values from HGLM
coeff1 <- hglm.fit$summary.fixed$mean

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


  # Make predictions on intensity of use from model
  br.hglm <- as.matrix(vars) %*% coeff1  #make predictions

  # Store results in raster stack
  terra::values(br.rast.hglm[[i]]) <- br.hglm  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 1 min


# Normalize predictions on 0-1 scale
br.rast.hglm2 <- normalize(br.rast.hglm)


# Assess model performance via Continuous Boyce Index
boyce.br.full.hglm <- vector("list", nlyr(br.rast.hglm2))
boyce.br.sub.hglm <- vector("list", nlyr(br.rast.hglm2))
tic()
for (i in 1:nlyr(br.rast.hglm2)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hglm[[i]] <- boyce(fit = br.rast.hglm2[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.hglm[[i]] <- boyce(fit = br.rast.hglm2[[i]],
                                  obs = obs_sub,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 5 sec



perc.use.br.full.hglm <- boyce.br.full.hglm %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.full.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #2.6 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.hglm %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



perc.use.br.sub.hglm <- boyce.br.sub.hglm %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.sub.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #1.3 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.hglm %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  ylim(0,1) +
  theme_bw()


boyce.br.full.hglm <- boyce.br.full.hglm %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGLM")

boyce.br.sub.hglm <- boyce.br.sub.hglm %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_sub",
             Method = "HGLM")





### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.hglm <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.hglm) <- my.ind.qa

# Define coeff values from HGLM
coeff1 <- hglm.fit$summary.fixed$mean

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


  # Make predictions on intensity of use from model
  qa.hglm <- as.matrix(vars) %*% coeff1  #make predictions

  # Store results in raster stack
  terra::values(qa.rast.hglm[[i]]) <- qa.hglm  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 1 sec


# Normalize predictions on 0-1 scale
qa.rast.hglm2 <- normalize(qa.rast.hglm)


# Assess model performance via Continuous Boyce Index
boyce.qa.hglm <- vector("list", nlyr(qa.rast.hglm2))
tic()
for (i in 1:nlyr(qa.rast.hglm2)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hglm[[i]] <- boyce(fit = qa.rast.hglm2[[i]],
                              obs = obs,
                              nbins = 10,
                              bin.method = "seq",
                              PEplot = FALSE,
                              rm.duplicate = FALSE,
                              method = "spearman")
  }
skrrrahh("khaled3")
toc()  #took 1 sec


perc.use.qa.hglm <- boyce.qa.hglm %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.qa.hglm, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #2.6 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.hglm %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()

boyce.qa.hglm <- boyce.qa.hglm %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGLM")





#####################
### Validate HGAM ###
#####################

### Brazil ###

my.ind.br <- names(cov_list_br$npp)
br.rast.hgam <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
names(br.rast.hgam) <- my.ind.br

tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log(),
                     id = hgam.fit$var.summary$id)


  # Make predictions on intensity of use from model
  br.hgam <- predict.bam(hgam.fit, newdata = vars, type = "terms", terms = c("s(log.bathym)","s(log.npp)","s(log.sst)"),
                         discrete = FALSE, na.action = na.pass)

  terra::values(br.rast.hgam[[i]]) <- rowSums(br.hgam + hgam.fit$coefficients[1])  #add intercept

}
skrrrahh('khaled2')
toc()  #took 1 min


# Normalize predictions on 0-1 scale
br.rast.hgam2 <- normalize(br.rast.hgam)


# Assess model performance via Continuous Boyce Index
boyce.br.full.hgam <- vector("list", nlyr(br.rast.hgam2))
boyce.br.sub.hgam <- vector("list", nlyr(br.rast.hgam2))
tic()
for (i in 1:nlyr(br.rast.hgam2)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hgam[[i]] <- boyce(fit = br.rast.hgam2[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.hgam[[i]] <- boyce(fit = br.rast.hgam2[[i]],
                                  obs = obs_sub,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 5 sec



perc.use.br.full.hgam <- boyce.br.full.hgam %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.full.hgam, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #4.4 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.hgam %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



perc.use.br.sub.hgam <- boyce.br.sub.hgam %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.sub.hgam, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #3.3 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.hgam %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()


boyce.br.full.hgam <- boyce.br.full.hgam %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGAM")

boyce.br.sub.hgam <- boyce.br.sub.hgam %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_sub",
             Method = "HGAM")





### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.hgam <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.hgam) <- my.ind.qa


# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log(),
                     id = hgam.fit$var.summary$id)


  # Make predictions on intensity of use from model
  qa.hgam <- predict.bam(hgam.fit, newdata = vars, type = "terms", terms = c("s(log.bathym)","s(log.npp)","s(log.sst)"),
                         discrete = FALSE, na.action = na.pass)

  terra::values(qa.rast.hgam[[i]]) <- rowSums(qa.hgam + hgam.fit$coefficients[1])  #add intercept

}
skrrrahh('khaled2')
toc()  #took 1 sec


# Normalize predictions on 0-1 scale
qa.rast.hgam2 <- normalize(qa.rast.hgam)


# Assess model performance via Continuous Boyce Index
boyce.qa.hgam <- vector("list", nlyr(qa.rast.hgam2))
tic()
for (i in 1:nlyr(qa.rast.hgam2)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hgam[[i]] <- boyce(fit = qa.rast.hgam2[[i]],
                              obs = obs,
                              nbins = 10,
                              bin.method = "seq",
                              PEplot = FALSE,
                              rm.duplicate = FALSE,
                              method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


perc.use.qa.hgam <- boyce.qa.hgam %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.qa.hgam, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #4.5 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.hgam %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()

boyce.qa.hgam <- boyce.qa.hgam %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGAM")






####################
### Validate BRT ###
####################

### Brazil ###

my.ind.br <- names(cov_list_br$npp)
br.rast.brt <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
terra::values(br.rast.brt) <- NA  # initially store all NAs for locs w/o predictions
names(br.rast.brt) <- my.ind.br

tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  print(paste0(i,"/",nlyr(cov_list_br$npp)))

  # Subset covars by month.year
  vars <- data.frame(bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs(),
                     npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])),
                     sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]]))) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(bathym, npp, sst)


  # Make predictions on intensity of use from model
  br.brt <- predict.gbm(brt.fit, newdata = vars[,-4], n.trees = brt.fit$n.trees)
  terra::values(br.rast.brt[[i]])[vars$row_id] <- br.brt  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 5.5 min


# Normalize predictions on 0-1 scale
br.rast.brt2 <- normalize(br.rast.brt)


# Assess model performance via Continuous Boyce Index
boyce.br.full.brt <- vector("list", nlyr(br.rast.brt2))
boyce.br.sub.brt <- vector("list", nlyr(br.rast.brt2))
tic()
for (i in 1:nlyr(br.rast.brt2)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.brt[[i]] <- boyce(fit = br.rast.brt2[[i]],
                                  obs = obs_full,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")

  boyce.br.sub.brt[[i]] <- boyce(fit = br.rast.brt2[[i]],
                                 obs = obs_sub,
                                 nbins = 10,
                                 bin.method = "seq",
                                 PEplot = FALSE,
                                 rm.duplicate = FALSE,
                                 method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 5 sec


perc.use.br.full.brt <- boyce.br.full.brt %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.full.brt, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #9.1 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.brt %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



perc.use.br.sub.brt <- boyce.br.sub.brt %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.sub.brt, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #8.1 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.brt %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  ylim(0,1) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



boyce.br.full.brt <- boyce.br.full.brt %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "BRT")

boyce.br.sub.brt <- boyce.br.sub.brt %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_sub",
             Method = "BRT")




### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.brt <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.brt) <- my.ind.qa

tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  vars <- data.frame(bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs(),
                     npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])),
                     sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]]))) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(bathym, npp, sst)


  # Make predictions on intensity of use from model
  qa.brt <- predict.gbm(brt.fit, newdata = vars[,-4], n.trees = brt.fit$n.trees)

  terra::values(qa.rast.brt[[i]]) <- NA  # initially store all NAs for locs w/o predictions
  terra::values(qa.rast.brt[[i]])[vars$row_id] <- qa.brt  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 2 sec


# Normalize predictions on 0-1 scale
qa.rast.brt2 <- normalize(qa.rast.brt)



boyce.qa.brt <- vector("list", nlyr(qa.rast.brt2))
tic()
for (i in 1:nlyr(qa.rast.brt2)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.brt[[i]] <- boyce(fit = qa.rast.brt2[[i]],
                             obs = obs,
                             nbins = 10,
                             bin.method = "seq",
                             PEplot = FALSE,
                             rm.duplicate = FALSE,
                             method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


perc.use.qa.brt <- boyce.qa.brt %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.qa.brt, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #6.5 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.brt %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()

boyce.qa.brt <- boyce.qa.brt %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "BRT")









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
alpha <- 2  #for calculating Matern covariance matrix


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
br.rast.hgpr <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
names(br.rast.hgpr) <- my.ind.br

# Define coeff values from HGPR
coeff1 <- hgpr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Define coeff values of fixed terms from HGPR
coeff2 <- hgpr.fit$summary.fixed$mean


# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log()) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(log.bathym, log.npp, log.sst)

  vars2 <- data.frame(Intercept = 1,
                      log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                        log(),
                      log.sst2 = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                        log() %>%
                        . ^ 2) %>%
    mutate(row_id = 1:nrow(.)) %>%
    filter(row_id %in% vars$row_id)



  # Generate matrices for covariate raster data (for prediction)
  A.mat <- vector("list", length(covars))
  for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = vars[[covars[[j]]]])
  }


  # Make predictions on intensity of use from model for GP terms
  br.hgpr <- A.mat %>%
    map2(.x = ., .y = coeff1,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols() %>%
    rowSums()  #sum up all predictions across covars

  # Make predictions using linear terms
  br.hgpr2 <- as.matrix(vars2[,1:3]) %*% coeff2

  # Store results in raster stack
  terra::values(br.rast.hgpr[[i]]) <- NA  # initially store all NAs for locs w/o predictions
  terra::values(br.rast.hgpr[[i]])[vars$row_id] <- br.hgpr + br.hgpr2[,1]
}
skrrrahh('khaled2')
toc()  #took 40 sec


# Normalize predictions on 0-1 scale
br.rast.hgpr2 <- normalize(br.rast.hgpr)


# Assess model performance via Continuous Boyce Index
boyce.br.full.hgpr <- vector("list", nlyr(br.rast.hgpr2))
boyce.br.sub.hgpr <- vector("list", nlyr(br.rast.hgpr2))
tic()
for (i in 1:nlyr(br.rast.hgpr2)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%
    dplyr::select(x, y)

  boyce.br.full.hgpr[[i]] <- boyce(fit = br.rast.hgpr2[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.hgpr[[i]] <- boyce(fit = br.rast.hgpr2[[i]],
                                  obs = obs_sub,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 5 sec



perc.use.br.full.hgpr <- boyce.br.full.hgpr %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.full.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #5 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.hgpr %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



perc.use.br.sub.hgpr <- boyce.br.sub.hgpr %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.sub.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #2.2 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.hgpr %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  ylim(0,1) +
  theme_bw()


boyce.br.full.hgpr <- boyce.br.full.hgpr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_all",
             Method = "HGPR")

boyce.br.sub.hgpr <- boyce.br.sub.hgpr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil_sub",
             Method = "HGPR")





### Qatar ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.hgpr <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.hgpr) <- my.ind.qa

# Define coeff values from HGPR
coeff1 <- hgpr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Define coeff values of fixed terms from HGPR
coeff2 <- hgpr.fit$summary.fixed$mean



# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log()) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(log.bathym, log.npp, log.sst)

  vars2 <- data.frame(Intercept = 1,
                      log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                        log(),
                      log.sst2 = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                        log() %>%
                        . ^ 2) %>%
    mutate(row_id = 1:nrow(.)) %>%
    filter(row_id %in% vars$row_id)



  # Generate matrices for covariate raster data (for prediction)
  A.mat <- vector("list", length(covars))
  for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = vars[[covars[[j]]]])
  }


  # Make predictions on intensity of use from model for GP terms
  qa.hgpr <- A.mat %>%
    map2(.x = ., .y = coeff1,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols() %>%
    rowSums()  #sum up all predictions across covars

  # Make predictions using linear terms
  qa.hgpr2 <- as.matrix(vars2[,1:3]) %*% coeff2

  # Store results in raster stack
  terra::values(qa.rast.hgpr[[i]]) <- NA  # initially store all NAs for locs w/o predictions
  terra::values(qa.rast.hgpr[[i]])[vars$row_id] <- qa.hgpr + qa.hgpr2[,1]

}
skrrrahh('khaled2')
toc()  #took 2 sec


# Normalize predictions on 0-1 scale
qa.rast.hgpr2 <- normalize(qa.rast.hgpr)


# Assess model performance via Continuous Boyce Index
boyce.qa.hgpr <- vector("list", nlyr(qa.rast.hgpr2))
tic()
for (i in 1:nlyr(qa.rast.hgpr2)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.hgpr[[i]] <- boyce(fit = qa.rast.hgpr2[[i]],
                              obs = obs,
                              nbins = 10,
                              bin.method = "seq",
                              PEplot = FALSE,
                              rm.duplicate = FALSE,
                              method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


perc.use.qa.hgpr <- boyce.qa.hgpr %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.qa.hgpr, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #2.6 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.hgpr %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()

boyce.qa.hgpr <- boyce.qa.hgpr %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar",
             Method = "HGPR")







####################################
### Summarize Validation Results ###
####################################

boyce.fit <- rbind(boyce.br.full.hglm, boyce.br.sub.hglm, boyce.qa.hglm,
                   boyce.br.full.hgam, boyce.br.sub.hgam, boyce.qa.hgam,
                   boyce.br.full.brt, boyce.br.sub.brt, boyce.qa.brt,
                   boyce.br.full.hgpr, boyce.br.sub.hgpr, boyce.qa.hgpr)

boyce.mean <- boyce.fit %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 5, position = position_dodge(width = 0.75)) +
  geom_violin(aes(color = Method), fill = "transparent", position = position_dodge(width = 0.75)) +
  geom_point(data = boyce.mean, aes(x = Region, y = mean, group = Method),
             size = 6, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24))





cum.perc <- list(Brazil_all.HGLM = perc.use.br.full.hglm,
                 Brazil_sub.HGLM = perc.use.br.sub.hglm,
                 Qatar.HGLM = perc.use.qa.hglm,
                 Brazil_all.HGAM = perc.use.br.full.hgam,
                 Brazil_sub.HGAM = perc.use.br.sub.hgam,
                 Qatar.HGAM = perc.use.qa.hgam,
                 Brazil_all.BRT = perc.use.br.full.brt,
                 Brazil_sub.BRT = perc.use.br.sub.brt,
                 Qatar.BRT = perc.use.qa.brt,
                 Brazil_all.HGPR = perc.use.br.full.hgpr,
                 Brazil_sub.HGPR = perc.use.br.sub.hgpr,
                 Qatar.HGPR = perc.use.qa.hgpr)

cum.perc.mean <- cum.perc %>%
  map(., ~apply(.x, 2, function(x) which(x >= 0.9)[1])) %>%
  map(., ~{data.frame(mean = mean(.x),
                sd = sd(.x)
                )}
      ) %>%
  bind_rows(.id = "id") %>%
  # pivot_longer(cols = everything(), names_to = "id", values_to = "cum.perc") %>%
  separate(col = id, sep = "\\.", into = c("Region","Method"))



ggplot(data = cum.perc.mean, aes(Region, mean)) +
  geom_linerange(aes(ymin = (mean - sd), ymax = (mean + sd), color = Method),
                 position = position_dodge(width = 0.75)) +
  geom_point(aes(group = Method, color = Method),
             size = 6, position = position_dodge(width = 0.75)) +
  # lims(y = c(0,10)) +
  labs(x="", y = "Avg # of bins accounting for 90% of obs") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 22))
