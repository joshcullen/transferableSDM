
### Fit RSF as Boosted Regression Tree ###

library(tidyverse)
library(gbm)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)

source('Scripts/helper functions.R')


#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
# rsf.pts_30 <- read_csv("Processed_data/GoM_Cm_RSFprep_30x.csv")
# rsf.pts_50 <- read_csv("Processed_data/GoM_Cm_RSFprep_50x.csv")

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



####################
### Fit BRT RSF ###
####################

# Remove rows w/ incomplete observations and log-transform skewed covars
rsf.pts_10s <- rsf.pts_10 %>%
  # drop_na(bathym, npp, sst) %>%
  mutate(bathym = abs(bathym),
         log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)



## Define params to optimize for model fitting
hyper_grid <- expand.grid(
  n.trees = c(1000, 5000),
  shrinkage = c(0.001, 0.01),
  interaction.depth = c(1, 3),
  bag.fraction = c(0.5, 0.75),
  min_RMSE = 0                     # a place to dump results
)


# Fit the model
tic()
for(i in 1:nrow(hyper_grid)) {
  print(i)

  # reproducibility
  set.seed(2023)

  gbm.fit <- gbm(
    formula = obs/wts2 ~ bathym + npp + sst,
    distribution = "poisson",
    data = rsf.pts_10s,
    n.trees = hyper_grid$n.trees[i],
    interaction.depth = hyper_grid$interaction.depth[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    shrinkage = hyper_grid$shrinkage[i],
    weights = rsf.pts_10s$wts2
  )

  # add min training error and trees to grid
  # hyper_grid$optimal_trees[i] <- which.min(gbm.fit$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.fit$train.error))
}
toc()  # took 29 min to run

hyper_grid %>%
  dplyr::arrange(min_RMSE)
# best model has 5000 trees, 0.01 shrinkage, depth of 3, and bag fraction of 0.5


# Fit model w/ optimized params
tic()
gbm.fit <- gbm(
  formula = obs/wts2 ~ bathym + npp + sst,
  distribution = "poisson",
  data = rsf.pts_10s,
  n.trees = 5000,
  interaction.depth = 3,
  bag.fraction = 0.5,
  shrinkage = 0.01,
  weights = rsf.pts_10s$wts2

)
toc()  #took 3.8 min to run

print(gbm.fit)
summary(gbm.fit)


pred.pop <- gbm.pdp(gbm.fit, continuous.resolution = 200)

ggplot() +
  geom_line(data = pred.pop, aes(x = x, y = exp(y)), linewidth = 1) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~covar, scales = "free", ncol = 1)







#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[!grepl(pattern = "Kd490", files)]  #remove Kd490 datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

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




##########################################
### Create example predictive surfaces ###
##########################################

# For pop mean in 2020-09
newdat <- data.frame(bathym = terra::values(cov_list$bathym) %>%
                       abs(),
                     npp = terra::values(cov_list$npp$`2020-09-01`),
                     sst = terra::values(cov_list$sst$`2020-09-01`))
names(newdat)[1:3] <- c('bathym','npp','sst')
summary(newdat)



mean.pred <- predict.gbm(gbm.fit, newdata = newdat, n.trees = gbm.fit$n.trees)
summary(mean.pred)

rast.pred <- cov_list$bathym
terra::values(rast.pred) <- exp(mean.pred)
plot(rast.pred)


rast.pred.df <- as.data.frame(rast.pred, xy = TRUE)
names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)

ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = log(pred))) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "Population Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20))








#####################
### Validate HGAM ###
#####################

dat.br <- read_csv("Processed_data/Brazil_Cm_Tracks_behav.csv")
dat.qa <- read_csv("Processed_data/Qatar_Cm_Tracks_behav.csv")

# Create indexing column "month.year"
dat.br <- dat.br %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  filter(behav == 'Resident')

dat.qa <- dat.qa %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  filter(behav == 'Resident')






### Brazil ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Brazil", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_br <- sapply(files, rast)
cov_list_br

names(cov_list_br) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_br[[var]]) <- gsub(names(cov_list_br[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list_br[[var]] <- terra::resample(cov_list_br[[var]], cov_list_br$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list_br <- map(cov_list_br, terra::project, 'EPSG:3395')






### Validate model using CBI ###

my.ind.br <- names(cov_list_br$npp)
br.rast.pred <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
terra::values(br.rast.pred) <- NA  # initially store all NAs for locs w/o predictions
names(br.rast.pred) <- my.ind.br

tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  vars <- data.frame(bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs(),
                     npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])),
                     sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]]))) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(bathym, npp, sst)


  # Make predictions on intensity of use from model
  br.pred <- predict.gbm(gbm.fit, newdata = vars[,-4], n.trees = gbm.fit$n.trees)
  terra::values(br.rast.pred[[i]]) <- br.pred  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 8 min



cbi.br <- vector("list", nlyr(br.rast.pred))
tic()
for (i in 1:nlyr(br.rast.pred)) {

  # Subset tracks by month.year
  obs <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  cbi.br[[i]] <- cbi(fit = br.rast.pred[[i]],
                     obs = obs,
                     nclass = 0,
                     window.w = "default",
                     res = 100,
                     PEplot = FALSE,
                     rm.duplicate = FALSE,  #after inspecting some of the PE plots, most occurrences have quadratic pattern at lower end of predictive range; corr should be more reflective by including all 0s
                     method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 16 sec


cbi.br <- cbi.br %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil")








### Qatar ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Qatar", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_qa <- sapply(files, rast)
cov_list_qa

names(cov_list_qa) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_qa[[var]]) <- gsub(names(cov_list_qa[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list_qa[[var]] <- terra::resample(cov_list_qa[[var]], cov_list_qa$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list_qa <- map(cov_list_qa, terra::project, 'EPSG:3395')






### Validate model using CBI ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.pred <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.pred) <- my.ind.qa

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
  qa.pred <- predict.gbm(gbm.fit, newdata = vars[,-4], n.trees = gbm.fit$n.trees, na.action = na.pass)

  terra::values(qa.rast.pred[[i]]) <- NA  # initially store all NAs for locs w/o predictions
  terra::values(qa.rast.pred[[i]])[vars$row_id] <- qa.pred  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 3 sec



cbi.qa <- vector("list", nlyr(qa.rast.pred))
tic()
for (i in 1:nlyr(qa.rast.pred)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  cbi.qa[[i]] <- cbi(fit = qa.rast.pred[[i]],
                     obs = obs,
                     nclass = 0,
                     window.w = "default",
                     res = 100,
                     PEplot = FALSE,
                     rm.duplicate = TRUE,
                     method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


cbi.qa <- cbi.qa %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar")






### Summarize Validation Results ###

cbi.fit <- rbind(cbi.br, cbi.qa)

cbi.mean <- cbi.fit %>%
  group_by(Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE))

set.seed(2023)
ggplot(data = cbi.fit, aes(Region, cor)) +
  # geom_boxplot(aes(fill = Region), outlier.color = NA) +
  geom_jitter(aes(fill = Region), pch = 21, height = 0, width = 0.2, alpha = 0.7, size = 5) +
  scale_fill_manual(values = c("#1DB100", "#00A2FF"), guide = "none") +
  geom_point(data = cbi.mean, aes(Region, mean), fill = "black", pch = 21, size = 6) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Continuous Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24))

