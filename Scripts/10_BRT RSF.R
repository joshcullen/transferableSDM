
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
  drop_na(bathym, npp, sst) %>%
  mutate(bathym = abs(bathym),
         log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)
# rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, 5000, 1)



## Define params to optimize for model fitting
hyper_grid <- expand.grid(
  n.trees = c(1000, 5000),
  interaction.depth = c(1, 3),
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
    shrinkage = 0.001,
    bag.fraction = 0.75,
    weights = rsf.pts_10s$wts2
  )

  # add min training error and trees to grid
  # hyper_grid$optimal_trees[i] <- which.min(gbm.fit$valid.error)
  hyper_grid$mean_RMSE[i] <- sqrt(mean(gbm.fit$train.error))
}
toc()  # took 4.5 min to run

hyper_grid %>%
  dplyr::arrange(mean_RMSE)
# best model has 5000 trees, 0.01 shrinkage, and depth of 5


# Fit model w/ optimized params
set.seed(2023)
tic()
gbm.fit <- gbm(
  formula = obs/wts2 ~ bathym + npp + sst,
  distribution = "poisson",
  data = rsf.pts_10s,
  n.trees = 5000,
  interaction.depth = 3,
  bag.fraction = 0.75,
  shrinkage = 0.001,
  weights = rsf.pts_10s$wts2
)
toc()  #took 2.5 min to run

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







###########################
### Export model object ###
###########################

saveRDS(gbm.fit, "Data_products/BRT_model_fit.rds")



