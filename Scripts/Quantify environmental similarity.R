
### Evaluate environmental similarity between training and testing regions ###
### Using Shape method from {flexsdm} ###

library(tidyverse)
library(terra)
library(flexsdm)
library(tictoc)
library(BRRR)



#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv") |>
  drop_na(bathym, npp, sst)



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






############################################
### Evaluate environmental dissimilarity ###
############################################


### Qatar ###
qa_sim <- vector("list", length = nlyr(cov_list_qa$sst)) |>
  set_names(names(cov_list_qa$sst))

tic()
for (i in seq_along(qa_sim)) {
  qa_sim[[i]] <- extra_eval(training_data = rsf.pts_10,
                            projection_data = rast(c(bathym = cov_list_qa[[1]],
                                                     npp = cov_list_qa[[2]][[i]],
                                                     sst = cov_list_qa[[3]][[i]])),
                            pr_ab = "obs",
                            n_cores = 5,
                            aggreg_factor = 1,
                            metric = "mahalanobis")

  cat("\nLayer", i, "of", length(qa_sim))
}
skrrrahh('khaled2')
toc()  #took 5.5 min w/ 5 cores

plot(qa_sim[[1]], main = "Qatar extrapolation pattern")



### Brazil ###
br_sim <- vector("list", length = nlyr(cov_list_br$sst)) |>
  set_names(names(cov_list_br$sst))

tic()
for (i in seq_along(br_sim)) {
  br_sim[[i]] <- extra_eval(training_data = rsf.pts_10,
                            projection_data = rast(c(bathym = cov_list_br[[1]],
                                                     npp = cov_list_br[[2]][[i]],
                                                     sst = cov_list_br[[3]][[i]])),
                            pr_ab = "obs",
                            n_cores = 20,
                            aggreg_factor = 1,
                            metric = "mahalanobis")

  cat("\nLayer", i, "of", length(br_sim))
}
skrrrahh('khaled2')
toc()  #took 5.75 hrs w/ 20 cores

plot(br_sim[[50]], main = "Brazil extrapolation pattern")
plot(rast(br_sim[1:12]), main = "Brazil extrapolation pattern")



##########################################
### Summarize similarity across months ###
##########################################

# Define index to subset cells that fall w/in neritic habitat (for Brazil only)
br_neritic_ind <- which(values(cov_list_br$bathym) > -200)


# Brazil
br_summ <- br_sim |>
  map(~values(.x) |>
        data.frame()) |>
  map(~summarize(.x,
                 mean = mean(extrapolation, na.rm = TRUE),
                 sd = sd(extrapolation, na.rm = TRUE),
                 mean_neritic = mean(extrapolation[br_neritic_ind], na.rm = TRUE),
                 sd_neritic = sd(extrapolation[br_neritic_ind], na.rm = TRUE))) |>
  bind_rows(.id = "month-year")


# Qatar
qa_summ <- qa_sim |>
  map(~values(.x) |>
        data.frame()) |>
  map(~summarize(.x,
                 mean = mean(extrapolation, na.rm = TRUE),
                 sd = sd(extrapolation, na.rm = TRUE))) |>
  bind_rows(.id = "month-year")




######################
### Export results ###
######################

saveRDS(qa_sim, "Data_products/Qatar_EnvironSim.rds")
saveRDS(br_sim, "Data_products/Brazil_EnvironSim.rds")
