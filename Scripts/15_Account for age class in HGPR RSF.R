
### Fit HGPR RSF at 5 km scale and accounting for life stage ###

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





####################
### Fit HGPR RSF ###
####################

# Remove rows w/ incomplete observations; log-transform covars
rsf.pts_10_5kms <- rsf.pts_10_5km %>%
  drop_na(bathym, npp, sst) %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))



# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
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
                 log.sst = c(12,38)) %>%
  map(log)


hgpr.age <- fit_hgpr(data = rsf.pts_10_5kms, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                     nbasis = 5, degree = 2, alpha = 2, age.class = TRUE, int.strategy = 'auto')  # took 3 min to run




# nbasis <- 5  #number of basis functions for approximating GP
# degree <- 2  #degree for defining 1D mesh of GP
# alpha <- 2  #for calculating Matern covariance matrix
#
#
#
# # Define weighted response variable
# y <- rsf.pts_10_5kms$obs / rsf.pts_10_5kms$wts
#
# # Define 1D mesh per covar
# mesh.list <- vector("list", length(covars))
# for (i in 1:length(covars)) {
#   mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
#                                      length.out = nbasis),
#                                  degree = degree,
#                                  boundary = 'free')
# }
#
# # Calculate Matern covariance matrix for SPDE (using PC priors)
# spde.list <- vector("list", length(covars))
# for (i in 1:length(covars)) {
#   spde.list[[i]] <-  inla.spde2.pcmatern(mesh.list[[i]],
#                                          alpha = alpha,
#                                          prior.range = c(pcprior[[i]][1], 0.05),
#                                          prior.sigma = c(pcprior[[i]][2], 0.95))
# }
#
#
#
# ############################
# ### Life stage-level GPs ###
# A.list.age <- vector("list", length(covars))
# index.list.age <- vector("list", length(covars))
#
# for (i in 1:length(covars)) {
#   A.list.age[[i]] <- inla.spde.make.A(mesh.list[[i]],
#                                       loc = rsf.pts_10_5kms[[covars[[i]]]],
#                                      group = rsf.pts_10_5kms$Age1,
#                                      n.group = ngroup.age
#   )
#   index.list.age[[i]] <-  inla.spde.make.index(paste(covars[i], "age", sep = "."),
#                                               n.spde = spde.list[[i]]$n.spde,
#                                               n.group = ngroup.age)
# }
# ############################
#
#
# ################################
# ### Include random GP slopes ###
# A.list.id <- vector("list", length(covars))
# index.list.id <- vector("list", length(covars))
#
# for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
#   A.list.id[[i]] <- inla.spde.make.A(mesh.list[[i]],
#                                      loc = rsf.pts_10_5kms[[covars[[i]]]],
#                                      group = rsf.pts_10_5kms$id1,
#                                      n.group = ngroup.id
#   )
#   index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
#                                               n.spde = spde.list[[i]]$n.spde,
#                                               n.group = ngroup.id)
# }
# ################################
#
#
# # Create INLA stack of A matrices and other data
# st.est <- inla.stack(data = list(y = y),
#                      A = c(A.list.age, A.list.id,
#                            1, 1, 1),
#                      effects = c(index.list.age, index.list.id,
#                                  list(Intercept = rep(1, nrow(rsf.pts_10_5kms)), id1 = rsf.pts_10_5kms$id1,
#                                       log.sst = rsf.pts_10_5kms$log.sst)))
#
# # Define formula for HGPR RSF model
# formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +  #fixed terms
#   # life stage-level terms
#   f(log.bathym.age, model=spde.list[[1]], group = log.bathym.age.group, control.group = list(model = 'iid')) +
#   f(log.npp.age, model=spde.list[[2]], group = log.npp.age.group, control.group = list(model = 'iid')) +
#   f(log.sst.age, model=spde.list[[3]], group = log.sst.age.group, control.group = list(model = 'iid')) +
#   # id-level terms
#   f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
#   f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
#   f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
#   f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))
#
#
# ## Run model
# stack.data <-  inla.stack.data(st.est)
# set.seed(2023)
# tic()
# hgpr.age <- inla(formula, data=stack.data, family="poisson", weights = rsf.pts_10_5kms$wts,
#                  control.predictor = list(A = inla.stack.A(st.est), compute = TRUE),
#                  control.fixed = list(  #physiologically-informed SST component
#                    mean = list(log.sst = 10*6.592, `I(log.sst^2)` = 10*-1),
#                    prec = list(log.sst = 0.1, `I(log.sst^2)` = 0.1)),
#                  num.threads = 1:1)  #for greater reproducibility
# toc()  #took 14 min to run

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


# Make predictions of linear SST terms
sst.newdata <- data.frame(sst = newdat.list$log.sst,
                          sst.2 = newdat.list$log.sst ^ 2) %>%
  as.matrix()

fixed.sst.coeff <- hgpr.age$summary.fixed$mean[-1]
fixed.sst.pred <- sst.newdata %*% fixed.sst.coeff %>%
  data.frame(pred = .) %>%
  mutate(sst = exp(sst.newdata[,1]))


# Facet of all IDs and covars
ggplot() +
  geom_line(data = pred.vals.age, aes(x = x, y = mean, color = age), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ covar, scales = "free")


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
  # lims(y = c(0,5000)) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ age, scales = "free")


# SST
ggplot() +
  geom_line(data = pred.vals.age %>%
              filter(covar == 'sst'), aes(x = x, y = mean, color = age), linewidth = 1) +
  theme_bw() +
  # lims(y = c(0,2.5e5)) +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ age, scales = "free")

ggplot() +
  geom_line(data = fixed.sst.pred, aes(x = sst, y = pred), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))



#######################################
### ID-level marginal effects plots ###
#######################################

# Generate matrices for covariate raster data (for prediction)
A.me.id <- vector("list", length(covars))
newdat.list <- map(mesh.seq, function(x) {seq(from = x[1], to = x[2], length.out = 500)})
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.me.id[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.list[[covars[[i]]]])
}

# Replicate list elements for mapping over lists
A.me.id <- rep(A.me.id, each = max(rsf.pts_10_5kms$id1))

# Store resulting GP coeffs per covar into a list
pred.coeffs.id <- hgpr.age$summary.random[4:6] %>%
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = rep(1:49, each = 6))) %>%
  map(., ~split(.x, .x$id)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10_5kms$id1)), rep(1:max(rsf.pts_10_5kms$id1), length(covars)), sep = "_"))

# Make predictions via linear algebra
pred.vals.id <- A.me.id %>%
  map2(.x = ., .y = pred.coeffs.id,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10_5kms$id1)), rep(1:max(rsf.pts_10_5kms$id1), length(covars)), sep = "_")) %>%
  bind_rows() %>%
  mutate(bathym = newdat.list$log.bathym,
         npp = newdat.list$log.npp,
         sst = newdat.list$log.sst) %>%
  mutate(across(everything(), exp)) %>%
  pivot_longer(cols = -c(bathym, npp, sst), names_to = "label", values_to = "mean") %>%
  separate(col = label, into = c('covar','id'), sep = "_")


# Wrangle data so that x and y values match up in long format
pred.vals.id <- pred.vals.id %>%
  mutate(x = case_when(str_detect(covar, "log.bathym") ~ bathym,
                       str_detect(covar, "log.npp") ~ npp,
                       str_detect(covar, "log.sst") ~ sst)) %>%
  dplyr::select(-c(bathym, npp, sst)) %>%
  mutate(covar = gsub(pattern = "log.", replacement = "", x = covar))


# Facet of all IDs and covars
ggplot() +
  geom_line(data = pred.vals.id, aes(x = x, y = log(mean), color = factor(id)), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ covar, scales = "free")


# Depth
ggplot() +
  geom_line(data = pred.vals.id %>%
              filter(covar == 'bathym'), aes(x = x, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(x = c(0,800), y = c(0,1e29)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


# NPP
ggplot() +
  geom_line(data = pred.vals.id %>%
              filter(covar == 'npp'), aes(x = x / 1000, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(y = c(0,5000)) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


# SST
ggplot() +
  geom_line(data = pred.vals.id %>%
              filter(covar == 'sst'), aes(x = x, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(y = c(0,2.5e5)) +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))







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


# rast.pred.df <- rbind(as.data.frame(rast.pred.juv, xy = TRUE) %>%
#                         mutate(age = "Juv"),
#                       as.data.frame(rast.pred.adult, xy = TRUE) %>%
#                         mutate(age = "Adult"))
# names(rast.pred.df)[3] <- "pred"
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






# For ID 181796 in 2020-08
rsf.pts_10_5kms %>%
  filter(id == 181796, obs == 1) %>%
  group_by(month.year) %>%
  count()  # most obs are in August

newdat.181796 <- data.frame(log.bathym = as.vector(terra::values(cov_list$bathym)) %>%
                              abs() %>%
                              log(),
                            log.npp = as.vector(terra::values(cov_list$npp$`2020-08-01`)) %>%
                              log(),
                            log.sst = as.vector(terra::values(cov_list$sst$`2020-08-01`)) %>%
                              log()) %>%
  mutate(row_id = 1:nrow(.)) %>%
  drop_na(log.bathym, log.npp, log.sst)
summary(newdat.181796)

A.181796 <- vector("list", length(covars))
for (i in 1:length(covars)) {
  A.181796[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.181796[[covars[[i]]]])
}


# Store resulting GP coeffs per covar into a list
ind <- rsf.pts_10_5kms[rsf.pts_10_5kms$id == 181796,]$id1[1]
pred.coeffs <- hgpr.age$summary.random[4:6] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = rep(1:49, each = 6))) %>%
  map(., ~dplyr::filter(.x, id == ind)) %>%
  map(., pull, mean)


# Make predictions via linear algebra
x181796.pred <- A.181796 %>%
  map2(.x = ., .y = pred.coeffs,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  bind_cols() %>%
  rowSums()  #sum up all predictions across covars


# Define dummy raster for storing predictions
rast.pred.181796 <- cov_list$bathym
terra::values(rast.pred.181796) <- NA  # initially store all NAs for locs w/o predictions
terra::values(rast.pred.181796)[newdat.181796$row_id] <- x181796.pred


# rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
# names(rast.pred.181796.df)[3] <- "pred"
bbox <- ext(rast.pred.181796)

ggplot() +
  geom_spatraster(data = rast.pred.181796) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_path(data = rsf.pts_10_5kms %>%
              filter(id == 181796, month.year == '2020-08-01', obs == 1), aes(x, y),
            color = "chartreuse", linewidth = 0.5) +
  labs(x="",y="", title = "ID 181796: August 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20))




###########################
### Export model object ###
###########################

saveRDS(hgpr.age, "Data_products/HGPR_model_fit_scale_age.rds")


