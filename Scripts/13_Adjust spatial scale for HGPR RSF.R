
### Fit HGLM RSF at different spatial scales ###

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

rsf.pts_10_5km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_10_10km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_10km.csv")
rsf.pts_10_20km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_20km.csv")
rsf.pts_10_40km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_40km.csv")
# rsf.pts_10_80km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_80km.csv")

rsf.list <- list(
  sc.5 = rsf.pts_10_5km,
  sc.10 = rsf.pts_10_10km,
  sc.20 = rsf.pts_10_20km,
  sc.40 = rsf.pts_10_40km#,
  # sc.80 = rsf.pts_10_80km
)

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10_5km)
summary(rsf.pts_10_5km)



#####################
### Fit HGLM RSFs ###
#####################

# Remove rows w/ incomplete observations; log-transform covars
rsf.list2 <- rsf.list %>%
  map(., ~{.x %>%
      drop_na(bathym, npp, sst) %>%
      mutate(log.bathym = log(abs(bathym)),
             log.npp = log(npp),
             log.sst = log(sst))}
      )


# Now explore transformed distributions
rsf.list2 %>%
  bind_rows(.id = "scale") %>%
  mutate(across(scale, factor, levels = c('sc.5','sc.10','sc.20','sc.40'))) %>%
  pivot_longer(cols = c(log.bathym, log.npp, log.sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(scale ~ covar, scales = "free")


### Down-weighted Poisson regression

# Define pixel area for each spatial scale (in m^2)
Area.list <- list(sc.5 = 4759.836 ^ 2,
                  sc.10 = 9532.72 ^ 2,
                  sc.20 = 19065.44 ^ 2,
                  sc.40 = 38130.88 ^ 2#,
                  # sc.80 = 76261.76 ^ 2
                  )
rsf.list2 <- rsf.list2 %>%
  map2(.x = ., .y = Area.list,
       ~{.x %>%
           mutate(wts = case_when(obs == 0 ~ .y / sum(obs == 0),
                                  obs == 1 ~ 1e-6))
           # mutate(wts = case_when(obs == 0 ~ 5000,
           #                        obs == 1 ~ 1))
         }
  )



# Add ID as integer for INLA
rsf.list2 <- rsf.list2 %>%
  map(.x = .,
       ~{.x %>%
           mutate(id1 = as.integer(factor(id)),
                  id2 = id1,
                  id3 = id1,
                  id4 = id1,
                  id5 = id1,
                  id6 = id1,
                  id7 = id1) #%>%
           # arrange(id1)
         }
  )

# # create vector of ID values
# id.vals <- unique(rsf.list2$sc.5$id1)
#
# RSF.formula <- obs/wts ~ log.bathym + I(log.bathym ^ 2) + log.npp + I(log.npp ^ 2) + log.sst + I(log.sst ^ 2) +
#   f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
#   f(id2, log.bathym, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id3, I(log.bathym ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id4, log.npp, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id5, I(log.npp ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id6, log.sst, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id7, I(log.sst ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))


# plan(multisession, workers = availableCores() - 2)
# tic()
# hglm.fit <- future_map(rsf.list2,
#                        ~fit_hglm(data = .x),
#                        .options = furrr_options(seed = TRUE))
# toc()
# plan(sequential)  #took 23 min to run


# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup <- n_distinct(rsf.list2$sc.5$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)


hgpr.mod_5km <- fit_hgpr(data = rsf.list2$sc.5, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                         nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto')
#took 4 min
summary(hgpr.mod_5km)

hgpr.mod_10km <- fit_hgpr(data = rsf.list2$sc.10, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                          nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto')
#took 5 min to run
summary(hgpr.mod_10km)

hgpr.mod_20km <- fit_hgpr(data = rsf.list2$sc.20, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                          nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'eb')
#took 4.3 min to run
summary(hgpr.mod_20km)

hgpr.mod_40km <- fit_hgpr(data = rsf.list2$sc.40, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                          nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto')
#took 3.5 min to run
summary(hgpr.mod_40km)

# hgpr.mod_80km <- fit_hgpr(data = rsf.list2$sc.80, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
#                           nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'eb')
# #took x min to run
# summary(hgpr.mod_80km)


# summary(hglm.fit$sc.5)
# summary(hglm.fit$sc.10)
# summary(hglm.fit$sc.20)
# summary(hglm.fit$sc.40)


# Store all model results in list
hgpr.fit <- list(sc.5 = hgpr.mod_5km,
                 sc.10 = hgpr.mod_10km,
                 sc.20 = hgpr.mod_20km,
                 sc.40 = hgpr.mod_40km)




###############################################
### Population-level marginal effects plots ###
###############################################

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
pred.coeffs.pop <- vector("list", length = length(hgpr.fit))
pred.vals.pop <- vector("list", length = length(hgpr.fit)) %>%
  set_names(names(hgpr.fit))

for (i in 1:length(pred.coeffs.pop)) {
pred.coeffs.pop[[i]] <- hgpr.fit[[i]]$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Make predictions via linear algebra
pred.vals.pop[[i]] <- A.me.pop %>%
  map2(.x = ., .y = pred.coeffs.pop[[i]],
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
       }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list)) %>%
  mutate(across(mean:x, exp))
}
pred.vals.pop <- pred.vals.pop %>%
  bind_rows(.id = "scale") %>%
  mutate(across(scale, \(x) factor(x, levels = c('sc.5','sc.10','sc.20','sc.40'))))


# Make predictions of linear SST terms
sst.newdata <- data.frame(sst = newdat.list$log.sst,
                          sst.2 = newdat.list$log.sst ^ 2) %>%
  as.matrix()

fixed.sst.pred <- vector("list", length = length(hgpr.fit)) %>%
  set_names(names(hgpr.fit))

for (i in 1:length(fixed.sst.pred)) {
  fixed.sst.coeff <- hgpr.fit[[i]]$summary.fixed$mean[-1]

  fixed.sst.pred[[i]] <- sst.newdata %*% fixed.sst.coeff %>%
    data.frame(pred = .) %>%
    mutate(sst = exp(sst.newdata[,1]))
}
fixed.sst.pred <- fixed.sst.pred %>%
  bind_rows(.id = "scale") %>%
  mutate(across(scale, factor, levels = c('sc.5','sc.10','sc.20','sc.40')))


# Depth
ggplot() +
  # geom_ribbon(data = pred.vals[[1]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.bathym"), aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ scale, scales = "free")

# NPP
ggplot() +
  # geom_ribbon(data = pred.vals[[2]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.npp"), aes(x = x / 1000, y = mean), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = expression(paste("NPP (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ scale, scales = "free")

# SST
ggplot() +
  # geom_ribbon(data = pred.vals[[3]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst"), aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ scale, scales = "free")

ggplot() +
  geom_line(data = fixed.sst.pred, aes(x = sst, y = pred), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ scale, scales = "free")








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

## Create coarsened raster layers
scale.fact <- c(1, 2, 4, 8)
cov_coarse_list <- vector("list", length(scale.fact))
names(cov_coarse_list) <- c("sc.5", "sc.10", "sc.20", "sc.40")

for (i in 1:length(cov_coarse_list)) {
  cov_coarse_list[[i]] <- cov_list %>%
    map(., ~terra::aggregate(.x, fact = scale.fact[i], fun = mean, na.rm = TRUE))
}



####################################
### Generate predictive surfaces ###
####################################
rast.pred <- vector("list", length = length(cov_coarse_list)) %>%
  set_names(names(cov_coarse_list))

for (i in 1:length(rast.pred)) {
  rast.sep.20 <- data.frame(log.bathym = log(abs(terra::values(cov_coarse_list[[i]]$bathym))) %>%
                              as.vector(),
                            log.npp = log(terra::values(cov_coarse_list[[i]]$npp$`2020-09-01`)) %>%
                              as.vector(),
                            log.sst = log(terra::values(cov_coarse_list[[i]]$sst$`2020-09-01`)) %>%
                              as.vector()) %>%
    mutate(row_id = 1:nrow(.)) %>%
    drop_na(log.bathym, log.npp, log.sst)


  # Generate matrices for covariate raster data (for prediction)
  A.pop.sep.20 <- vector("list", length(covars))
  for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.pop.sep.20[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = rast.sep.20[[covars[[j]]]])
  }

  # Store resulting GP coeffs per covar into a list
  pred.coeffs <- hgpr.fit[[i]]$summary.random[1:3] %>%
    map(., ~pull(.x, mean))

  # Make predictions via linear algebra
  pred.sep.20 <- A.pop.sep.20 %>%
    map2(.x = ., .y = pred.coeffs,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols() %>%
    rowSums()  #sum up all predictions across covars


  # Define dummy raster for storing predictions
  rast.pred[[i]] <- cov_coarse_list[[i]]$bathym
  terra::values(rast.pred[[i]]) <- NA  # initially store all NAs for locs w/o predictions
  terra::values(rast.pred[[i]])[rast.sep.20$row_id] <- pred.sep.20
  }

bbox <- ext(rast.pred[[1]])


# Generate predictive surface for GP at pop-level
pred.5 <- ggplot() +
  geom_spatraster(data = rast.pred$sc.5) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "5 km") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)

pred.10 <- ggplot() +
  geom_spatraster(data = rast.pred$sc.10) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "10 km") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)

pred.20 <- ggplot() +
  geom_spatraster(data = rast.pred$sc.20) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "20 km") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)

pred.40 <- ggplot() +
  geom_spatraster(data = rast.pred$sc.40) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "40 km") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE)


((pred.5 + pred.10) / (pred.20 + pred.40))







############################
### Export model objects ###
############################

saveRDS(hgpr.fit, "Data_products/HGPR_model_fit_scale.rds")
