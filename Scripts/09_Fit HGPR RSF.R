
### Fit RSF as Hierarchical Gaussian Process regression ###

library(tidyverse)
library(INLA)
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
### Fit HGPR RSF ###
####################

# Remove rows w/ incomplete observations; log-transform covars
rsf.pts_10s <- rsf.pts_10 %>%
  drop_na(bathym, npp, sst) %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))



# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)




# Add ID as integer for INLA
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(id1 = as.integer(factor(id)))

# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')
y <- rsf.pts_10s$obs / rsf.pts_10s$wts2

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup <- n_distinct(rsf.pts_10s$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)

nbasis <- 5  #number of basis functions for approximating GP
degree <- 2  #degree for defining 1D mesh of GP
alpha <- 2  #for calculating Matern covariance matrix


# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = nbasis),
                                 degree = degree,
                                 boundary = 'free')
}

# Calculate Matern covariance matrix for SPDE (using PC priors)
spde.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  spde.list[[i]] <-  inla.spde2.pcmatern(mesh.list[[i]],
                                         alpha = alpha,
                                         prior.range = c(pcprior[[i]][1], 0.05),
                                         prior.sigma = c(pcprior[[i]][2], 0.05))
}



############################
### Population-level GPs ###
A.list.pop <- vector("list", length(covars))
index.list.pop <- vector("list", length(covars))

for (i in 1:length(covars)) {
  A.list.pop[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                      loc = rsf.pts_10s[[covars[[i]]]]
  )
  index.list.pop[[i]] <-  inla.spde.make.index(paste(covars[i], "pop", sep = "."),
                                               n.spde = spde.list[[i]]$n.spde)
}
############################


################################
### Include random GP slopes ###
A.list.id <- vector("list", length(covars))
index.list.id <- vector("list", length(covars))

for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.list.id[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                     loc = rsf.pts_10s[[covars[[i]]]],
                                     group = rsf.pts_10s$id1,
                                     n.group = ngroup,
  )
  index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
                                              n.spde = spde.list[[i]]$n.spde,
                                              n.group = ngroup)
}
################################


# Create INLA stack of A matrices and other data
st.est <- inla.stack(data = list(y = y),
                     A = c(A.list.pop, A.list.id,
                           1, 1, 1),
                     effects = c(index.list.pop, index.list.id,
                                 list(Intercept = rep(1, nrow(rsf.pts_10s)), id1 = rsf.pts_10s$id1,
                                      log.sst = rsf.pts_10s$log.sst)))

# Define formula for HGPR RSF model
formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +  #fixed terms
  # pop-level terms
    f(log.bathym.pop, model=spde.list[[1]]) +
    f(log.npp.pop, model=spde.list[[2]]) +
    f(log.sst.pop, model=spde.list[[3]]) +
  # id-level terms
    f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
    f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
    f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
    f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))


## Run model
data <-  inla.stack.data(st.est)
set.seed(2023)
tic()
hgpr.fit <- inla(formula, data=data, family="poisson", weights = rsf.pts_10s$wts2,
                 control.predictor = list(A = inla.stack.A(st.est), compute = TRUE),
                 control.fixed = list(  #physiologically-informed SST component
                   mean = c(6.592, -1),
                   prec = c(100, 100)),
                 num.threads = 1:1)  #for greater reproducibility
toc()  #took 12.5 min to run single-threaded

summary(hgpr.fit)




###############################################
### Population-level marginal effects plots ###
###############################################

# Generate matrices for sequences of covariate ranges
A.me.pop <- vector("list", length(covars))
newdat.list <- map(mesh.seq, function(x) {seq(from = x[1], to = x[2], length.out = 500)})
for (i in 1:length(covars)) { #make marginal predictions to viz response by covar
  A.me.pop[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.list[[covars[[i]]]])
}


# Store resulting GP coeffs per covar into a list
pred.coeffs.pop <- hgpr.fit$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Make predictions via linear algebra
pred.vals.pop <- A.me.pop %>%
  map2(.x = ., .y = pred.coeffs.pop,
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
           }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list)) %>%
  mutate(across(mean:x, exp))



# Depth
ggplot() +
  # geom_ribbon(data = pred.vals[[1]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.bathym"), aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

# NPP
ggplot() +
  # geom_ribbon(data = pred.vals[[2]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.npp"), aes(x = x / 1000, y = mean), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = expression(paste("NPP (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

# SST
ggplot() +
  # geom_ribbon(data = pred.vals[[3]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst"), aes(x = x, y = mean), linewidth = 1.5) +
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
A.me.id <- rep(A.me.id, each = max(rsf.pts_10s$id1))

# Store resulting GP coeffs per covar into a list
pred.coeffs.id <- hgpr.fit$summary.random[4:6] %>%
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list.id[[1]]$log.bathym.id.group)) %>%
  map(., ~split(.x, .x$id)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10s$id1)), rep(1:max(rsf.pts_10s$id1), length(covars)), sep = "_"))

# Make predictions via linear algebra
pred.vals.id <- A.me.id %>%
  map2(.x = ., .y = pred.coeffs.id,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  set_names(paste(rep(covars, each = max(rsf.pts_10s$id1)), rep(1:max(rsf.pts_10s$id1), length(covars)), sep = "_")) %>%
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
  lims(x = c(0,800), y = c(0,1.5e6)) +
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
A.pop.sep.20 <- vector("list", length(covars))
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.pop.sep.20[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = rast.sep.20[[covars[[i]]]])
}

# Store resulting GP coeffs per covar into a list
pred.coeffs <- hgpr.fit$summary.random[1:3] %>%
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
rast.pred <- cov_list$bathym
terra::values(rast.pred) <- NA  # initially store all NAs for locs w/o predictions
terra::values(rast.pred)[rast.sep.20$row_id] <- pred.sep.20


rast.pred.df <- as.data.frame(rast.pred, xy = TRUE)
names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)


# Generate predictive surface for GP at pop-level
ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = pred)) +
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







# For ID 181796 in 2020-08
rsf.pts_10s %>%
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
ind <- rsf.pts_10s[rsf.pts_10s$id == 181796,]$id1[1]
pred.coeffs <- hgpr.fit$summary.random[4:6] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list.id[[1]]$log.bathym.id.group)) %>%
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


rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
names(rast.pred.181796.df)[3] <- "pred"
bbox <- ext(rast.pred.181796)

ggplot() +
  geom_raster(data = rast.pred.181796.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_path(data = rsf.pts_10s %>%
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

saveRDS(hgpr.fit, "Data_products/HGPR_model_fit.rds")


