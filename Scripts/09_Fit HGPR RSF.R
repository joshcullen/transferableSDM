
### Fit RSF as Hierarchical Gaussian Process regression ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(patchwork)

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
rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)
# rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, 5000, 1)




# Add ID as integer for INLA
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(id1 = as.integer(factor(id))) #%>%
  # arrange(id1)



# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup <- n_distinct(rsf.pts_10s$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)

hgpr.fit <- fit_hgpr(data = rsf.pts_10s, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                         nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto')  #took 3.25 min

summary(hgpr.fit)




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


# Make predictions of linear SST terms
sst.newdata <- data.frame(sst = newdat.list$log.sst,
                          sst.2 = newdat.list$log.sst ^ 2) %>%
  as.matrix()

fixed.sst.coeff <- hgpr.fit$summary.fixed$mean[-1]
# fixed.sst.coeff <- c(10*6.592, 10*-1)

fixed.sst.pred <- sst.newdata %*% fixed.sst.coeff %>%
  data.frame(pred = .) %>%
  mutate(sst = exp(sst.newdata[,1]))


# Depth
p.depth <- ggplot() +
  # geom_ribbon(data = pred.vals[[1]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.bathym"), aes(x = x, y = mean), linewidth = 1, lineend = "round") +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank())

# NPP
p.npp <- ggplot() +
  # geom_ribbon(data = pred.vals[[2]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.npp"), aes(x = x / 1000, y = mean), linewidth = 1, lineend = "round") +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = expression(paste("NPP (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank())

# SST
ggplot() +
  # geom_ribbon(data = pred.vals[[3]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst"), aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

ggplot() +
  geom_line(data = fixed.sst.pred, aes(x = sst, y = pred), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

# Joint plot for SST
p.sst <- ggplot() +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst"), aes(x = x, y = mean + fixed.sst.pred$pred), linewidth = 1,
            lineend = "round") +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank())








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
  map(., ~mutate(.x, id = rep(1:49, each = 6))) %>%
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
p.pred_map <- ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="") +
  theme_void() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) #+
  # theme(axis.text = element_text(size = 12),
  #       plot.title = element_text(size = 26, face = "bold"),
  #       legend.title = element_text(size = 18),
  #       legend.text = element_text(size = 16)) +
  # guides(fill = guide_colourbar(barwidth = 2, barheight = 20))



## Create composite plot w/ marginal effects and mapped prediction
plot_spacer() + plot_spacer() + p.depth + p.npp + p.sst + p.pred_map +
  plot_layout(ncol = 2, nrow = 3, heights = c(0.5, 5, 5)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag.position = c(0.08, 1),
        plot.tag = element_text(size = 18, hjust = 0, vjust = -0.4, face = 'bold'))

ggsave("Tables_Figs/Figure 6.png", width = 7, height = 5, units = "in", dpi = 400)




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


