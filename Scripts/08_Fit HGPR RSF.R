
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
A <- 2345557  #area of study region (Gulf of Mexico) in km^2; from region used to generate 'available' pts
rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)




# Add ID as integer for INLA
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(id1 = as.integer(factor(id)))


# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores range and SD, respectively
ngroup <- n_distinct(rsf.pts_10s$id1)  #number of IDs
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)  #range of covars on which to estimate model (on log-scale)


# Fit hybrid HGPR
hgpr.fit_hybrid <- fit_hgpr(data = rsf.pts_10s, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                            nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto',
                            method = "hybrid")
#took 1 min

summary(hgpr.fit_hybrid)


# Fit correlative HGPR
hgpr.fit_corr <- fit_hgpr(data = rsf.pts_10s, covars = covars, pcprior = pcprior, mesh.seq = mesh.seq,
                          nbasis = 5, degree = 2, alpha = 2, age.class = FALSE, int.strategy = 'auto',
                          method = "corr")
#took 1 min

summary(hgpr.fit_corr)




###############################################
### Population-level marginal effects plots ###
###############################################

# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),  #5 basis functions for 1D GP
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
pred.coeffs.pop_hybrid <- hgpr.fit_hybrid$summary.random[1:3] %>%
  map(., ~pull(.x, mean))
pred.coeffs.pop_corr <- hgpr.fit_corr$summary.random[1:3] %>%
  map(., ~pull(.x, mean))

# Make predictions via linear algebra
pred.vals.pop_hybrid <- A.me.pop %>%
  map2(.x = ., .y = pred.coeffs.pop_hybrid,
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
           }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list),
         method = "hybrid") %>%
  mutate(across(mean:x, exp))  #transform back to response/original scale

pred.vals.pop_corr <- A.me.pop %>%
  map2(.x = ., .y = pred.coeffs.pop_corr,
       ~{.x %*% .y %>%
           as.vector() %>%
           data.frame(mean = .)
       }
  ) %>%
  set_names(covars) %>%
  bind_rows(.id = 'covar') %>%
  mutate(x = unlist(newdat.list),
         method = "corr") %>%
  mutate(across(mean:x, exp))  #transform back to response/original scale


# Bind both model predictions together
pred.vals.pop <- rbind(pred.vals.pop_corr,
                       pred.vals.pop_hybrid)


### Make predictions of linear SST terms

# Subset only data range for SST vals
sst.newdata <- data.frame(sst = newdat.list$log.sst,
                          sst.2 = newdat.list$log.sst ^ 2) %>%
  as.matrix()

# Store fixed effect coeffs from hybrid HGPR (w/o intercept)
fixed.sst.coeff <- hgpr.fit_hybrid$summary.fixed$mean[-1]

# Make predictions via linear algebra
fixed.sst.pred <- sst.newdata %*% fixed.sst.coeff %>%
  data.frame(pred = .) %>%
  mutate(sst = exp(sst.newdata[,1]),
         pred = exp(pred))


### Viz marginal effects for corr HGPR model ###

# Depth
p.depth <- ggplot() +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.bathym",
                     method == 'corr'), aes(x = x, y = mean), linewidth = 1, lineend = "round") +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank()) #+
  # facet_wrap(~ method)

# NPP
p.npp <- ggplot() +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.npp",
                     method == 'corr'), aes(x = x / 1000, y = mean), linewidth = 1, lineend = "round") +
  theme_bw() +
  labs(x = expression(paste("NPP (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank()) #+
  # facet_wrap(~ method)

# SST
ggplot() +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst"), aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ method, scales = "free_y")

ggplot() +
  geom_line(data = fixed.sst.pred, aes(x = sst, y = pred), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

p.sst <- ggplot() +
  geom_line(data = pred.vals.pop %>%
              filter(covar == "log.sst",
                     method == 'corr'),
            aes(x = x, y = mean),
                linewidth = 1, lineend = "round") +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(#axis.title = element_text(size = 30),
        # axis.text = element_text(size = 24),
        panel.grid = element_blank()) #+
    # facet_wrap(~ method)











#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
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

# stores raster values from Setpember 2020 into data.frame
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
pred.coeffs <- hgpr.fit_corr$summary.random[1:3] %>%
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

# Convert from raster to data.frame (for ggplot2)
rast.pred.df <- as.data.frame(rast.pred, xy = TRUE)
names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)


### Generate predictive surface for GP at pop-level

# select points co-occurring during September 2020
tmp.pts <- rsf.pts_10s %>%
  filter(obs == 1, month.year == "2020-09-01")

p.pred_map <- ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_point(data = tmp.pts, aes(x, y), fill = "chartreuse", alpha = 0.8, size = 2,
             shape = 21, stroke = 0.25) +
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

# ggsave("Tables_Figs/Figure 4.png", width = 7, height = 5, units = "in", dpi = 400)







###########################
### Export model object ###
###########################

saveRDS(hgpr.fit_corr, "Data_products/HGPR_corr_model_fit.rds")
saveRDS(hgpr.fit_hybrid, "Data_products/HGPR_hybrid_model_fit.rds")


