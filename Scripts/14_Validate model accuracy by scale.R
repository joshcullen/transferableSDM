
### Compare model transferability by scale ###

library(tidyverse)
library(INLA)
library(terra)
library(sf)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(tidyterra)
library(patchwork)
library(future)
library(furrr)
library(wesanderson)

source('Scripts/helper functions.R')


##########################
### Load fitted models ###
##########################

hgpr.fit <- readRDS("Data_products/HGPR_model_fit_scale.rds")




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

# Create coarsened raster layers
scale.fact <- c(1, 2, 4, 8)
cov_coarse_br <- vector("list", length(scale.fact))
names(cov_coarse_br) <- c("sc.5", "sc.10", "sc.20", "sc.40")

for (i in 1:length(cov_coarse_br)) {
  cov_coarse_br[[i]] <- cov_list_br %>%
    map(., ~terra::aggregate(.x, fact = scale.fact[i], fun = mean, na.rm = TRUE))
}






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

# Create coarsened raster layers
scale.fact <- c(1, 2, 4, 8)
cov_coarse_qa <- vector("list", length(scale.fact))
names(cov_coarse_qa) <- c("sc.5", "sc.10", "sc.20", "sc.40")

for (i in 1:length(cov_coarse_qa)) {
  cov_coarse_qa[[i]] <- cov_list_qa %>%
    map(., ~terra::aggregate(.x, fact = scale.fact[i], fun = mean, na.rm = TRUE))
}






######################################
### Validate HGPR by spatial scale ###
######################################


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







### Brazil ###

tic()
br.rast.hgpr <- map2(.x = cov_coarse_br, .y = hgpr.fit,
                            ~predict.hgpr(cov_list = .x, model_fit = .y, covars = covars,
                                          mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE)
                     )
skrrrahh('khaled2')
toc()  #took 1.5 min to run


# Normalize predictions on 0-1 scale
br.rast.hgpr2 <- map(br.rast.hgpr, normalize)



# Assess model performance via Continuous Boyce Index
boyce.br.full <- boyce.br.sub <- vector("list", length(br.rast.hgpr2)) %>%
  set_names(names(cov_coarse_br))

my.ind.br <- names(cov_coarse_br$sc.5$npp)

tic()
for (j in 1:length(br.rast.hgpr2)) {
  for (i in 1:nlyr(br.rast.hgpr2[[j]])) {

    # Subset tracks by month.year
    obs_full <- dat.br %>%
      filter(month.year == my.ind.br[i]) %>%
      dplyr::select(x, y)

    obs_sub <- dat.br %>%
      filter(month.year == my.ind.br[i], x < -3800000) %>%
      dplyr::select(x, y)

    boyce.br.full[[j]][[i]] <- boyce(fit = br.rast.hgpr2[[j]][[i]],
                                     obs = obs_full,
                                     nbins = 10,
                                     bin.method = "seq",
                                     PEplot = FALSE,
                                     rm.duplicate = FALSE,
                                     method = "spearman")

    boyce.br.sub[[j]][[i]] <- boyce(fit = br.rast.hgpr2[[j]][[i]],
                                    obs = obs_sub,
                                    nbins = 10,
                                    bin.method = "seq",
                                    PEplot = FALSE,
                                    rm.duplicate = FALSE,
                                    method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 10 sec



perc.use.br.full <- boyce.br.full %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
      ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      # t() %>%
      data.frame()}
      )

# check fewest bins that contain >=90% of all obs
map(perc.use.br.full, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()}
  )  #slightly increases w/ increasing scale

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full %>%
  map(~mutate(.x, bin = factor(10:1, levels = 10:1))) %>%
  bind_rows(.id = "scale") %>%
  mutate(scale = factor(scale, levels = c('sc.5','sc.10','sc.20','sc.40'))) %>%
  pivot_longer(cols = -c(scale, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Brazil_all") +
  facet_wrap(~scale)



perc.use.br.sub <- boyce.br.sub %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
  ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      # t() %>%
      data.frame()}
  )

# check fewest bins that contain >=90% of all obs
map(perc.use.br.sub, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
    mean()}
)  #1.5 bins for each scale; shows impact of FDN locs

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub %>%
  map(~mutate(.x, bin = factor(10:1, levels = 10:1))) %>%
  bind_rows(.id = "scale") %>%
  mutate(scale = factor(scale, levels = c('sc.5','sc.10','sc.20','sc.40'))) %>%
  pivot_longer(cols = -c(scale, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Brazil_sub") +
  facet_wrap(~scale)



boyce.br.full2 <- boyce.br.full %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
      ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Brazil_all")}
  ) %>%
  bind_rows(.id = "scale")

boyce.br.sub2 <- boyce.br.sub %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
  ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Brazil_main")}
  ) %>%
  bind_rows(.id = "scale")





### Qatar ###

tic()
qa.rast.hgpr <- map2(.x = cov_coarse_qa, .y = hgpr.fit,
                     ~predict.hgpr(cov_list = .x, model_fit = .y, covars = covars,
                                   mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE)
)
skrrrahh('khaled2')
toc()  #took 8 sec to run


# Normalize predictions on 0-1 scale
qa.rast.hgpr2 <- map(qa.rast.hgpr, normalize)



# Assess model performance via Continuous Boyce Index
boyce.qa <- vector("list", length(qa.rast.hgpr2)) %>%
  set_names(names(cov_coarse_qa))

my.ind.qa <- names(cov_coarse_qa$sc.5$npp)

tic()
for (j in 1:length(qa.rast.hgpr2)) {
  for (i in 1:nlyr(qa.rast.hgpr2[[j]])) {

    # Subset tracks by month.year
    obs <- dat.qa %>%
      filter(month.year == my.ind.qa[i]) %>%
      dplyr::select(x, y)

    boyce.qa[[j]][[i]] <- boyce(fit = qa.rast.hgpr2[[j]][[i]],
                                obs = obs,
                                nbins = 10,
                                bin.method = "seq",
                                PEplot = FALSE,
                                rm.duplicate = FALSE,
                                method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 1 sec



perc.use.qa <- boyce.qa %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
  ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      data.frame()}
  )

# check fewest bins that contain >=90% of all obs
map(perc.use.qa, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
    mean()}
)  #slightly increases w/ increasing scale

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa %>%
  map(~mutate(.x, bin = factor(10:1, levels = 10:1))) %>%
  bind_rows(.id = "scale") %>%
  mutate(scale = factor(scale, levels = c('sc.5','sc.10','sc.20','sc.40'))) %>%
  pivot_longer(cols = -c(scale, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Qatar") +
  facet_wrap(~scale)




boyce.qa2 <- boyce.qa %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
  ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Qatar")}
  ) %>%
  bind_rows(.id = "scale")








####################################
### Summarize Validation Results ###
####################################

boyce.fit <- rbind(boyce.br.full2, boyce.br.sub2, boyce.qa2) %>%
  mutate(across(scale, \(x) factor(x, levels = c('sc.5','sc.10','sc.20','sc.40'))))
levels(boyce.fit$scale) <- c(5, 10, 20, 40)

boyce.mean <- boyce.fit %>%
  group_by(scale, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = scale), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(scale, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.mean, aes(x = Region, y = mean, group = scale),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_viridis_d(option = "mako", guide = "none", begin = 0.5, end = 0.95, direction = -1) +
  scale_fill_viridis_d("Scale (km)", option = "mako", begin = 0.5, end = 0.95, direction = -1) +
  scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))

# ggsave("Tables_Figs/Figure 8.png", width = 7, height = 5, units = "in", dpi = 400)




# Create example prediction maps per method

bbox <- ext(qa.rast.hgpr2$sc.5)

# Break rasters into bins used for Boyce Index
qa.rast.hgpr.d <- map(qa.rast.hgpr2,
                      ~classify(.x, seq(0, 1, by = 0.1)))
qa.rast.hgpr.d <- map(qa.rast.hgpr.d, function(x) x + 1)
qa.rast.hgpr.df <- map(qa.rast.hgpr.d,
                       ~{as.data.frame(.x, xy = TRUE) %>%
                           mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1)))
                         }) %>%
  bind_rows(.id = 'Scale') %>%
  mutate(across(Scale, \(x) factor(x, levels = c("sc.5","sc.10","sc.20","sc.40"))))
levels(qa.rast.hgpr.df$Scale) <- c('5 km', '10 km', '20 km', '40 km')


# time-matched observations
tmp.qa <- dat.qa %>%
  filter(month.year == "2014-03-01") %>%
  st_as_sf(., coords = c('x','y'), crs = 3395, remove = FALSE) %>%
  .[lengths(st_intersects(., qa.sf)) == 0,]  #remove all pts overlapping land

ggplot() +
  geom_raster(data = qa.rast.hgpr.df, aes(x, y, fill = `2014-03-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = qa.sf) +
  geom_point(data = tmp.qa, aes(x, y), fill = "chartreuse", alpha = 0.4, size = 1,
             shape = 21, stroke = 0.25) +
  labs(x="",y="", title = "March 2014") +
  geom_text(data = data.frame(x = c(5700000,0,0,0),
                              y = c(2870000,0,0,0),
                              Scale = factor(c('5 km','10 km','20 km','40 km'),
                                              levels = c('5 km','10 km','20 km','40 km')),
                              lab = c("Qatar","","","")), aes(x, y, label = lab),
            size = 5, fontface = "italic") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        strip.text = element_text(size = 10, face = "bold", hjust = 0)) +
  facet_wrap(~ Scale)

# ggsave("Tables_Figs/Figure 7.png", width = 4, height = 5, units = "in", dpi = 400)


# p.5km.qa <- ggplot() +
#   geom_spatraster(data = qa.rast.hgpr$sc.5, aes(fill = `2014-03-01`)) +
#   scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
#   geom_sf(data = qa.sf) +
#   labs(x="",y="", title = "5 km") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.10km.qa <- ggplot() +
#   geom_spatraster(data = qa.rast.hgpr$sc.10, aes(fill = `2014-03-01`)) +
#   scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
#   geom_sf(data = qa.sf) +
#   labs(x="",y="", title = "10 km") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.20km.qa <- ggplot() +
#   geom_spatraster(data = qa.rast.hgpr$sc.20, aes(fill = `2014-03-01`)) +
#   scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
#   geom_sf(data = qa.sf) +
#   labs(x="",y="", title = "20 km") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(plot.title = element_text(size = 16, face = "bold"),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#
#
# p.40km.qa <- ggplot() +
#   geom_spatraster(data = qa.rast.hgpr$sc.40, aes(fill = `2014-03-01`)) +
#   scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
#   geom_sf(data = qa.sf) +
#   # geom_point(data = tmp.pts, aes(x, y), color = "blue", alpha = 0.7, size = 1) +
#   labs(x="",y="", title = "40 km") +
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
# p.5km.qa + p.10km.qa + p.20km.qa + p.40km.qa +
#   plot_layout(ncol = 2)

# ggsave("Tables_Figs/Figure 7.png", width = 6, height = 5, units = "in", dpi = 400)







# perc.use.br.full <- perc.use.br.full %>%
#   bind_rows(.id = "scale") %>%
#   mutate(Region = "Brazil_all")
# perc.use.br.sub <- perc.use.br.sub %>%
#   bind_rows(.id = "scale") %>%
#   mutate(Region = "Brazil_sub")
# perc.use.qa <- perc.use.qa %>%
#   bind_rows(.id = "scale") %>%
#   mutate(Region = "Qatar")

# cum.perc <- list(Brazil_all = perc.use.br.full,
#                  Brazil_sub = perc.use.br.sub,
#                  Qatar = perc.use.qa)
#
# cum.perc.mean <- cum.perc %>%
#   map_depth(., 2, ~apply(.x, 2, function(x) which(x >= 0.9)[1])) %>%
#   map_depth(., 2, ~{data.frame(mean = mean(.x),
#                       sd = sd(.x)
#   )}
#   ) %>%
#   map(., bind_rows, .id = "scale") %>%
#   bind_rows(.id = "Region") %>%
#   mutate(across(scale, factor, levels = c('sc.5','sc.10','sc.20','sc.40')))
#
#
#
# ggplot(data = cum.perc.mean, aes(Region, mean)) +
#   geom_linerange(aes(ymin = (mean - sd), ymax = (mean + sd), color = scale),
#                  position = position_dodge(width = 0.75)) +
#   geom_point(aes(group = scale, color = scale),
#              size = 6, position = position_dodge(width = 0.75)) +
#   # lims(y = c(0,10)) +
#   labs(x="", y = "Avg # of bins accounting for 90% of obs") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 22))





###############################################################
### Create composite plot of marginal effects across models ###
###############################################################


# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),
                                 degree = 2,
                                 boundary = 'free')
}




# Define predictive seq of covars
vars <- data.frame(log.bathym = seq(mesh.seq$log.bathym[1], mesh.seq$log.bathym[2], length.out = 500),
                   log.npp = seq(mesh.seq$log.npp[1], mesh.seq$log.npp[2], length.out = 500),
                   log.sst = seq(mesh.seq$log.sst[1], mesh.seq$log.sst[2], length.out = 500))

vars2 <- data.frame(Intercept = 1,
                    log.sst = seq(mesh.seq$log.sst[1], mesh.seq$log.sst[2], length.out = 500),
                    log.sst2 = seq(mesh.seq$log.sst[1], mesh.seq$log.sst[2], length.out = 500) ^ 2)


# Generate A matrices for prediction
A.mat <- vector("list", length(covars))
for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = vars[[covars[[j]]]])
}
names(A.mat) <- covars




### Predict across models w/ different scales ###

marg.eff.list <- vector("list", length = length(hgpr.fit)) %>%
  set_names(names(hgpr.fit))


for (i in seq_along(hgpr.fit)) {
  # Define coeff values from HGPR
  coeff1 <- hgpr.fit[[i]]$summary.random[1:3] %>%
    map(., ~pull(.x, mean))

  # Define coeff values of fixed terms from HGPR
  coeff2 <- hgpr.fit[[i]]$summary.fixed$mean

  # Make predictions on intensity of use from model for GP terms
  hgpr.pred <- A.mat %>%
    map2(.x = ., .y = coeff1,
         ~{.x %*% .y %>%
             as.vector()}
    ) %>%
    bind_cols() %>%
    pivot_longer(cols = everything(), names_to = "covar", values_to = "mean") %>%
    arrange(covar) %>%
    mutate(x = c(vars$log.bathym, vars$log.npp, vars$log.sst)) %>%
    mutate(across(c(mean, x), \(z) exp(z))) %>%
    mutate(covar = case_when(covar == "log.bathym" ~ "depth",
                             covar == "log.npp" ~ "npp",
                             covar == "log.sst" ~ "sst")) %>%
    filter(!(covar == "depth" & x > 300))  #to constrain plot to only show depth up to 300 m

  # Make predictions using linear terms
  hgpr.pred2 <- as.matrix(vars2) %*% coeff2 %>%
    exp() %>%
    as.vector()


  # Add linear SST predictions to GP SST predictions
  hgpr.pred[hgpr.pred$covar == "sst",]$mean <- hgpr.pred[hgpr.pred$covar == "sst",]$mean +
    hgpr.pred2


  marg.eff.list[[i]] <- hgpr.pred
}



marg.eff <- bind_rows(marg.eff.list, .id = "Scale") %>%
  mutate(Scale = case_when(Scale == "sc.5" ~ "5 km",
                           Scale == "sc.10" ~ "10 km",
                           Scale == "sc.20" ~ "20 km",
                           Scale == "sc.40" ~ "40 km"),
         x = case_when(covar == "npp" ~ x / 1000,
                       TRUE ~ x),
         covar = case_when(covar == "depth" ~ "Depth (m)",
                           covar == "npp" ~ "NPP (g C m^-2 d^-1)",
                           covar == "sst" ~ "SST (°C)")) %>%
  mutate(Scale = factor(Scale, levels = unique(Scale)))





# Plot pop-level marginal effects
ggplot() +
  geom_line(data = marg.eff, aes(x = x, y = mean, color = Scale), linewidth = 1) +
  scale_color_viridis_d(option = "mako", begin = 0.5, end = 0.95, direction = -1) +
  theme_bw() +
  labs(x = "", y = "Relative Intensity of Use") +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold")
  ) +
  facet_wrap(Scale ~ covar, ncol = 3, scales = "free", strip.position = "bottom")

# ggsave("Tables_Figs/Figure S4.png", width = 11, height = 9, units = "in", dpi = 400)




### Export Boyce Index results ###

# write_csv(boyce.fit, "Data_products/boyce_scale_results.csv")
