
### Aggregate environmental raster layers to coarser scales ###

library(tidyverse)
library(terra)
library(tidyterra)
library(future)
library(furrr)
library(cmocean)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(patchwork)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10_5km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_10 <- rsf.pts_10_5km %>%
  dplyr::select(-c(bathym, npp, sst))

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



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





################################################
### Adjust scale of environmental covariates ###
################################################

scale.fact <- c(2, 4, 8, 16)
cov_coarse_list <- vector("list", length(scale.fact))
names(cov_coarse_list) <- c("sc.10", "sc.20", "sc.40", "sc.80")

for (i in 1:length(cov_coarse_list)) {
  cov_coarse_list[[i]] <- cov_list %>%
    map(., ~terra::aggregate(.x, fact = scale.fact[i], fun = mean, na.rm = TRUE))
}


# Viz comparison of SST at different scales
sst.sc <- cov_coarse_list %>%
  map_depth(1, pluck, "sst") %>%
  append(cov_list["sst"], .) %>%
  set_names(c(5, 10, 20, 40, 80))


p.sst.5 <- ggplot() +
  geom_spatraster(data = sst.sc$`5`$`2011-10-01`$`2011-10-01`) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf) +
  coord_sf(xlim = c(ext(sst.sc$`5`$`2011-10-01`)[1], ext(sst.sc$`5`$`2011-10-01`)[2]),
           ylim = c(ext(sst.sc$`5`$`2011-10-01`)[3], ext(sst.sc$`5`$`2011-10-01`)[4]),
           expand = FALSE) +
  theme_bw()

p.sst.10 <- ggplot() +
  geom_spatraster(data = sst.sc$`10`$`2011-10-01`) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf) +
  coord_sf(xlim = c(ext(sst.sc$`5`$`2011-10-01`)[1], ext(sst.sc$`5`$`2011-10-01`)[2]),
           ylim = c(ext(sst.sc$`5`$`2011-10-01`)[3], ext(sst.sc$`5`$`2011-10-01`)[4]),
           expand = FALSE) +
  theme_bw()

p.sst.20 <- ggplot() +
  geom_spatraster(data = sst.sc$`20`$`2011-10-01`) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf) +
  coord_sf(xlim = c(ext(sst.sc$`5`$`2011-10-01`)[1], ext(sst.sc$`5`$`2011-10-01`)[2]),
           ylim = c(ext(sst.sc$`5`$`2011-10-01`)[3], ext(sst.sc$`5`$`2011-10-01`)[4]),
           expand = FALSE) +
  theme_bw()

p.sst.40 <- ggplot() +
  geom_spatraster(data = sst.sc$`40`$`2011-10-01`) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf) +
  coord_sf(xlim = c(ext(sst.sc$`5`$`2011-10-01`)[1], ext(sst.sc$`5`$`2011-10-01`)[2]),
           ylim = c(ext(sst.sc$`5`$`2011-10-01`)[3], ext(sst.sc$`5`$`2011-10-01`)[4]),
           expand = FALSE) +
  theme_bw()

p.sst.80 <- ggplot() +
  geom_spatraster(data = sst.sc$`80`$`2011-10-01`) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf) +
  coord_sf(xlim = c(ext(sst.sc$`5`$`2011-10-01`)[1], ext(sst.sc$`5`$`2011-10-01`)[2]),
           ylim = c(ext(sst.sc$`5`$`2011-10-01`)[3], ext(sst.sc$`5`$`2011-10-01`)[4]),
           expand = FALSE) +
  theme_bw()

(p.sst.5 | ((p.sst.10 + p.sst.20) / (p.sst.40 + p.sst.80))) +
  plot_layout(guides = 'collect') & theme(legend.position = "top")





########################################
### Extract environmental covariates ###
########################################

# Extract environ covars by month.year
plan(multisession, workers = availableCores() - 2)
rsf.pts_10_10km <- extract.covars(data = rsf.pts_10, layers = cov_coarse_list$sc.10, dyn_names = c('npp', 'sst'),
                             along = FALSE, ind = "month.year", imputed = FALSE)
#takes 1 min to run on desktop (18 cores)

rsf.pts_10_20km <- extract.covars(data = rsf.pts_10, layers = cov_coarse_list$sc.20, dyn_names = c('npp', 'sst'),
                                  along = FALSE, ind = "month.year", imputed = FALSE)
#takes 9 sec to run on desktop (18 cores)

rsf.pts_10_40km <- extract.covars(data = rsf.pts_10, layers = cov_coarse_list$sc.40, dyn_names = c('npp', 'sst'),
                                  along = FALSE, ind = "month.year", imputed = FALSE)
#takes 7 sec to run on desktop (18 cores)

rsf.pts_10_80km <- extract.covars(data = rsf.pts_10, layers = cov_coarse_list$sc.80, dyn_names = c('npp', 'sst'),
                                  along = FALSE, ind = "month.year", imputed = FALSE)
#takes 7 sec to run on desktop (18 cores)
plan(sequential)






# Viz example of available point covar values by month.year
x181796 <- list(`5km` = rsf.pts_10_5km %>%
  filter(id == 181796) %>%
    mutate(across(id, as.character)),
  `10km` = rsf.pts_10_10km %>%
    filter(id == 181796),
  `20km` = rsf.pts_10_20km %>%
    filter(id == 181796),
  `40km` = rsf.pts_10_40km %>%
    filter(id == 181796),
  `80km` = rsf.pts_10_80km %>%
    filter(id == 181796)
) %>%
  bind_rows(.id = "scale") %>%
  mutate(across(scale, factor, levels = c('5km','10km','20km','40km','80km')))

ggplot() +
  geom_point(data = x181796 %>%
               filter(obs == 0), aes(x, y, color = bathym)) +
  geom_sf(data = gom.sf) +
  geom_point(data = x181796 %>%
               filter(obs == 1), aes(x, y), color = 'red') +
  theme_bw() +
  coord_sf(xlim = c(min(x181796$x), max(x181796$x)),
           ylim = c(min(x181796$y), max(x181796$y))) +
  labs(x="",y="") +
  scale_color_cmocean("Depth (m)", name = 'deep',
                               direction = -1, breaks = c(-250,-200,-150,-100,-50,0)) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ scale)


ggplot() +
  geom_point(data = x181796 %>%
               filter(id == 181796, obs == 0, month.year == "2020-08-01"), aes(x, y, color = sst)) +
  geom_sf(data = gom.sf) +
  geom_point(data = x181796 %>%
               filter(id == 181796, obs == 1, month.year == "2020-08-01"), aes(x, y), color = 'red') +
  scale_color_cmocean("SST (°C)", name = 'thermal') +
  labs(x="",y="") +
  theme_bw() +
  coord_sf(xlim = c(min(x181796$x), max(x181796$x)),
           ylim = c(min(x181796$y), max(x181796$y))) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ scale)







###########################
### Export prepped data ###
###########################

# write_csv(rsf.pts_10_10km, "Processed_data/GoM_Cm_RSFprep_10x_10km.csv")
# write_csv(rsf.pts_10_20km, "Processed_data/GoM_Cm_RSFprep_10x_20km.csv")
# write_csv(rsf.pts_10_40km, "Processed_data/GoM_Cm_RSFprep_10x_40km.csv")
# write_csv(rsf.pts_10_80km, "Processed_data/GoM_Cm_RSFprep_10x_80km.csv")
