
### Prepare tracks for RSF ###

library(tidyverse)
library(lubridate)
library(bayesmove)
library(terra)
library(future)
library(furrr)
library(sf)
library(sfarrow)
library(tictoc)
library(amt)
library(tidyterra)
library(cmocean)
library(patchwork)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

dat <- read_csv("Processed_data/GoM_Cm_Tracks_behav.csv")
gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(dat)
summary(dat)


dat <- dat %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))

dat <- prep_data(dat = dat, coord.names = c('x','y'), id = "id")

# Define bounding box of study extent for population
bbox <- dat %>%
  st_as_sf(coords = c('x','y'), crs = 3395) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(5000)  #buffer each side by 5 km

ggplot() +
  geom_sf(data = gom.sf) +
  geom_point(data = dat, aes(x, y, color = behav), alpha = 0.6) +
  geom_sf(data = bbox, color = "red", fill = "transparent", linewidth = 1) +
  scale_color_viridis_d("State", end = 0.98) +
  labs(x="",y="") +
  coord_sf(expand = FALSE) +
  theme_bw()
# ggsave("../../Conference Presentations/SERSTM 2023/behav_gom_map.png", width = 8, height = 6,
#        units = "in", dpi = 400)


#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
# files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
# files <- files[!grepl(pattern = "Kd490", files)]  #remove Kd490 datasets
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



## Map example of environmental covariates

# Define map extent
bbox_p <- ext(cov_list$bathym)

p.gom.bathym <- ggplot() +
  geom_spatraster(data = cov_list$bathym, aes(fill = `GoM bathymetry`)) +
  scale_fill_cmocean("Depth (m)", name = "deep", direction = -1, breaks = c(0, -2000, -4000)) +
  geom_sf(data = gom.sf, linewidth = 0.25) +
  coord_sf(xlim = c(bbox_p[1], bbox_p[2]),
           ylim = c(bbox_p[3], bbox_p[4]),
           expand = FALSE,
           label_axes = "----") +
  theme_bw() +
  theme(legend.position = "top")

p.gom.npp <- ggplot() +
  geom_spatraster(data = cov_list$npp / 1000, aes(fill = `2019-10-01`)) +
  scale_fill_cmocean(expression(paste("NPP (", g~C~m^-2~d^-1, ")")), name = "algae", direction = 1) +
  geom_sf(data = gom.sf, linewidth = 0.25) +
  coord_sf(xlim = c(bbox_p[1], bbox_p[2]),
           ylim = c(bbox_p[3], bbox_p[4]),
           expand = FALSE,
           label_axes = "----") +
  theme_bw() +
  theme(legend.position = "top")

p.gom.sst <- ggplot() +
  geom_spatraster(data = cov_list$sst, aes(fill = `2019-10-01`)) +
  scale_fill_cmocean("SST (°C)", name = "thermal", direction = 1) +
  geom_sf(data = gom.sf, linewidth = 0.25) +
  coord_sf(xlim = c(bbox_p[1], bbox_p[2]),
           ylim = c(bbox_p[3], bbox_p[4]),
           expand = FALSE,
           label_axes = "----") +
  theme_bw() +
  theme(legend.position = "top")

# Make composite plot
p.gom.bathym + p.gom.npp + p.gom.sst +
  plot_layout(nrow = 1)

# ggsave("Tables_Figs/Figure S4.png", width = 8, height = 3, units = "in", dpi = 400)


#####################################################################
### Generate random (available) points and extract environ covars ###
#####################################################################

# Remove points that intersect land
tic()
dat2 <- dat %>%
  split(.$id) %>%
  map(~{.x %>%
      st_as_sf(coords = c('x','y'), crs = 3395) %>%
      st_mask(gom.sf) %>%
      mutate(x = st_coordinates(.)[,1],
             y = st_coordinates(.)[,2]) %>%
      st_drop_geometry()}) %>%
  bind_rows()
toc()  #took 10 sec

# How many points per ID?
table(dat2$id)  #min of 38; max of 1961


# Mask the bbox by the land layer to generate available pts in water
bbox_mask <- st_mask(bbox, gom.sf)

# Check polygon of availability
ggplot() +
  geom_sf(data = bbox, fill = "lightblue") +
  geom_sf(data = gom.sf) +
  geom_point(data = dat2, aes(x, y, color = behav)) +
  geom_sf(data = bbox_mask, color = "red", fill = "transparent") +
  scale_color_viridis_d() +
  theme_bw()


# Convert to class for use of {amt} functions
dat3 <- make_track(dat2, .x = x, .y = y, .t = date, id = id, month.year = month.year,
                   behav = behav, crs = 3395) %>%
  nest(data = -id)
# dat3 <- dat2 %>%
#   nest(data = -id)

# Create KDE_href for each ID and add as nested column
tic()
dat4 <- dat3 %>%
  mutate(ud = map(data, ~hr_akde(.x, levels = 0.99)))
toc()  #took 40 sec to run

# Add buffer to extend available area
mean(dat2$step, na.rm = TRUE) * 6  #check avg distance that could be covered in a day (based on 4 hr time step)

dat4 <- dat4 %>%
  mutate(ud_buff = map(ud, ~{.x %>%
      hr_isopleths() %>%
      slice(2) %>%
      sf::st_buffer(dist = 1e4)}  #add 10 km buffer (avg daily distance)
      ))

# Clip UDs by land mask
dat4 <- dat4 %>%
  mutate(ud_buff = map(ud_buff, ~st_mask(.x, gom.sf)))


# Subset only 'Resident' locations
dat4 <- dat4 %>%
  mutate(data_res = map(data, ~{.x %>%
      filter(behav == 'Resident')}
      ))
# dat4 <- dat3 %>%
#   mutate(data_res = map(data, ~{.x %>%
#       filter(behav == 'Resident')}
#   ))


# Generate available points and randomly assign to month.year per ID
set.seed(2023)
tic()
dat4 <- dat4 %>%
  mutate(avail_pts10 = map2(.x = ud_buff,
                          .y = data_res,
                          ~{data.frame(geometry = st_sample(.x, size = 10 * nrow(.y))) %>%
                              mutate(month.year = sample(unique(.y$month.year), size = n(), replace = T)) %>%
                              st_sf()}
                          )
)
toc()  #takes 3.5 sec to run


# set.seed(2023)
# tic()
# dat4 <- dat4 %>%
#   mutate(avail_pts10 = map2(.x = cov_list$bathym,
#                             .y = data_res,
#                             ~{sample_rast_gradient(.x, n_pts = 10 * nrow(.y)) %>%
#                                 mutate(month.year = sample(unique(.y$month.year), size = n(), replace = T))
#                             }
#   )
#   )
# toc()


# Viz example of available point spread by month.year
ggplot() +
  geom_sf(data = dat4$avail_pts10[[2]], aes(color = month.year)) +
  geom_point(data = dat4$data_res[[2]], aes(x_, y_)) +
  theme_bw() +
  facet_wrap(~ month.year)


# Join all used and available points together
used <- dat4 %>%
  dplyr::select(id, data_res) %>%
  unnest(cols = data_res) %>%
  mutate(obs = 1) %>%
  rename(x = x_, y = y_) %>%
  dplyr::select(id, month.year, x, y, obs) %>%
  data.frame()

table(used$id)  #check N per ID again; min = 34, max = 1961


avail_10 <- dat4 %>%
  dplyr::select(id, avail_pts10) %>%
  mutate(avail_pts10 = map(avail_pts10, ~{.x %>%
      mutate(x = st_coordinates(.)[,1],
             y = st_coordinates(.)[,2]) %>%
      st_drop_geometry()}
      )) %>%
  unnest(cols = avail_pts10) %>%
  mutate(obs = 0)



rsf.pts_10 <- rbind(used, avail_10)


# Remove ID 169273 due to unrealistic track w/ many anomalous locations
rsf.pts_10 <- rsf.pts_10 %>%
  filter(!id %in% c(169273))


# Extract environ covars by month.year
plan(multisession, workers = availableCores() - 2)
rsf.pts_10 <- extract.covars(data = rsf.pts_10, layers = cov_list, dyn_names = c('npp', 'sst'),
                           along = FALSE, ind = "month.year", imputed = FALSE)
#takes 1 min to run on desktop (24 cores)

plan(sequential)



# Viz example of available point covar values by month.year
x181796 <- rsf.pts_10 %>%
  filter(id == 181796)

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
  cmocean::scale_color_cmocean("Depth (m)", name = 'deep',
                               direction = -1, breaks = c(-250,-200,-150,-100,-50,0)) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ month.year)
# ggsave("../../Conference Presentations/SERSTM 2023/use_avail_depth_map.png", width = 6, height = 6,
#        units = "in", dpi = 400)

ggplot() +
  geom_point(data = rsf.pts_10 %>%
               filter(id == 181796, obs == 0), aes(x, y, color = sst)) +
  geom_sf(data = gom.sf) +
  geom_point(data = rsf.pts_10 %>%
               filter(id == 181796, obs == 1), aes(x, y), color = 'red') +
  cmocean::scale_color_cmocean("SST (°C)", name = 'thermal') +
  labs(x="",y="") +
  theme_bw() +
  coord_sf(xlim = c(min(x181796$x), max(x181796$x)),
           ylim = c(min(x181796$y), max(x181796$y))) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ month.year)
# ggsave("../../Conference Presentations/SERSTM 2023/use_avail_sst_map.png", width = 6, height = 6,
#        units = "in", dpi = 400)



###########################
### Export prepped data ###
###########################

# write_csv(rsf.pts_10, "Processed_data/GoM_Cm_RSFprep_10x.csv")
# write_csv(rsf.pts_30, "Processed_data/GoM_Cm_RSFprep_30x.csv")
# write_csv(rsf.pts_50, "Processed_data/GoM_Cm_RSFprep_50x.csv")
