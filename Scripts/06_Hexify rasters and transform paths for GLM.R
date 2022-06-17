

### Hexify environmental rasters and modify data for analysis by GLM ###

library(tidyverse)
library(lubridate)
library(sf)
library(terra)
library(tictoc)
library(furrr)
library(future)

source('Scripts/hexify2.R')
source('Scripts/make_glm_data2.R')



## Load data
load('Data_products/discretized_tracks.RData')


## Load in environ rasters and quadrature times
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'KdPAR', 'NPP', 'SST')

# Change names for NPP and SST to match KdPAR (YYYY-MM-01)
for (var in c('NPP', 'SST')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}

# Create list of quadrature times (names of each raster layer) per covariate
quad_list <- cov_list %>%
  map(names)
quad_list[[1]] <- quad_list$SST[1]  #need to provide a date for static bathymetry layer
names(cov_list$bathym) <- quad_list[[1]]  #also provide as name for raster layer


## Check CRS of spatial layers
st_crs(hex_grid$`128352`$poly)  #EPSG:3395
map(cov_list, crs)  #EPSG:4326
#as long as all rasters have the same CRS, hex_grid will be re-projected during extraction step




## Hexify rasters

# Create empty list to store results
hex_cov_df <- vector("list", length(hex_grid)) %>%
  set_names(names(hex_grid))

bounds_list <- disc_path %>%
  map(., ~{.x %>%
      pull(date) %>%
      range()
    })

tic()
for (i in 1:length(hex_grid)) {
  print(paste("Hexifying rasters for", names(hex_grid)[i]))
  hex_cov_df[[i]] <- hexify_df(bounds = bounds_list[[i]], raster_list = cov_list,
                               hex_list = hex_grid[[i]], quad_list = quad_list)
}
toc()  #takes 2 min on laptop




## As example, viz hexified rasters for PTT 128352 for 1st month of tracking for each covar

poly <- hex_grid$`181805`$poly %>%
  st_sf(sf_column_name = 'geometry') %>%
  mutate(hex = 1:n())
poly <- poly %>%
  left_join(.,
            hex_cov_df$`181805` %>%
                filter(Time == "2020-07-01") %>%
                dplyr::select(hex, bathym:SST),
            by = "hex")

#bathymetry
ggplot() +
  geom_sf(data = poly, aes(fill = bathym), size = 0.25) +
  cmocean::scale_fill_cmocean(name = "deep") +
  theme_bw()


#KdPAR
ggplot() +
  geom_sf(data = poly, aes(fill = KdPAR), size = 0.25) +
  cmocean::scale_fill_cmocean(name = "turbid") +
  theme_bw()


#NPP
ggplot() +
  geom_sf(data = poly, aes(fill = NPP), size = 0.25) +
  cmocean::scale_fill_cmocean(name = "algae") +
  theme_bw()


#SST
ggplot() +
  geom_sf(data = poly, aes(fill = SST), size = 0.25) +
  cmocean::scale_fill_cmocean(name = "thermal") +
  theme_bw()





### Transform data for analysis by Poisson GLM ###

# create empty list to store results
glm_data <- vector("list", length(disc_path)) %>%
  set_names(names(disc_path))



plan(multisession, workers = availableCores() - 2)
tic()
for (i in 1:length(disc_path)) {

  message("Starting ", names(disc_path)[i],"...")

  # define aesthetics of progress bar
  progressr::handlers(progressr::handler_progress(incomplete = ".", complete = "*",
                                                  current = "o", clear = FALSE))

  disc_path_list <- disc_path[[i]] %>%
    split(.$rep)

  progressr::with_progress({
    #set up progress bar
    p<- progressr::progressor(steps = length(disc_path_list))

    glm_data_list <- future_map(disc_path_list,
                           ~make_glm_data(disc_path = ., neighbor_df = hex_grid[[i]]$neighbor_df,
                                          hex_cov_df = hex_cov_df[[i]], p = p),
                           .options = furrr_options(seed = 2022)
    )
  })

  glm_data[[i]] <- bind_rows(glm_data_list)
}
toc()  # takes 5.5 min to run on laptop
plan(sequential)  #return to single core



### Export hexified rasters ###

# write.csv(hex_cov_df, "Data_products/Hexified environ covariates.csv", row.names = FALSE)
