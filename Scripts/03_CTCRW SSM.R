
### Fit CTCRW to turtle data via {foieGras} ###

library(tidyverse)
library(lubridate)
library(foieGras)  #v1.1
library(sf)  #v1.0.9
library(sfarrow)
library(tictoc)
library(plotly)
library(pathroutr)
library(patchwork)
library(trelliscopejs)
library(future)
library(furrr)


#### Load data ####

dat <- read_csv("Processed_data/Prefiltered_Cm_Tracks.csv")

glimpse(dat)
summary(dat)


#### Wrangle data for analysis using {foieGras} ####

# Convert all 'Quality' values to "G" for FastGPS data and 'Date' to datetime format
dat <- dat %>%
  mutate(Quality = ifelse(Type == 'FastGPS', 'G', Quality))


# Rename columns for {foieGras}
dat2<- dat %>%
  rename(id = Ptt, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)  #reorders and subsets the columns

glimpse(dat2)




## Define bouts per PTT where gap > 7 days (for splitting imputed tracks)
int_tbl <- dat2 %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x, bout = crawlUtils::cu_get_bouts(x = .$date, gap = 3, time_unit = "days"))) %>%
  purrr::map(., function (x)
  {
    date <- NULL
    x <- x %>%
      group_by(.data[["bout"]]) %>%
      summarize(start = min(date), end = max(date))
  })



#### Inspect time steps of transmissions for making predictions ####

tmp <- dat2 %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x,
                        dt = difftime(c(date[-1], NA),
                                      date,
                                      units = "hours") %>%
                          as.numeric())
  ) %>%
  bind_rows()


ggplot(tmp, aes(date, dt)) +
  geom_point() +
  theme_bw() +
  facet_trelliscope(~id, nrow = 5, ncol = 5, scales = "free_x")

# Determine primary time step by ID
tmp %>%
  group_by(id) %>%
  summarize(mean = mean(dt, na.rm = TRUE),
            median = median(dt, na.rm = TRUE)) %>%
  data.frame()
# Mean/median time step is ~1-2 hrs for most IDs




#################
#### Run SSM ####
#################

#### Account for location error at observed irregular time interval ####

# Estimate 'true' locations regularized at 2 hr time step
set.seed(2022)
tic()
fit_crw <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = 2,
                          control = ssm_control(verbose = 1))
toc()  #took 5.5 min

print(fit_crw, n = nrow(fit_crw))  #all indiv. models converged
summary(fit_crw)

# Viz the filtered outliers (black x's), raw observations (blue points), and estimated locations (orange points), along with associated uncertainty (orange shading)
plot(fit_crw, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)



################################################################
### Filter out long time gaps and re-route paths around land ###
################################################################

# Load land layers
gom.sf <- st_read_parquet('Environ_data/GoM_land.parquet')
br.sf <- st_read_parquet('Environ_data/Brazil_land.parquet')
qa.sf <- st_read_parquet('Environ_data/Qatar_land.parquet')

# Grab results and create data.frame
res_crw <- grab(fit_crw, what = "predicted")



# Separate tracks by region
gom.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'GoM',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)
br.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'Brazil',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)
qa.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'Qatar',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)



# Create visibility graphs
tic()
gom_vis_graph <- st_buffer(gom.tracks, dist = 50000) %>%
  st_union() %>%
  st_convex_hull() %>%
  st_intersection(gom.sf) %>%
  st_collection_extract('POLYGON') %>%
  st_sf() %>%
  prt_visgraph()
toc()  #took 28 sec

tic()
br_vis_graph <- st_buffer(br.tracks, dist = 50000) %>%
  st_union() %>%
  st_convex_hull() %>%
  st_intersection(br.sf) %>%
  st_collection_extract('POLYGON') %>%
  st_sf() %>%
  prt_visgraph()
toc()  #took 17 sec

tic()
qa_vis_graph <- st_buffer(qa.tracks, dist = 50000) %>%
  st_union() %>%
  st_convex_hull() %>%
  st_intersection(qa.sf) %>%
  st_collection_extract('POLYGON') %>%
  st_sf() %>%
  prt_visgraph()
toc()  #took 2 sec




# Filter by the defined bout periods
gom.tracks2 <- gom.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'GoM',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows() #%>%
  # filter(!is.na(bout))  #remove predictions that don't belong to a bout period

br.tracks2 <- br.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'Brazil',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()

qa.tracks2 <- qa.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'Qatar',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()


# Remove predictions that begin or end on land
gom.tracks3 <- prt_trim(gom.tracks2, gom.sf)
br.tracks3 <- prt_trim(br.tracks2, br.sf)
qa.tracks3 <- prt_trim(qa.tracks2, qa.sf)




# Convert tracks to list for re-routing
gom.tracks3.list <- gom.tracks3 %>%
  split(.$id)
br.tracks3.list <- br.tracks3 %>%
  split(.$id)
qa.tracks3.list <- qa.tracks3 %>%
  split(.$id)


# Re-route GoM tracks
plan(multisession)
tic()
gom.tracks.route <- gom.tracks3.list %>%
  future_map(., ~{.x %>%
      prt_trim(gom.sf) %>%
      prt_reroute(gom.sf, gom_vis_graph, blend=FALSE)},
             .options = furrr_options(packages = c("sf","dplyr"),
                                      seed = 2022)
             )
toc()  #took 1.5 min


gom.tracks.route <- gom.tracks.route %>%
  future_map2(.x = ., .y = gom.tracks3.list, ~prt_update_points(.x, .y),
              .options = furrr_options(seed = 2022)
  )
plan(sequential)


# Re-route Brazil tracks
plan(multisession)
tic()
br.tracks.route <- br.tracks3.list %>%
  future_map(., ~{.x %>%
      prt_trim(br.sf) %>%
      prt_reroute(br.sf, br_vis_graph, blend=FALSE)},
      .options = furrr_options(packages = c("sf","dplyr"),
                               seed = 2022)
  )
toc()  #took 1.5 min


br.tracks.route <- br.tracks.route %>%
  future_map2(.x = ., .y = br.tracks3.list, ~prt_update_points(.x, .y),
              .options = furrr_options(seed = 2022)
  )
plan(sequential)


# Re-route Qatar tracks
plan(multisession)
tic()
qa.tracks.route <- qa.tracks3.list %>%
  future_map(., ~{.x %>%
      prt_trim(qa.sf) %>%
      prt_reroute(qa.sf, qa_vis_graph, blend=FALSE)},
      .options = furrr_options(packages = c("sf","dplyr"),
                               seed = 2022)
  )
toc()  #took 4 sec


qa.tracks.route <- qa.tracks.route %>%
  future_map2(.x = ., .y = qa.tracks3.list, ~prt_update_points(.x, .y),
              .options = furrr_options(seed = 2022)
  )
plan(sequential)



# Change from list of sf objects back to single data.frame
gom.tracks.route <- bind_rows(gom.tracks.route) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

br.tracks.route <- bind_rows(br.tracks.route) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

qa.tracks.route <- bind_rows(qa.tracks.route) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()


#############################
### Inspect mapped tracks ###
#############################

# GoM
ggplot() +
  geom_sf(data = gom.sf) +
  geom_path(data = gom.tracks.route, aes(x, y, group = id, color = id)) +  #modeled tracks
  theme_bw()

# Brazil
ggplot() +
  geom_sf(data = br.sf) +
  geom_path(data = br.tracks.route, aes(x, y, group = id, color = id)) +  #modeled tracks
  theme_bw()

# Qatar
ggplot() +
  geom_sf(data = qa.sf) +
  geom_path(data = qa.tracks.route, aes(x, y, group = id, color = id)) +  #modeled tracks
  theme_bw()


# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf) +
    geom_path(data = gom.tracks.route, aes(x, y, group = id, color = id), alpha = 0.7) +  #modeled tracks
    theme_bw()
)

plotly::ggplotly(
  ggplot() +
    geom_sf(data = br.sf) +
    geom_path(data = br.tracks.route, aes(x, y, group = id, color = id), alpha = 0.7) +  #modeled tracks
    theme_bw()
)

plotly::ggplotly(
  ggplot() +
    geom_sf(data = qa.sf) +
    geom_path(data = qa.tracks.route, aes(x, y, group = id, color = id), alpha = 0.7) +  #modeled tracks
    theme_bw()
)




### Export fitted tracks ###

write.csv(gom.tracks.route, "Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foieGras.csv", row.names = FALSE)
write.csv(br.tracks.route, "Processed_data/Processed_Brazil_Cm_Tracks_SSM_2hr_foieGras.csv", row.names = FALSE)
write.csv(qa.tracks.route, "Processed_data/Processed_Qatar_Cm_Tracks_SSM_2hr_foieGras.csv", row.names = FALSE)
