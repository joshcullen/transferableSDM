
### Fit CTCRW to turtle data via {foieGras} ###

library(tidyverse)
library(lubridate)
library(aniMotum)  #v1.1
library(sf)  #v1.0.9
library(sfarrow)
library(tictoc)
library(plotly)
library(trelliscopejs)

source('Scripts/helper functions.R')


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




## Define bouts per PTT where gap > 3 days (for splitting fitted tracks)
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
fit_crw <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = 4,
                          control = ssm_control(verbose = 1))
toc()  #took 5 min

print(fit_crw, n = nrow(fit_crw))  #all indiv. models converged
summary(fit_crw)

# Viz the filtered outliers (black x's), raw observations (blue points), and estimated locations (orange points), along with associated uncertainty (orange shading)
plot(fit_crw, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)



################################################################
### Filter out long time gaps and re-route paths around land ###
################################################################

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


# Generate hi-res coastline spatial layers
tic()
gom.sf <- ptolemy::extract_gshhg(data = gom.tracks, resolution = 'f', buffer = 500000)
toc()  #takes 1.5 min to run

tic()
br.sf <- ptolemy::extract_gshhg(data = br.tracks, resolution = 'f', buffer = 200000)
toc()  #takes 1 min to run

tic()
qa.sf <- ptolemy::extract_gshhg(data = qa.tracks, resolution = 'f', buffer = 100000)
toc()  #takes 1 min to run



# Filter by the defined bout periods
gom.tracks2 <- gom.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'GoM',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()

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






# Change from list of sf objects back to single data.frame
gom.tracks3 <- gom.tracks2 %>%
  drop_na(bout) %>%  #remove predictions that fall w/in 3-day time gaps
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()

br.tracks3 <- br.tracks2 %>%
  drop_na(bout) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()

qa.tracks3 <- qa.tracks2 %>%
  drop_na(bout) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()


#############################
### Inspect mapped tracks ###
#############################



# GoM
ggplotly(
  ggplot() +
    geom_sf(data = gom.sf) +
    geom_path(data = gom.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = gom.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")


# Brazil
ggplotly(
  ggplot() +
    geom_sf(data = br.sf) +
    geom_path(data = br.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = br.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")


# Qatar
ggplotly(
  ggplot() +
    geom_sf(data = qa.sf) +
    geom_path(data = qa.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = qa.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")






### Export fitted tracks ###

write.csv(gom.tracks3, "Processed_data/Processed_GoM_Cm_Tracks_SSM_12hr_aniMotum.csv", row.names = FALSE)
write.csv(br.tracks3, "Processed_data/Processed_Brazil_Cm_Tracks_SSM_12hr_aniMotum.csv", row.names = FALSE)
write.csv(qa.tracks3, "Processed_data/Processed_Qatar_Cm_Tracks_SSM_12hr_aniMotum.csv", row.names = FALSE)


## Write files for coastline spatial layers
st_write_parquet(gom.sf, 'Environ_data/GoM_land.parquet')
st_write_parquet(br.sf, 'Environ_data/Brazil_land.parquet')
st_write_parquet(qa.sf, 'Environ_data/Qatar_land.parquet')
