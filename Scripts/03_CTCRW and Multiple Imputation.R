

### Fit CTCRW to turtle data via {crawl} and impute paths ###

library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(sf)
library(ptolemy)
library(crawlUtils)  #v0.1.07; needs to be used w/ crawl v2.2.1 otherwise it crashes
library(trip)
library(doRNG)  #only needed when running sourced function cu_crw_sample()
library(tictoc)
library(progressr)
library(arrow)
library(sfarrow)

# load in modified {crawlUtils} functions (that I know work)
source("Scripts/crawlUtils functions.R")


### Load turtle track data ###

dat<- read.csv("Processed_data/Prefiltered_Cm_Tracks.csv") %>%
  mutate(date = as_datetime(date))



## Apply speed and distance-angle filter on tracks before analyzing w/ CTCRW SSM
dat$keep <- trip(obj = dat %>%
                   relocate(any_of(c("Longitude", "Latitude", "date", "Ptt")), .before = 'Source'),
                 TORnames = c("datetime", "Ptt")) %>%
  sda(x = ., smax = 3 * 3600 / 1000, ang = c(15, 25), distlim = c(2.5, 5), pre = NULL)


## Viz filtered relocations
ggplot(dat %>% filter(Ptt == 169273), aes(Longitude, Latitude, color = keep)) +
  geom_point(size = 0.5) +
  theme_bw() +
  coord_equal()


## Filter data to only retain locs that passed sda filter
dat.filt <- dat %>%
  filter(keep == TRUE) %>%
  dplyr::select(-keep)


## Separate out tracks by region
dat.gom <- dat.filt %>%
  filter(Region == 'GoM')
dat.br <- dat.filt %>%
  filter(Region == 'Brazil')
dat.qa <- dat.filt %>%
  filter(Region == 'Qatar')

dat.gom.sf <- dat.gom %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326, remove = FALSE) %>%
  st_transform(3395)
dat.br.sf <- dat.br %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326, remove = FALSE) %>%
  st_transform(3395)
dat.qa.sf <- dat.qa %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326, remove = FALSE) %>%
  st_transform(3395)


## Generate hi-res coastline spatial layers
gom.sf <- ptolemy::extract_gshhg(data = dat.gom.sf, resolution = 'f', buffer = 50000)
br.sf <- ptolemy::extract_gshhg(data = dat.br.sf, resolution = 'f', buffer = 50000)
qa.sf <- ptolemy::extract_gshhg(data = dat.qa.sf, resolution = 'f', buffer = 50000)



## Viz spatial layers and tracks
ggplot() +
  geom_sf(data = gom.sf, fill = "grey60", size = 0.2) +
  geom_sf(data = dat.gom.sf, alpha = 0.1, color = 'blue')

ggplot() +
  geom_sf(data = br.sf, fill = "grey60", size = 0.2) +
  geom_sf(data = dat.br.sf, alpha = 0.1, color = 'blue')

ggplot() +
  geom_sf(data = qa.sf, fill = "grey60", size = 0.2) +
  geom_sf(data = dat.qa.sf, alpha = 0.1, color = 'blue')



dat.gom.sf %>%
  group_by(Ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING") %>%
  mapview::mapview(zcol = "Ptt", map.types = c("Esri.WorldImagery"))


plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf,
            fill = "grey60", size = 0.2) +
    geom_sf(data = dat.gom.sf %>% filter(Ptt == 169273), aes(color = Ptt),
            alpha = 0.1)
)



## Reformat data for reading by {crawl}/{crawlUtils} functions

dat.gom.tracks <- dat.gom.sf %>%
  # filter(Ptt == 181800) %>%
  # filter(Ptt %in% c(159776, 175692, 181796, 181800, 181807)) %>%
  janitor::clean_names() %>%
  # mutate(x = st_coordinates(.)[,1],
  #        y = st_coordinates(.)[,2]) %>%
  # st_drop_geometry() %>%
  st_sf() %>%  #helps make use of dplyr::rename work below
  rename(datetime = date) %>%
  group_by(ptt) %>%
  arrange(datetime) %>%
  ungroup()

dat.br.tracks <- dat.br.sf %>%
  janitor::clean_names() %>%
  st_sf() %>%  #helps make use of dplyr::rename work below
  rename(datetime = date) %>%
  group_by(ptt) %>%
  arrange(datetime) %>%
  ungroup()

dat.qa.tracks <- dat.qa.sf %>%
  janitor::clean_names() %>%
  st_sf() %>%  #helps make use of dplyr::rename work below
  rename(datetime = date) %>%
  group_by(ptt) %>%
  arrange(datetime) %>%
  ungroup()


## Define bouts per PTT where gap > 7 days (for splitting imputed tracks)
gom_int_tbl <- dat.gom.tracks %>%
  split(.$ptt) %>%
  map(., ~mutate(.x, bout = cu_get_bouts(x = .$datetime, gap = 7, time_unit = "days"))) %>%
  map(., ~cu_bout_summary(x = .x, bout = "bout"))

br_int_tbl <- dat.br.tracks %>%
  split(.$ptt) %>%
  map(., ~mutate(.x, bout = cu_get_bouts(x = .$datetime, gap = 7, time_unit = "days"))) %>%
  map(., ~cu_bout_summary(x = .x, bout = "bout"))

qa_int_tbl <- dat.qa.tracks %>%
  split(.$ptt) %>%
  map(., ~mutate(.x, bout = cu_get_bouts(x = .$datetime, gap = 7, time_unit = "days"))) %>%
  map(., ~cu_bout_summary(x = .x, bout = "bout"))



#######################
### Fit CTCRW model ###
#######################

### GoM ###

set.seed(123)
dat.gom.tracks.crawl <- cu_add_argos_cols(dat.gom.tracks) %>%
  mutate(type = case_when(type == "Argos_kf" ~ "Argos_ls",
                          TRUE ~ type)) %>%
  split(.$ptt)

tic()
# with_progress({
  dat.gom.tracks.crw <- crawlUtils::cu_crw_argos(data_list = dat.gom.tracks.crawl, bm = FALSE)
# })
toc()  #takes 2 min


## Create Visibility Graph
gom_vis_graph <- pathroutr::prt_visgraph(gom.sf)


# Make predictions at 2 hr time interval
tic()
dat.gom.tracks.pred <- cu_crw_predict(fit_list = dat.gom.tracks.crw, predTime = "2 hours",
                                      barrier = gom.sf, vis_graph = gom_vis_graph)
toc()  #took 5 min on desktop


preds.gom <- bind_rows(dat.gom.tracks.pred)

preds.gom.sf <- preds.gom %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_sf(data = gom.sf) +
  geom_sf(data = dat.gom.tracks %>% filter(ptt %in% preds.gom$ptt), aes(color = ptt), size = 0.5) +
  geom_sf(data = preds.gom, aes(color = ptt), size = 0.5) +
  theme_bw()



# Perform multiple imputation on tracks
doFuture::registerDoFuture()
plan("multisession", workers = availableCores() - 2)
set.seed(2022)

nsims <- 20
tic()
progressr::with_progress({
  dat.gom.tracks.sims <- cu_crw_sample(fit_list = dat.gom.tracks.crw, predTime = "2 hours",
                                       size = nsims, barrier = gom.sf, vis_graph = gom_vis_graph)
})
toc()  #took 1.85 hrs to run w/ 20 imputations on desktop

plan(sequential)


# Viz multiple imputations for tracks
dat.gom.tracks.sims2 <- dat.gom.tracks.sims %>%
  set_names(unique(preds.gom$ptt)) %>%
  map(., set_names, 1:nsims) %>%
  # map_depth(2, ~{.$alpha.sim %>%
  #     as.data.frame()}) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "ptt") %>%
  mutate(rep = paste(ptt, rep, sep = "_"))

# filter by the defined bout periods
dat.gom.tracks.sims2 <- dat.gom.tracks.sims2 %>%
  split(.$ptt) %>%
  future_map2(., gom_int_tbl, ~cu_join_interval_tbl(x = .x, int_tbl = .y)) %>%
  bind_rows() %>%
  filter(!is.na(bout))  #remove predictions that don't belong to a bout period


dat.gom.tracks.sims.sf <- dat.gom.tracks.sims2 %>%
  group_by(ptt, rep, bout) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf) +
    geom_sf(data = dat.gom.tracks.sims.sf, aes(color = ptt), size = 0.15, alpha = 0.25) +
    geom_sf(data = preds.gom, aes(color = ptt), size = 0.5) +
    theme_bw() +
    coord_sf()
)






### Brazil ###

set.seed(123)
dat.br.tracks.crawl <- cu_add_argos_cols(dat.br.tracks) %>%
  mutate(type = case_when(type == "Argos_kf" ~ "Argos_ls",
                          TRUE ~ type)) %>%
  split(.$ptt)

tic()
# with_progress({
dat.br.tracks.crw <- crawlUtils::cu_crw_argos(data_list = dat.br.tracks.crawl, bm = FALSE)
# })
toc()  #takes 1.5 min


## Create Visibility Graph
br_vis_graph <- pathroutr::prt_visgraph(br.sf)


# Make predictions at 2 hr time interval
tic()
dat.br.tracks.pred <- cu_crw_predict(fit_list = dat.br.tracks.crw, predTime = "2 hours",
                                      barrier = br.sf, vis_graph = br_vis_graph)
toc()  #took 5.5 min on laptop


preds.br <- bind_rows(dat.br.tracks.pred)

preds.br.sf <- preds.br %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_sf(data = br.sf) +
  geom_sf(data = dat.br.tracks %>% filter(ptt %in% preds.br$ptt), aes(color = ptt), size = 0.5) +
  geom_sf(data = preds.br, aes(color = ptt), size = 0.5) +
  theme_bw()



# Perform multiple imputation on tracks
doFuture::registerDoFuture()
plan("multisession", workers = availableCores() - 2)
set.seed(2022)

nsims <- 20
tic()
progressr::with_progress({
  dat.br.tracks.sims <- cu_crw_sample(fit_list = dat.br.tracks.crw, predTime = "2 hours",
                                       size = nsims, barrier = br.sf, vis_graph = br_vis_graph)
})
toc()  #took 71 min to run w/ 20 imputations on laptop

plan(sequential)


# Viz multiple imputations for tracks
dat.br.tracks.sims2 <- dat.br.tracks.sims %>%
  set_names(unique(preds.br$ptt)) %>%
  map(., set_names, 1:nsims) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "ptt") %>%
  mutate(rep = paste(ptt, rep, sep = "_"))

# filter by the defined bout periods
dat.br.tracks.sims2 <- dat.br.tracks.sims2 %>%
  split(.$ptt) %>%
  future_map2(., br_int_tbl, ~cu_join_interval_tbl(x = .x, int_tbl = .y)) %>%
  bind_rows() %>%
  filter(!is.na(bout))  #remove predictions that don't belong to a bout period


dat.br.tracks.sims.sf <- dat.br.tracks.sims2 %>%
  group_by(ptt, rep, bout) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

plotly::ggplotly(
  ggplot() +
    geom_sf(data = br.sf) +
    geom_sf(data = dat.br.tracks.sims.sf, aes(color = ptt), size = 0.15, alpha = 0.25) +
    geom_sf(data = preds.br, aes(color = ptt), size = 0.5) +
    theme_bw() +
    coord_sf()
)






### Qatar ###

set.seed(123)
dat.qa.tracks.crawl <- dat.qa.tracks %>%
  dplyr::select(-c(error_radius:error_ellipse_orientation)) %>%
  cu_add_argos_cols() %>%
  # mutate(type = case_when(type == "Argos_kf" ~ "Argos_ls",
  #                         TRUE ~ type)) %>%
  split(.$ptt)

tic()
# with_progress({
dat.qa.tracks.crw <- crawlUtils::cu_crw_argos(data_list = dat.qa.tracks.crawl, bm = FALSE)
# })
toc()  #takes 7.5 sec


## Create Visibility Graph
qa_vis_graph <- pathroutr::prt_visgraph(qa.sf)


# Make predictions at 2 hr time interval
tic()
dat.qa.tracks.pred <- cu_crw_predict(fit_list = dat.qa.tracks.crw, predTime = "2 hours",
                                     barrier = qa.sf, vis_graph = qa_vis_graph)
toc()  #took 7.5 sec on laptop


preds.qa <- bind_rows(dat.qa.tracks.pred)

preds.qa.sf <- preds.qa %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_sf(data = qa.sf) +
  geom_sf(data = dat.qa.tracks %>% filter(ptt %in% preds.qa$ptt), aes(color = ptt), size = 0.5) +
  geom_sf(data = preds.qa, aes(color = ptt), size = 0.5) +
  theme_bw()



# Perform multiple imputation on tracks
doFuture::registerDoFuture()
plan("multisession", workers = availableCores() - 2)
set.seed(2022)

nsims <- 20
tic()
progressr::with_progress({
  dat.qa.tracks.sims <- cu_crw_sample(fit_list = dat.qa.tracks.crw, predTime = "2 hours",
                                       size = nsims, barrier = qa.sf, vis_graph = qa_vis_graph)
})
toc()  #took 2.5 min to run w/ 20 imputations on laptop

plan(sequential)


# Viz multiple imputations for tracks
dat.qa.tracks.sims2 <- dat.qa.tracks.sims %>%
  set_names(unique(preds.qa$ptt)) %>%
  map(., set_names, 1:nsims) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "ptt") %>%
  mutate(rep = paste(ptt, rep, sep = "_"))

# filter by the defined bout periods
dat.qa.tracks.sims2 <- dat.qa.tracks.sims2 %>%
  split(.$ptt) %>%
  future_map2(., qa_int_tbl, ~cu_join_interval_tbl(x = .x, int_tbl = .y)) %>%
  bind_rows() %>%
  filter(!is.na(bout))  #remove predictions that don't belong to a bout period


dat.qa.tracks.sims.sf <- dat.qa.tracks.sims2 %>%
  group_by(ptt, rep, bout) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

plotly::ggplotly(
  ggplot() +
    geom_sf(data = qa.sf) +
    geom_sf(data = dat.qa.tracks.sims.sf, aes(color = ptt), size = 0.15, alpha = 0.25) +
    geom_sf(data = preds.qa, aes(color = ptt), size = 0.5) +
    theme_bw() +
    coord_sf()
)







### Export results ###

# Remove geometry col from sf objects
preds.gom.df <- preds.gom %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
dat.gom.tracks.sims.df <- dat.gom.tracks.sims2 %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

preds.br.df <- preds.br %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
dat.br.tracks.sims.df <- dat.br.tracks.sims2 %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

preds.qa.df <- preds.qa %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()
dat.qa.tracks.sims.df <- dat.qa.tracks.sims2 %>%
  mutate(mu.x = st_coordinates(.)[,1],
         mu.y = st_coordinates(.)[,2]) %>%
  st_drop_geometry()


## Write files for best-fit and imputed tracks
# write.csv(preds.gom.df, "Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write.csv(dat.gom.tracks.sims.df, "Processed_data/Imputed_GoM_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write_parquet(preds.gom.df, "Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr.parquet")
# write_parquet(dat.gom.tracks.sims.df, "Processed_data/Imputed_GoM_Cm_Tracks_SSM_2hr.parquet")

# write.csv(preds.br.df, "Processed_data/Processed_Brazil_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write.csv(dat.br.tracks.sims.df, "Processed_data/Imputed_Brazil_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write_parquet(preds.br.df, "Processed_data/Processed_Brazil_Cm_Tracks_SSM_2hr.parquet")
# write_parquet(dat.br.tracks.sims.df, "Processed_data/Imputed_Brazil_Cm_Tracks_SSM_2hr.parquet")

# write.csv(preds.qa.df, "Processed_data/Processed_Qatar_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write.csv(dat.qa.tracks.sims.df, "Processed_data/Imputed_Qatar_Cm_Tracks_SSM_2hr.csv", row.names = FALSE)
# write_parquet(preds.qa.df, "Processed_data/Processed_Qatar_Cm_Tracks_SSM_2hr.parquet")
# write_parquet(dat.qa.tracks.sims.df, "Processed_data/Imputed_Qatar_Cm_Tracks_SSM_2hr.parquet")


## Write files for coastline spatial layers
# st_write_parquet(gom.sf, 'Environ_data/GoM_land.parquet')
# st_write_parquet(br.sf, 'Environ_data/Brazil_land.parquet')
# st_write_parquet(qa.sf, 'Environ_data/Qatar_land.parquet')

