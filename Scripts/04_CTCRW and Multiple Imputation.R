

### Fit CTCRW to turtle data via {crawl} and impute paths ###

library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(terra)
library(sf)
library(ptolemy)
library(crawlUtils)
library(trip)
library(doRNG)  #only needed when running sourced function cu_crw_sample()
library(tictoc)
library(progressr)

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


## Filter data to only retain locs from GoM and that passed sda filter

dat.filt <- dat %>%
  filter(keep == TRUE) %>%
  dplyr::select(-keep)

dat.gom <- dat.filt %>%
  filter(Region == 'GoM' & Longitude < -74 & Longitude > -100)

dat.gom.sf <- dat.gom %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326, remove = FALSE) %>%
  st_transform(3395)


gom.sf <- ptolemy::extract_gshhg(data = dat.gom.sf, resolution = 'f')


ggplot() +
  geom_sf(data = gom.sf, fill = "grey60", size = 0.2) +
  geom_sf(data = dat.gom.sf, alpha = 0.1, color = 'blue')

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



## Use PTT 159776, 181807 and 181800 as examples

dat.tracks <- dat.gom.sf %>%
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


## Define bouts per PTT where gap > 7 days (for splitting imputed tracks)
int_tbl <- dat.tracks %>%
  split(.$ptt) %>%
  map(., ~mutate(.x, bout = cu_get_bouts(x = .$datetime, gap = 7, time_unit = "days"))) %>%
  map(., ~cu_bout_summary(x = .x, bout = "bout"))




# Fit CTCRW model
doFuture::registerDoFuture()
plan("multisession")

set.seed(123)
dat.tracks.crawl <- cu_add_argos_cols(dat.tracks) %>%
  mutate(type = case_when(type == "Argos_kf" ~ "Argos_ls",
                          TRUE ~ type)) %>%
  split(.$ptt)

tic()
dat.tracks.crw <- crawlUtils::cu_crw_argos(data_list = dat.tracks.crawl, bm = FALSE)
toc()  #takes 2 min


## Create Visibility Graph
vis_graph <- pathroutr::prt_visgraph(gom.sf)


# Make predictions at 2 hr time interval
tic()
dat.tracks.pred <- cu_crw_predict(fit_list = dat.tracks.crw, predTime = "2 hours",
                                      barrier = gom.sf, vis_graph = vis_graph)
toc()  #took 5 min to run


preds <- bind_rows(dat.tracks.pred) %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  # geom_sf(data = gom.sf) +
  geom_sf(data = dat.tracks %>% filter(ptt %in% preds$ptt), aes(color = ptt), size = 0.5) +
  geom_sf(data = preds, aes(color = ptt), size = 0.5) +
  theme_bw()



# Perform multiple imputation on tracks
doFuture::registerDoFuture()
plan("multisession")
set.seed(123)

nsims <- 5
tic()
progressr::with_progress({
  dat.tracks.sims <- cu_crw_sample(fit_list = dat.tracks.crw[1], predTime = "2 hours", size = nsims,
                                   barrier = gom.sf, vis_graph = vis_graph)
})
toc()



# Viz multiple imputations for tracks

foo <- dat.tracks.sims %>%
  set_names(preds$ptt[1]) %>%
  map(., set_names, 1:nsims) %>%
  # map_depth(2, ~{.$alpha.sim %>%
  #     as.data.frame()}) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "ptt") %>%
  mutate(rep = paste(ptt, rep, sep = "_"))

# filter by the defined bout periods
foo.bout <- foo %>%
  split(.$ptt) %>%
  map2(., int_tbl[1], ~cu_join_interval_tbl(x = .x, int_tbl = .y)) %>%
  bind_rows()


test <- foo.bout %>%
  filter(!is.na(bout)) %>%  #remove predictions that don't belong to a bout period
  group_by(ptt, rep, bout) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf) +
    geom_sf(data = test, aes(color = ptt), size = 0.15, alpha = 0.25) +
    geom_sf(data = preds, aes(color = ptt), size = 0.5) +
    theme_bw() +
    coord_sf()
)


