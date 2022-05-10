

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


### Load turtle track data ###

dat<- read.csv("Processed_data/Prefiltered_Cm_Tracks.csv") %>%
  mutate(date = as_datetime(date))



## Apply speed and distance-angle filter on tracks before analyzing w/ CTCRW SSM

dat$keep <- trip(obj = dat %>%
                   relocate(any_of(c("Longitude", "Latitude", "date", "Ptt")), .before = 'Source'),
                 TORnames = c("datetime", "Ptt")) %>%
  sda(x = ., smax = 3 * 3600 / 1000, ang = c(15, 25), distlim = c(2.5, 5), pre = NULL)


## Viz filtered relocations

ggplot(dat, aes(Longitude, Latitude, color = keep)) +
  geom_point(size = 0.5) +
  theme_bw() +
  coord_equal()



dat.filt <- dat %>%
  filter(keep == TRUE) %>%
  dplyr::select(-keep)



dat.gom <- dat.filt %>%
  filter(Region == 'GoM' & Longitude < -80 & Longitude > -100)

dat.gom.sf <- dat.gom %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326) %>%
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


# plotly::ggplotly(
#   ggplot() +
#     geom_sf(data = gom.sf,
#             fill = "grey60", size = 0.2) +
#     geom_sf(data = dat.gom.sf, aes(color = Ptt),
#             alpha = 0.1)
# )



## Use PTT 159776, 181807 and 181800 as examples

dat.tracks <- dat.gom.sf %>%
  filter(Ptt %in% c(159776, 175692, 181796, 181800, 181807)) %>%
  janitor::clean_names() %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  rename(datetime = date) %>%
  group_by(ptt) %>%
  arrange(datetime) %>%
  ungroup()





# Fit CTCRW model
doFuture::registerDoFuture()
plan("multisession")

set.seed(123)
dat.tracks.crawl <- cu_add_argos_cols(dat.tracks) %>%
  mutate(type = case_when(type == "Argos_kf" ~ "Argos_ls",
                          TRUE ~ type)) %>%
  bayesmove::df_to_list(., ind = "ptt")

progressr::with_progress({
  dat.tracks.crw <- cu_crw_argos(data_list = dat.tracks.crawl, bm = FALSE)
})


## Create Visibility Graph
vis_graph <- pathroutr::prt_visgraph(gom.sf)


# Make predictions at 2 hr time interval
progressr::with_progress({
  dat.tracks.pred <- cu_batch_predict(fit_list = dat.tracks.crw, predTime = "2 hours",
                                      barrier = gom.sf, vis_graph = vis_graph, crs = 3395)
})

preds <- bind_rows(dat.tracks.pred) %>%
  group_by(ptt) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")

ggplot() +
  geom_sf(data = gom.sf) +
  geom_point(data = dat.tracks, aes(x, y, color = ptt), size = 0.5) +
  scale_color_viridis_d(option = "turbo") +
  geom_sf(data = preds, aes(color = ptt), size = 0.5) +
  theme_bw()



# Perform multiple imputation on tracks
nsims <- 20
progressr::with_progress({
  dat.tracks.sims <- cu_crw_sample(fit_list = dat.tracks.crw, predTime = "2 hours", size = nsims,
                                   barrier = gom.sf, vis_graph = vis_graph, crs = 3395)
})

## create own functions to generate multiple imputations to see if this fixes anything
# proc_imp <- function(crw_fit, iter, predTime) {
#
#   simObj <- crw_fit %>%
#     crawl::crwSimulator(predTime = predTime, parIS = 0)
#
#   sim_tracks = list()
#   for (i in 1:iter) {
#     sim_tracks[[i]] <- crawl::crwPostIS(simObj, fullPost = FALSE)
#   }
#
#   sim_tracks.df <- sim_tracks %>%
#     set_names(1:iter) %>%
#     map(~{.$alpha.sim %>%
#         as.data.frame()}) %>%
#     bind_rows(.id = "rep")
#
#   return(sim_tracks.df)
# }
#
# sims <- proc_imp(crw_fit = dat.tracks.crw[[2]], iter = 20, predTime = "1 hour")


# Viz multiple imputations for tracks

foo <- dat.tracks.sims %>%
  set_names(unique(dat.tracks$ptt)) %>%
  map(., set_names, 1:nsims) %>%
  # map_depth(2, ~{.$alpha.sim %>%
  #     as.data.frame()}) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "ptt") %>%
  mutate(rep = paste(ptt, rep, sep = "_"))


test <- foo %>%
  # filter(id == 181800) %>%
  group_by(ptt, rep) %>%
  # slice(1:1000)
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

# ggplot() +
#   geom_sf(data = gom.sf) +
#   scale_color_viridis_d(option = "turbo") +
#   geom_path(data = foo, aes(mu.x, mu.y, group = rep, color = ptt),
#             size = 0.5, alpha = 0.5) +
#   geom_path(data = preds, aes(mu.x, mu.y, group = ptt, color = ptt),
#             size = 0.5) +
#   geom_point(data = dat.tracks, aes(x, y, color = ptt), size = 0.5) +
#   theme_bw() +
#   coord_sf()



## Check what's going on with imputed tracks for PTT 175692

ggplot(dat.tracks %>%
         filter(ptt == 175692), aes(datetime, y)) +
  geom_point() +
  theme_bw()
## large gap towards beginning of track; try using tools from {crawlUtils} to investigate further and potentially split track before analysis

ptt.175692 <- dat.tracks %>%
  filter(ptt == 175692) %>%
  st_as_sf(., coords = c('x','y'), crs = 3395)
ptt.175692$bout <- cu_get_bouts(x = ptt.175692$datetime, gap = 7, time_unit = "days")



cu_bout_summary <- function(x,bout){
  datetime <- NULL
  x <- x %>% group_by(.data[[bout]]) %>% st_drop_geometry() %>%
    summarize(start=min(datetime), end=max(datetime))

  return(x)
}
cu_bout_summary(x = ptt.175692, bout = "bout")


test <- dat.tracks %>%
  filter(ptt == 175692)
test2 <- cu_join_interval_tbl(x = test, int_tbl = x)
