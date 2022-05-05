

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
  geom_sf(data = gom.sf,
          fill = "grey60", size = 0.2) +
  geom_sf(data = dat.gom.sf,
          alpha = 0.1, color = 'blue')


# plotly::ggplotly(
#   ggplot() +
#     geom_sf(data = gom.sf,
#             fill = "grey60", size = 0.2) +
#     geom_sf(data = dat.gom.sf, aes(color = Ptt),
#             alpha = 0.1)
# )



## Use PTT 181807 and 181800 as examples

dat.tracks <- dat.gom.sf %>%
  filter(Ptt %in% c(159776, 181800, 181807)) %>%
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
plan("multisession", workers = 3)

set.seed(123)
dat.tracks.crawl <- cu_add_argos_cols(dat.tracks) %>%
  bayesmove::df_to_list(., ind = "ptt")
dat.tracks.crw <- cu_crw_argos(data_list = dat.tracks.crawl, bm = FALSE)


# Make predictions at 1 hr time interval
dat.tracks.pred <- cu_batch_predict(fit_list = dat.tracks.crw, predTime = "1 hour",
                                    barrier = gom.sf)

ggplot() +
  geom_sf(data = gom.sf) +
  geom_sf(data = dat.gom.sf %>% filter(Ptt == 159776), aes(color = date),
          size = 0.15) +
  scale_color_viridis_c(option = "turbo") +
  geom_path(data = dat.tracks.pred[[1]], aes(mu.x, mu.y),
            size = 0.5) +
  theme_bw() +
  coord_sf(xlim = range(dat.tracks.pred[[1]]$mu.x),
           ylim = range(dat.tracks.pred[[1]]$mu.y))



nsims <- 20
dat.tracks.sims <- cu_crw_sample(fit_list = dat.tracks.crw, predTime = "1 hour", barrier = gom.sf,
                                 size = nsims)

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
  map_depth(2, ~{.$alpha.sim %>%
      as.data.frame()}) %>%
  map(~bind_rows(.x, .id = "rep")) %>%
  bind_rows(.id = "id") %>%
  mutate(rep = paste(id, rep, sep = "_"))


test <- foo %>%
  filter(id == 181800) %>%
  group_by(rep) %>%
  slice(1:700)

# ggplot() +
#   # geom_sf(data = gom.sf) +
#   geom_path(data = test, aes(mu.x, mu.y, group = rep, color = rep),
#             size = 0.15) +
#   theme_bw() #+
#   # coord_sf()

ggplot() +
  # geom_sf(data = gom.sf) +
  geom_sf(data = dat.gom.sf %>% filter(Ptt == 181800), aes(color = date),
            size = 0.15) +
  scale_color_viridis_c(option = "turbo") +
  geom_path(data = test, aes(mu.x, mu.y, group = rep),
            size = 0.5) +
  theme_bw() #+
# coord_sf()
