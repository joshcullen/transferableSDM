

### Fit CTCRW to turtle data via {crawl} and impute paths ###

library(tidyverse)
library(lubridate)
library(future)
library(furrr)
library(terra)
library(sf)
library(ptolemy)
library(crawlUtils)


### Load turtle track data ###

dat<- read.csv("Raw_data/Master Sat Tag Dataset.csv") %>%
  mutate(Date = as_datetime(Date))


dat.gom <- dat %>%
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


plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf,
            fill = "grey60", size = 0.2) +
    geom_sf(data = dat.gom.sf, aes(color = Ptt),
            alpha = 0.1)
)



## Use PTT 181807 and 181800 as examples

dat.tracks <- dat.gom.sf %>%
  filter(Ptt %in% c(181800, 181807)) %>%
  janitor::clean_names() %>%
  rename(datetime = date) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  group_by(ptt) %>%
  arrange(datetime) %>%
  ungroup()

# Fit CTCRW model
plan("multisession", workers = 2)

dat.tracks.crawl <- cu_add_argos_cols(dat.tracks) %>%
  bayesmove::df_to_list(., ind = "ptt")
dat.tracks.crw <- cu_crw_argos(data_list = dat.tracks.crawl, bm = FALSE)


# Make predictions at 1 hr time interval
# dat.tracks.pred <- cu_batch_predict(fit_list = dat.tracks.crw, predTime = "1 hour",
#                                     barrier = gom.sf)
nsims <- 10
dat.tracks.sims <- cu_crw_sample(fit_list = dat.tracks.crw, predTime = "4 hours", barrier = gom.sf,
                                 size = nsims)
# simulator <- crawl::crwSimulator(dat.tracks.crw[[1]], predTime="1 hour", parIS=0)



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
  slice(1:400)

ggplot() +
  # geom_sf(data = gom.sf) +
  geom_path(data = test, aes(mu.x, mu.y, group = rep, color = rep),
            size = 0.15) +
  theme_bw() #+
  # coord_sf()

ggplot() +
  # geom_sf(data = gom.sf) +
  geom_sf(data = dat.gom.sf %>% filter(Ptt == 181800), aes(color = Date),
            size = 0.15) +
  scale_color_viridis_c(option = "turbo") +
  # geom_path(data = test, aes(mu.x, mu.y, group = rep, color = rep),
  #           size = 0.5) +
  theme_bw() #+
# coord_sf()
