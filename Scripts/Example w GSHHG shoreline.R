
## Load in GSHHG shoreline data to pre-specified bbox and for a selected reoslution ##


library(sf)
library(ptolemy)
library(tidyverse)
library(lubridate)



## Load turtle data
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
  geom_sf(data = dat.gom.sf,
          alpha = 0.1, color = 'blue')
)
