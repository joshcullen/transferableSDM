
#############################################################
#### Consolidate Domit Brazil Green Turtle Tracking Data ####
#############################################################


library(tidyverse)
library(lubridate)
library(bayesmove)


### Load files ###

# 'Locations' files containing Argos data with or without FastlocGPS records
dat.loc<- list.files(path = "./Raw_data/Raw/Domit_Brazil",
                     pattern = "*Locations.csv",
                     full.names = T) %>%
  map_df(~read.csv(.))





### Clean and filter data ###

# Check date formats for each ID
dat.loc %>%
  group_by(Ptt) %>%
  dplyr::select(Date) %>%
  slice(1)


dat.loc2<- dat.loc %>%
  mutate(Date = as_datetime(Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")) %>%
  drop_na(Date, Latitude, Longitude) %>%  #remove rows w/ missing coords
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()




### Explore and compare datasets ###

summary(dat.loc2)
glimpse(dat.loc2)
table(dat.loc2$Quality)
table(dat.loc2$Ptt)



### View in {bayesmove} Shiny app ###

dat.loc2 %>%
  sf::st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>%
  mapview::mapview(., zcol = 'Ptt')

dat.loc2 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  # Filter(x = ., f = function(x) !all(is.na(x))) %>%  #remove all cols w/ ONLY NAs
  dplyr::select(id, date, x, y) %>%
  shiny_tracks(., 4326)

# checks out




### Inspect summary stats per ID

dat.loc2 %>%
  group_by(Ptt) %>%
  summarize(start = min(Date),
            end = max(Date),
            DaysAtLib = end - start,
            n = n())


### Export data

glimpse(dat.loc2)
write.csv(dat.loc2, "Raw_data/Brazil_Domit_Cm_2016_2020.csv", row.names = FALSE)
