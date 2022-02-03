
############################################################
#### Consolidate Fuentes GoM Green Turtle Tracking Data ####
############################################################


library(tidyverse)
library(lubridate)
library(bayesmove)


### Load files ###

# 'Locations' files containing Argos data with or without FastlocGPS records
dat.loc<- list.files(path = "./Raw_data/Raw/Fuentes_GoM",
                     pattern = "*Locations.csv",
                     full.names = T) %>%
  map_df(~read.csv(.))

# FastGPS file for which to append 'Residual' column to FastGPS observations from dat.loc
dat.fast<- list.files(path = "./Raw_data/Raw/Fuentes_GoM",
                      pattern = "*FastGPS.csv",
                      full.names = T) %>%
  map_df(~read.csv(., skip = 3)) %>%
  rename(Ptt = Name)

# Remove any FastGPS locations with Residual > 30
dat.fast2<- filter(dat.fast, Residual < 31 | is.na(Residual))





## Merge 'Satellites' and 'Residual' columns from dat.fast with dat.loc
dat.loc2<- left_join(dat.loc,
                     dat.fast2[,c("Ptt","Latitude","Longitude","Satellites","Residual")],
                     by = c('Ptt','Latitude','Longitude'))
# The number of obs increases by 103, but this is just due to duplicates being inserted from dat.fast2 for some reason when using left_join(); these can be filtered later during data cleaning




### Clean and filter data ###

# Remove fractional seconds from GPS locs
dat.loc2$Date<- gsub("\\..([^ ]*)", "", dat.loc2$Date)
dat.loc2$Date<- as_datetime(dat.loc2$Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")

dat.loc3<- dat.loc2 %>%
  drop_na(Latitude, Longitude) %>%  #remove rows w/ missing coords
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()




### Explore and compare datasets ###

summary(dat.loc3)
glimpse(dat.loc3)
table(dat.loc3$Quality)
table(dat.loc3$Ptt)



### View in {bayesmove} Shiny app ###

dat.loc3 %>%
  sf::st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>%
  mapview::mapview(., zcol = 'Ptt')

dat.loc3 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  Filter(x = ., f = function(x) !all(is.na(x))) %>%  #remove all cols w/ ONLY NAs
  shiny_tracks(., 4326)

# checks out




### Inspect summary stats per ID

dat.loc3 %>%
  group_by(Ptt) %>%
  summarize(start = min(Date),
            end = max(Date),
            DaysAtLib = end - start,
            n = n())

# Check which IDs have Argos ONLY vs those with Argos and FastGPS
dat.loc3 %>%
  group_by(Ptt, Type) %>%
  summarize(n = n()) %>%
  data.frame()


### Export data

glimpse(dat.loc3)
write.csv(dat.loc3, "Raw_data/GoM_Fuentes_Cm_2016_2017.csv", row.names = FALSE)
