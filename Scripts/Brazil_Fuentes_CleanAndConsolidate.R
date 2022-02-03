
###############################################################
#### Consolidate Fuentes Brazil Green Turtle Tracking Data ####
###############################################################


library(tidyverse)
library(lubridate)
library(bayesmove)


### Load files ###

# 'Locations' files containing Argos data with or without FastlocGPS records
dat.loc<- list.files(path = "./Raw_data/Raw/Fuentes_Brazil",
                     pattern = "*Locations.csv",
                     full.names = T) %>%
  map(~read.csv(.)) %>%
  map(., ~mutate(., across(c('DeployID', 'Ptt'), as.character))) %>%  #need DeployID and Ptt to be character for all datasets since some were given names (i.e., character strings) and I have two 160105 IDs to deal with
  bind_rows()

# FastGPS file for which to append 'Residual' column to FastGPS observations from dat.loc
dat.fast<- list.files(path = "./Raw_data/Raw/Fuentes_Brazil",
                      pattern = "*FastGPS.csv",
                      full.names = T) %>%
  map_df(~read.csv(., skip = 3)) %>%
  rename(DeployID = Name, Date = InitTime)

# Remove any FastGPS locations with Residual > 30
dat.fast2<- filter(dat.fast, Residual < 31 | is.na(Residual))





## Merge 'Satellites' and 'Residual' columns from dat.fast with dat.loc
dat.loc2<- left_join(dat.loc,
                     dat.fast2[,c("DeployID","Date","Satellites","Residual")],
                     by = c('DeployID','Date'))
# The number of obs increases by 748, but this is just due to duplicates being inserted from dat.fast2 for some reason when using left_join(); these can be filtered later during data cleaning




### Clean and filter data ###

# Check date formats for each ID
dat.loc2 %>%
  group_by(Ptt) %>%
  dplyr::select(Date) %>%
  slice(1)

# Remove fractional seconds from GPS locs
dat.loc2$Date<- gsub("\\..([^ ]*)", "", dat.loc2$Date)

# Need to also account for different datetime format of 160105_1 and 160105_2
dat.loc2$Date<- parse_date_time(dat.loc2$Date,
                                orders = c("%H:%M:%S %d-%b-%Y", "%m/%d/%y %H:%M"), tz = "UTC")

dat.loc3<- dat.loc2 %>%
  drop_na(Date, Latitude, Longitude) %>%  #remove rows w/ missing coords
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()




### Explore and compare datasets ###

summary(dat.loc3)
glimpse(dat.loc3)
table(dat.loc3$Quality)
table(dat.loc3$Ptt)

# Change Type = 'User' to 'FastGPS'
dat.loc3$Type<- ifelse(dat.loc3$Type != 'User', dat.loc3$Type, 'FastGPS')


### View in {bayesmove} Shiny app ###

dat.loc3 %>%
  sf::st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>%
  mapview::mapview(., zcol = 'Ptt')

dat.loc3 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  Filter(x = ., f = function(x) !all(is.na(x))) %>%  #remove all cols w/ ONLY NAs
  # dplyr::select(id, date, x, y) %>%
  dplyr::select(-c(Comment)) %>%
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
write.csv(dat.loc3, "Raw_data/Brazil_Fuentes_Cm_2016_2021.csv", row.names = FALSE)
