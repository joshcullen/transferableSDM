
######################################################
#### Consolidate Sasso Green Turtle Tracking Data ####
######################################################


library(tidyverse)
library(lubridate)
library(bayesmove)


### Load files ###

## From Wildlife Computers website

# 'Locations' files containing Argos data with or without FastlocGPS records
dat.loc<- list.files(path = "./Raw_data/Raw/Sasso_GoM1",
                 pattern = "*Locations.csv",
                 full.names = T) %>%
  map_df(~read.csv(.))

# FastGPS file for which to append 'Residual' column to FastGPS observations from dat.loc
dat.fast<- list.files(path = "./Raw_data/Raw/Sasso_GoM1",
                       pattern = "*FastGPS.csv",
                       full.names = T) %>%
  map_df(~read.csv(., skip = 3)) %>%
  rename(Ptt = Name, Date = InitTime, Type = InitType)

# Remove any FastGPS locations with Residual > 30
dat.fast<- filter(dat.fast, Residual < 31 | is.na(Residual))



## From Chris Sasso email

# Remove IDs from other species (in Chazz 2016 files)
filt.id<- c(132793, 128354, 162057, 162056, 162054)

dat.chazz.loc<- list.files(path = "./Raw_data/Raw/Sasso_GoM2",
                           pattern = "*Locations.csv",
                           full.names = T) %>%
  map_df(~read.csv(.)) %>%
  filter(!(Ptt %in% filt.id))

dat.chazz.fast<- list.files(path = "./Raw_data/Raw/Sasso_GoM2",
                            pattern = "*revised.csv",
                            full.names = T) %>%
  map_df(~read.csv(.)) %>%
  rename(Ptt = Name)

# Remove any FastGPS locations with Residual > 30
dat.chazz.fast<- filter(dat.chazz.fast, Residual < 31 | is.na(Residual))





## Merge 'Satellites' and 'Residual' columns from dat.fast with dat.loc
dat.loc2<- left_join(dat.loc, dat.fast[,c("Ptt","Date","Type","Satellites","Residual")],
                     by = c('Ptt','Date','Type'))
# The number of obs increases by 5, but this is just due to duplicates being inserted from dat.fast for some reason when using left_join(); these can be filtered later during data cleaning

dat.chazz.loc2<- left_join(dat.chazz.loc, dat.chazz.fast[,c("Ptt","Latitude","Longitude",
                                                            "Satellites","Residual")],
                           by = c('Ptt','Latitude','Longitude'))
# The number of obs increases by 24, but this is just due to duplicates being inserted from dat.chazz.fast for some reason when using left_join(); these can be filtered later during data cleaning




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


# Remove fractional seconds from GPS locs
dat.chazz.loc2$Date<- gsub("\\..([^ ]*)", "", dat.chazz.loc2$Date)
dat.chazz.loc2$Date<- as_datetime(dat.chazz.loc2$Date, format = "%H:%M:%S %d-%b-%Y", tz = "UTC")

dat.chazz.loc3<- dat.chazz.loc2 %>%
  drop_na(Latitude, Longitude) %>%  #remove rows w/ missing coords
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()




### Explore and compare datasets ###

summary(dat.loc3); summary(dat.chazz.loc3)
glimpse(dat.loc3); glimpse(dat.chazz.loc3)
table(dat.loc3$Quality); table(dat.chazz.loc3$Quality)
table(dat.loc3$Ptt); table(dat.chazz.loc3$Ptt)

# Since chazz files have GPS data for duplicate IDs (as opposed to data downloaded directly from Wildlife Computers), keep this set of data and remove the rest from dat.loc3


ind<- unique(dat.loc3$Ptt)[which(unique(dat.loc3$Ptt) %in% unique(dat.chazz.loc3$Ptt))]
dat.loc4<- dat.loc3 %>%
  filter(!(Ptt %in% ind)) %>%
  dplyr::select(-Count) %>%
  rbind(., dat.chazz.loc3)

dat.loc4$Type<- ifelse(dat.loc4$Type != "GPS", dat.loc4$Type, "FastGPS")



### View in {bayesmove} Shiny app ###

dat.loc4 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  Filter(x = ., f = function(x) !all(is.na(x))) %>%  #remove all cols w/ ONLY NAs
  shiny_tracks(., 4326)

# checks out




### Inspect summary stats per ID

dat.loc4 %>%
  group_by(Ptt) %>%
  summarize(start = min(Date),
            end = max(Date),
            DaysAtLib = end - start,
            n = n())

# Check which IDs have Argos ONLY vs those with Argos and FastGPS
dat.loc4 %>%
  group_by(Ptt, Type) %>%
  summarize(n = n()) %>%
  data.frame()


### Export data

glimpse(dat.loc4)
write.csv(dat.loc4, "Raw_data/GoM_Sasso_Cm_2016_2017.csv", row.names = FALSE)
