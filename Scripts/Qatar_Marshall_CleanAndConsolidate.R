

#########################################################################################
#### Convert Marshall Green Turtle Data into Proper Format (DMS -> DD; .xls -> .csv) ####
#########################################################################################


library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(mapview)
library(readxl)


### Load files

#Specify columns to keep (since not all files have same columns)
cols<- c('Platform ID', 'Latitude', 'Longitude', 'Loc. quality', 'Loc. date', 'Error radius',
         'Semi-major axis', 'Semi-minor axis', 'Ellipse orientation')


dat <- list.files(path = './Raw_data/Raw/Marshall_Qatar/', pattern = "*.xls", full.names = T) %>%
  map(~read_excel(.) %>%  #read .xls files
        dplyr::select(all_of(cols))  #only select pre-specified columns
      ) %>%
  bind_rows()  #convert from list to data.frame

dat<- dat %>%
  rename(Ptt = `Platform ID`,  #modify names
         Quality = `Loc. quality`,
         Date = `Loc. date`,
         Error.radius = `Error radius`,
         Error.Semi.major.axis = `Semi-major axis`,
         Error.Semi.minor.axis = `Semi-minor axis`,
         Error.Ellipse.orientation = `Ellipse orientation`)



### Clean and filter out only green turtle IDs

id<- c(135131:135136, 135138:135139, 144320:144321)

dat2<- dat %>%
  filter(Ptt %in% id) %>%  #filter by ID
  drop_na(Latitude, Longitude) %>%  #remove rows w/ missing coords
  mutate(Date = as_datetime(Date, format = "%d.%m.%y %H:%M:%S", tz = 'UTC')) %>%
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()



### Convert from DMS to DD

dat3<- dat2 %>%
  mutate(Latitude = sp::char2dms(Latitude, chd = '°', chm = "'", chs = '"') %>%
           sp::as.numeric.DMS(),
         Longitude = sp::char2dms(Longitude, chd = '°', chm = "'", chs = '"') %>%
           sp::as.numeric.DMS())

# Change 'Error' values from class character to numeric
dat3<- dat3 %>%
  mutate(across(Error.radius:Error.Ellipse.orientation, as.numeric))


### Verify coordinates using map
dat3 %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>%
  mapview(., zcol = 'Ptt')

dat3 %>%
  rename(id = Ptt, x = Longitude, y = Latitude, date = Date) %>%
  bayesmove::shiny_tracks(., epsg = 4326)

#checks out



### Inspect summary stats per ID

dat3 %>%
  group_by(Ptt) %>%
  summarize(start = min(Date),
            end = max(Date),
            DaysAtLib = end - start,
            n = n())


### Export data

glimpse(dat3)
write.csv(dat3, "Raw_data/Qatar_Marshall_Cm_2014_2015.csv", row.names = FALSE)
