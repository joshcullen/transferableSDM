
#######################################################
#### Consolidate Lamont Green Turtle Tracking Data ####
#######################################################


library(tidyverse)
library(janitor)
library(lubridate)
library(sf)
library(rnaturalearth)
library(bayesmove)


### Load files

dat<- list.files(path = "./Raw_data/Raw/Lamont_GoM",
                 pattern = "*.csv",
                 full.names = T) %>%
  map(~read.csv(.))


# Check names of cols; 1st data file is different from last 4
map(dat, names)

# Change names of list elements 1, 2, 3, 7, and 9 to match others and then merge pertinent data together
dat[[1]]<- dat[[1]] %>%
  rename(Ptt = tag_id, Date = utc, Quality = lc, , Longitude = lon1, Latitude = lat1,
         Error.radius = error_radius, Error.Semi.major.axis = semi_major_axis,
         Error.Semi.minor.axis = semi_minor_axis,
         Error.Ellipse.orientation = ellipse_orientation)
dat[c(2:3,7,9)]<- dat[c(2:3,7,9)] %>%
  map({. %>%
      rename(Ptt = `Platform.ID.No.`, Quality = `Loc..quality`, Date = `Loc..date`,
             Error.Semi.major.axis = Semi.major.axis,
             Error.Semi.minor.axis = Semi.minor.axis,
             Error.Ellipse.orientation = Ellipse.orientation)})

# Check date formats
dat %>%
  map(., {. %>%
      select(Date) %>%
      slice(1)})
View(dat[[9]])  #need to double-check element 9

#elements 1, 2, 6, 8 are reported as "%m/%d/%Y %H:%M"
#elements 3, 7, 9 are reported as "%m/%d/%Y %H:%M:%S"
#elements 4, 5, 10-18 are reported as "%H:%M:%S %d-%b-%Y"


# Change 'Date' column to datetime format
dat[c(1:2,6,8)]<- dat[c(1:2,6,8)] %>%
  map({. %>%
      mutate(Date = as_datetime(Date, format = "%m/%d/%Y %H:%M", tz = 'UTC'))
  })

dat[c(3,7,9)]<- dat[c(3,7,9)] %>%
  map({. %>%
      mutate(Date = as_datetime(Date, format = "%m/%d/%Y %H:%M:%S", tz = 'UTC'))
  })

dat[c(4:5,10:18)]<- dat[c(4:5,10:18)] %>%
  map({. %>%
      mutate(Date = as_datetime(Date, format = "%H:%M:%S %d-%b-%Y", tz = 'UTC'))
  })



#Specify columns to keep (since not all files have same columns)
cols<- c('Ptt', 'Date', 'Latitude', 'Longitude', 'Quality', 'Error.radius',
         'Error.Semi.major.axis', 'Error.Semi.minor.axis', 'Error.Ellipse.orientation')

dat2<- dat %>%
  map_df(., ~dplyr::select(., all_of(cols)))  #only select pre-specified columns)




### Clean and filter data

dat3<- dat2 %>%
  drop_na(Latitude, Longitude, Date) %>%  #remove rows w/ missing coords or date
  group_by(Ptt) %>%
  arrange(Date, .by_group = TRUE) %>%  #reorder by ID by date
  distinct(Date, .keep_all = TRUE) %>%  #remove duplicate observations
  ungroup()


### Explore data
summary(dat3)
glimpse(dat3)
table(dat3$Quality)
table(dat3$Ptt)



### View in {bayesmove} Shiny app

dat3 %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  shiny_tracks(., 4326)

# checks out



### Inspect summary stats per ID

dat3 %>%
  group_by(Ptt) %>%
  summarize(start = min(Date),
            end = max(Date),
            DaysAtLib = end - start,
            n = n())



### Export data

glimpse(dat3)
write.csv(dat3, "Raw_data/GoM_Lamont_Cm_2011_2020.csv", row.names = FALSE)
