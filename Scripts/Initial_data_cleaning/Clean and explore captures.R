

#########################################################
#### Clean and Explore All Green Turtle Capture Data ####
#########################################################


library(tidyverse)
library(janitor)
library(lubridate)
library(bayesmove)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)


## Load data and clean column names
dat<- read.csv("Raw_data/2021-10 Master Turtle Captures.csv") %>%
  janitor::clean_names()


# Explore data
summary(dat)  #97 variables
glimpse(dat)


# Only retain relevant columns for further analysis
# 'tags column' indicates that there were 89 recaptures and 626 new recruits
dat2<- dat %>%
  filter(species == "Cm") %>%  #only keep green turtles (N = 274)
  dplyr::select(permit_number:days_at_large, flipper_tag_1:sex, satellite_tag:release_time,
                acoustic_tag:acoustic_tag_serial_no)


# Make sure that date format is consistent across all observations
dat2$date<- parse_date_time(x = dat2$date,
                  orders = c("m/d/y", "m/d/Y"),
                  tz = "UTC") #%>%
  # as.character()



# Update data frame to read into Shiny app
dat2<- dat2 %>%
  rename(id = turtle_id, x = gps_x, y = gps_y) #%>%
  # mutate(date = paste(date, time)) %>%  #need to create datetime (to get Shiny app to work)
  # mutate_at("date", ~strptime(., format = "%Y/%m/%d %H:%M") %>%
  #             as.POSIXct())



# Viz in ggplot
states<- ne_states(country = "United States of America", returnclass = "sf")
countries<- ne_countries(scale = 10, continent = c("North America", "South America"),
                         returnclass = "sf")

# All sites
ggplot() +
  geom_sf(data = countries) +
  geom_sf(data = states) +
  geom_point(data = dat2, aes(x, y, color = site), alpha = 0.75) +
  scale_color_viridis_d() +
  xlim(-100, -40) +
  ylim(-35, 35) +
  theme_bw()

# FL and Bahamas
ggplot() +
  geom_sf(data = countries) +
  geom_sf(data = states) +
  geom_point(data = dat2 %>%
               filter(y > 0), aes(x, y, color = site), alpha = 0.75) +
  scale_color_viridis_d() +
  xlim(-86, -76) +
  ylim(22, 32) +
  theme_bw()

# Brazil
ggplot() +
  geom_sf(data = countries) +
  geom_sf(data = states) +
  geom_point(data = dat2 %>%
               filter(y < 0), aes(x, y, color = site), alpha = 0.75) +
  scale_color_viridis_d() +
  xlim(-49, -48) +
  ylim(-26, -25) +
  theme_bw()


# View in Shiny app 274 total turtles (90% from Bimini and Crystal River)
dat3<- dat2 %>%
  drop_na(x) %>%
  dplyr::select(id, date, x, y, site)
shiny_tracks(dat3, 4326)

