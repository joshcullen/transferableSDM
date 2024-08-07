
#### Merge and Filter All Green Turtle Tracking Data ####


library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(bayesmove)
library(wesanderson)


########################################
### Load data and check column names ###
########################################

# Read in separate data files
dat<- list.files(path = "./Raw_data",
             pattern = "(Cm)",
             full.names = T) %>%
  map(~read_csv(.))

map(dat, names)
cols<- c('Ptt', 'Date', 'Longitude', 'Latitude', 'Type', 'Quality', 'Error.radius',
         'Error.Semi.major.axis', 'Error.Semi.minor.axis', 'Error.Ellipse.orientation')

# Add 'Type' column to Meg Lamont's and Chris Marshall's data (elements 4 and 6)
dat[c(4,6)]<- dat[c(4,6)] %>%
  map(., ~mutate(., Type = "Argos"))


# Convert all Ptt cols to character and bind all list elements together
dat2<- dat %>%
  map(., ~mutate(., across(Ptt, as.character))) %>%
  map_df(., dplyr::select, all_of(cols))  #select only pertinent columns


# Load in metadata for tags
turt.meta<- read_csv("Raw_data/Turtle tag metadata.csv")


# Join metadata w/ turtle obs
dat3<- left_join(turt.meta, dat2, by = 'Ptt') %>%
  relocate(Ptt, .before = Date)

# Load in metadata for deployment dates and locations per ID
deploy <- read_csv("Raw_data/Turtle deploy metadata.csv") %>%
  mutate(Deploy_Date = parse_date_time(Deploy_Date, orders = "mdyHM", truncated = 2))


# Load land spatial layers
states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))
north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')
brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')
mid.east<- ne_countries(scale = 10, country = c("Qatar", "Iran", "Bahrain", "Saudi Arabia",
                                                "United Arab Emirates"), returnclass = 'sf')
world<- ne_countries(scale = 110, returnclass = 'sf')


####################
### Explore data ###
####################

summary(dat3)
glimpse(dat3)

table(dat3$Type, useNA = 'ifany')
table(dat3$Ptt, useNA = 'ifany')
table(dat3$Quality, useNA = 'ifany')

# Look at no. of obs by region and age class (before filtering by behavioral state)
dat3 %>%
  group_by(Region, Age) %>%
  tally()

## Brazil, Adult: 59386
## Brazil, Juv:   15120
## GoM, Adult:    20283
## GoM, Adult:    27758
## Qatar, Juv:     1582



#############################
### Explore mapped tracks ###
#############################

## Overloads app with n >50k
## View by Region

# GoM
dat3 %>%
  filter(Region == 'GoM') %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  shiny_tracks(., epsg = 4326)

# Brazil
dat3 %>%
  filter(Region == 'Brazil') %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  shiny_tracks(., epsg = 4326)

# Qatar
dat3 %>%
  filter(Region == 'Qatar') %>%
  rename(id = Ptt, date = Date, x = Longitude, y = Latitude) %>%
  shiny_tracks(., epsg = 4326)


# Global View
dat3.sf<- dat3 %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326)
pal<- wes_palette('Zissou1', n = 5)[c(1,5)]

ggplot() +
  geom_sf(data = world) +
  geom_sf(data = dat3.sf, aes(color = Age), alpha = 0.3) +
  scale_color_manual(values = pal) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(crs = st_crs("+proj=moll"))






###########################
### Export cleaned data ###
###########################

write_csv(dat3, "Raw_data/Master Sat Tag Dataset.csv", row.names = FALSE)
