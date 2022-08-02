
###############################################################
#### Create Master File for All Green Turtle Tracking Data ####
###############################################################


library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(bayesmove)
library(wesanderson)


### Load data and check column names ###
dat<- list.files(path = "./Raw_data",
             pattern = "(Cm)",
             full.names = T) %>%
  map(~read.csv(.))

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
turt.meta<- read.csv("Raw_data/Turtle tag metadata.csv")


# Join metadata w/ turtle obs
dat3<- left_join(turt.meta, dat2, by = 'Ptt') %>%
  relocate(Ptt, .before = Date)





### Explore data ###

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

dat3$Date<- as_datetime(dat3$Date)
# dat3$Quality<- ifelse(dat3$Type == "FastGPS", NA, dat3$Quality)

glimpse(dat3)
table(dat3$Quality, useNA = 'ifany')
dat3[which(dat3$Quality == ""),]

#Replace missing Quality w/ NA
dat3[which(dat3$Quality == ""), "Quality"]<- NA


### Viz in bayesmove Shiny app ###

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




### Create summary maps per region and age class ###

states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))

north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')
south.am<- ne_countries(scale = 10, continent = "South America", returnclass = 'sf')
asia<- ne_countries(scale = 10, continent = "Asia", returnclass = 'sf')
world<- ne_countries(scale = 110, returnclass = 'sf')

pal<- wes_palette('Zissou1', n = 5)[c(1,5)]

# Create ggplot2 theme for maps
theme_map <- function() {

  font<- "Arial"

  theme_bw() %+replace%    #replace elements we want to change

  theme(
    axis.text = element_text(size = 14, family = font),
    axis.title = element_text(size = 18, family = font),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
}


# GoM

ggplot() +
  geom_sf(data = north.am) +
  geom_sf(data = states) +
  geom_point(data = dat3 %>%
               filter(Region == 'GoM'), aes(Longitude, Latitude, color = Age),
             alpha = 0.3) +
  scale_color_manual(values = pal) +
  theme_map() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(-97, -76), ylim = c(18, 31))


# Brazil

ggplot() +
  geom_sf(data = south.am %>%
            filter(name == 'Brazil')) +
  geom_point(data = dat3 %>%
               filter(Region == 'Brazil'), aes(Longitude, Latitude, color = Age),
             alpha = 0.1) +
  annotate(geom = "text", label = "italic(Brasil)", x = -45, y = -10, size = 10,
           parse = TRUE) +
  scale_color_manual(values = pal) +
  theme_map() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(-50, -30), ylim = c(-30, 0))


# Qatar

ggplot() +
  geom_sf(data = asia %>%
            filter(name %in% c('Qatar','Saudi Arabia','United Arab Emirates','Bahrain',
                               'Iran'))) +
  geom_point(data = dat3 %>%
               filter(Region == 'Qatar'), aes(Longitude, Latitude, color = Age),
             alpha = 0.3) +
  annotate(geom = "text", label = "italic(Qatar)", x = 51.15, y = 25, size = 6,
           parse = TRUE) +
  annotate(geom = "text", label = "italic(KSA)", x = 50.5, y = 24.25, size = 6,
           parse = TRUE) +
  scale_color_manual(values = pal) +
  theme_map() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(50, 54), ylim = c(24, 27))


# Global View
dat3.sf<- dat3 %>%
  st_as_sf(., coords = c('Longitude','Latitude'), crs = 4326)

ggplot() +
  geom_sf(data = world) +
  geom_sf(data = dat3.sf, aes(color = Age), alpha = 0.3) +
  scale_color_manual(values = pal) +
  theme_map() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(crs = st_crs("+proj=moll"))



### Export master data.frame

write.csv(dat3, "Raw_data/Master Sat Tag Dataset.csv", row.names = FALSE)
