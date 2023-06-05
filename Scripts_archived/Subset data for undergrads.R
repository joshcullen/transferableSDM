

### Subset data for undergrads (Spring 2022) ###

## This is an abridged script following steps used in the 'Clean and explore tracks.R' file.

library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(bayesmove)
library(wesanderson)


### Load data and check column names ###
dat<- read.csv("Raw_data/Master Sat Tag Dataset.csv")



### Prep data ###
dat$Date<- as_datetime(dat$Date)



### Load spatial data ###

states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))

north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')
south.am<- ne_countries(scale = 10, continent = "South America", returnclass = 'sf')

pal<- wes_palette('Zissou1', n = 5)[c(1,5)]


### Subset data ###

# Colby Anderson
## Just sent the entire CSV for the Qatar data


# Chanti Max
chanti.df<- dat %>%
  filter(Region == 'Brazil' & Age == 'Adult')

ggplot() +
  geom_sf(data = south.am %>%
            filter(name == 'Brazil')) +
  geom_path(data = chanti.df, aes(Longitude, Latitude, color = Ptt), alpha = 0.7) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(-40, -30), ylim = c(-10, 0))

# write.csv(chanti.df, "../../Undergrad Docs/Independent Research Projects/ChantiMax/Brazil_adult_Cm_data.csv", row.names = FALSE)


# Margaret Rivas
## Just copied same file from Chanti's folder

# Natalie Tosto

nat.df<- dat %>%
  filter(Region == 'GoM' & Source == 'Fuentes')

ggplot() +
  geom_sf(data = north.am) +
  geom_sf(data = states) +
  geom_point(data = nat.df, aes(Longitude, Latitude, color = Ptt),
             alpha = 0.3) +
  scale_color_viridis_d() +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(-85, -76), ylim = c(24, 31))

# write.csv(nat.df, "../../Undergrad Docs/Independent Research Projects/NatalieTosto/CrystalRiver_Bimini_Cm_data.csv", row.names = FALSE)


# Lauren Hearn

lauren.df<- dat %>%
  filter(Region == 'Brazil')

ggplot() +
  geom_sf(data = south.am) +
  geom_path(data = lauren.df, aes(Longitude, Latitude, color = Ptt), alpha = 0.7) +
  scale_color_viridis_d() +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.position = "top",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_sf(xlim = c(-50, -30), ylim = c(-30, 0))

# write.csv(lauren.df, "../../Undergrad Docs/Independent Research Projects/LaurenHearn/Brazil_Cm_data.csv", row.names = FALSE)
