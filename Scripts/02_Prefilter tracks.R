
##########################################################################
### Prefilter tracks for analysis by continuous-time state-space model ###
##########################################################################


library(tidyverse)
library(lubridate)
library(sf)
library(rnaturalearth)
library(plotly)



## Load turtle data
dat<- read.csv("Raw_data/Master Sat Tag Dataset.csv") %>%
  mutate(Date = as_datetime(Date))

deploy <- read.csv("Raw_data/Turtle deploy metadata.csv") %>%
  mutate(Deploy_Date = parse_date_time(Deploy_Date, orders = "mdyHM", truncated = 2))


## Load land spatial layers
states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))
north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')

brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')

mid.east<- ne_countries(scale = 10, country = c("Qatar", "Iran", "Bahrain", "Saudi Arabia",
                                                "United Arab Emirates"), returnclass = 'sf')


## Remove all LC Z locs
dat2<- dat %>%
  rename(date = Date) %>%
  filter(Quality != "Z")  #remove Z LCs




## Remove observations before deployment; easiest to do this per region and source

###############
# Fuentes-GoM #
###############

deploy.fuentes.gom <- deploy %>%
  filter(Source == "Fuentes" & Region == "GoM") %>%
  split(.$PTT)

fuentes.gom<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'GoM') %>%
  split(.$Ptt) %>%
  map2(.x = ., .y = deploy.fuentes.gom,
       ~{.x %>%
           filter(., date >= .y$Deploy_Date)}) %>%
  bind_rows()

ggplot() +
  geom_sf(data = north.am) +
  geom_sf(data = states) +
  geom_path(data = fuentes.gom, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-85, -70), ylim = c(24, 31))

# Filter by study extent
fuentes.gom2<- fuentes.gom %>%
  filter(Latitude < 30 & Longitude > -84 & Longitude < -74)


# Check if data have been sufficiently filtered
fuentes.gom2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)
#potential issues w/ 159774 and 159783

ggplot() +
  geom_path(data = fuentes.gom2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Create interactive map
plotly::ggplotly(
  ggplot() +
    geom_sf(data = north.am) +
    geom_sf(data = states) +
    geom_path(data = fuentes.gom2, aes(Longitude, Latitude, color = factor(Ptt)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(-85, -70), ylim = c(24, 31))
)




# Filter by longitude for Crystal River pre-deployment locations
cr.filt.dates <- fuentes.gom2 %>%
  filter(str_detect(Ptt, "^159")) %>%
  group_by(Ptt) %>%
  filter(Longitude < -82.7) %>%
  slice(1) %>%
  dplyr::select(date) %>%
  group_split()

fuentes.gom3 <- fuentes.gom2 %>%
  split(.$Ptt)

fuentes.gom3[1:5] <- fuentes.gom3[1:5] %>%
  map2(.x = ., .y = cr.filt.dates, ~{
    .x %>%
      filter(date >= .y$date)
  })

# Convert from list to data.frame and remove aberrant locs from end of 159774 track
fuentes.gom3 <- bind_rows(fuentes.gom3) %>%
  filter(!(Ptt == '159774' & date > as_datetime('2016-07-04 12:00:05')))

# Check if data have been sufficiently filtered
ggplot() +
  geom_path(data = fuentes.gom3, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

fuentes.gom3 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

#original length is 8556; after filetering w/ deploy dates, new length is 8527


##############
# Lamont-GoM #
##############

deploy.lamont.gom <- deploy %>%
  filter(Source == 'Lamont' & Region == 'GoM')

lamont.gom<- dat2 %>%
  filter(Source == 'Lamont' & Region == 'GoM') %>%
  split(.$Ptt)

deploy.ind <- which(names(lamont.gom) %in% deploy.lamont.gom$PTT)

lamont.gom[deploy.ind] <- lamont.gom[deploy.ind] %>%
  map2(.x = ., .y = deploy.lamont.gom %>% split(.$PTT),
       ~{.x %>%
           filter(., date >= .y$Deploy_Date)})

lamont.gom <- bind_rows(lamont.gom)

ggplot() +
  geom_sf(data = north.am) +
  geom_sf(data = states) +
  geom_path(data = lamont.gom, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-96, -70), ylim = c(20, 31))

# Filter by study extent
lamont.gom2<- lamont.gom %>%
  filter(Latitude < 31 & Longitude > -96 & Longitude < -81)


# Check if data have been sufficiently filtered
lamont.gom2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)
#Likely undeployed locs in PTTs 142658, 142659, 161459, 172677, 175692

ggplot() +
  geom_path(data = lamont.gom2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Check possible issues w/ 161459
lamont.gom2 %>%
  filter(Ptt == 161459) %>%
  slice(1:100)  #remove obs before 2017-08-03 (roughly)


# Filter by date
lamont.gom3 <- lamont.gom2 %>%
  filter(!(Ptt %in% 161459 & date < "2017-08-03"))

# Check if data have been sufficiently filtered
ggplot() +
  geom_path(data = lamont.gom3, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Create interactive map
plotly::ggplotly(
  ggplot() +
    geom_sf(data = north.am) +
    geom_sf(data = states) +
    geom_path(data = lamont.gom3, aes(Longitude, Latitude, color = factor(Ptt)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(-96, -70), ylim = c(20, 31))
)


lamont.gom3 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

#original has length 19512; filtered by deploy date has length 19444





#############
# Sasso-GoM #
#############

sasso.gom<- dat2 %>%
  filter(Source == 'Sasso' & Region == 'GoM')

ggplot() +
  geom_sf(data = north.am) +
  geom_sf(data = states) +
  geom_path(data = sasso.gom, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-90, -78), ylim = c(18, 30))

# Filter by study extent
sasso.gom2<- sasso.gom %>%
  filter(Latitude > 20 & Latitude < 29 & Longitude > -88 & Longitude < -80)


# Check if data have been sufficiently filtered
sasso.gom2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

ggplot() +
  geom_path(data = sasso.gom2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Create interactive map
plotly::ggplotly(
  ggplot() +
    geom_sf(data = north.am) +
    geom_sf(data = states) +
    geom_path(data = sasso.gom2, aes(Longitude, Latitude, color = factor(Ptt)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(-90, -78), ylim = c(18, 30))
)






###############
# Fuentes-FDN #
###############

deploy.fuentes.fdn <- deploy %>%
  filter(str_detect(PTT, '^205')) %>%
  split(.$PTT)

fuentes.fdn<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'Brazil' & Age == 'Adult') %>%
  split(.$Ptt) %>%
  map2(.x = ., .y = deploy.fuentes.fdn,
       ~{.x %>%
           filter(., date >= .y$Deploy_Date)}) %>%
  bind_rows()

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = fuentes.fdn, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -1))

# Check if data have been sufficiently filtered
fuentes.fdn %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)


# don't need to remove any (pre-deployment) locations

#original has length 24162; filtered by deploy date has same length





##################
# Fuentes-Parana #
##################

deploy.fuentes.par <- deploy %>%
  filter(str_detect(PTT, "^160")) %>%
  split(.$PTT)

fuentes.par<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'Brazil' & Age == 'Juv') %>%
  split(.$Ptt)

deploy.ind <- which(names(fuentes.par) %in% names(deploy.fuentes.par))

fuentes.par[deploy.ind] <- fuentes.par[deploy.ind] %>%
  map2(.x = ., .y = deploy.fuentes.par,
       ~{.x %>%
           filter(., date >= .y$Deploy_Date)})

fuentes.par <- bind_rows(fuentes.par)


ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = fuentes.par, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-49, -38), ylim = c(-27, -15))

# Check if data have been sufficiently filtered
fuentes.par %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

#remove points before 2016-05-22 12:00 for PTT 160105_1
#remove points before 2018-03-16 for PTT 160105_2


# Filter by date
fuentes.par2 <- fuentes.par %>%
  filter(!(Ptt %in% '160105_1' & date < "2016-05-22 12:00:00")) %>%
  filter(!(Ptt %in% '160105_2' & date < "2018-03-16"))

# Check if data have been sufficiently filtered
fuentes.par2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

#original has length of 7530; filtered by deploy date has same length



################
# Domit-Parana #
################

deploy.domit.par <- deploy %>%
  filter(Source == "Domit") %>%
  split(.$PTT)

domit.par<- dat2 %>%
  filter(Source == 'Domit' & Region == 'Brazil') %>%
  split(.$Ptt) %>%
  map2(.x = ., .y = deploy.domit.par,
       ~{.x %>%
           filter(., date > .y$Deploy_Date)}) %>%
  bind_rows()

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = domit.par, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-49, -46), ylim = c(-27, -23))



# Check if data have been sufficiently filtered
domit.par %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)
#remove first 2 points for each PTT since these are generally associated w/ same same location on land


domit.par2 <- domit.par %>%
  group_by(Ptt) %>%
  slice(3:n())

ggplot() +
  geom_path(data = domit.par2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Create interactive map
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = domit.par2, aes(Longitude, Latitude, color = factor(Ptt)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(-49, -46), ylim = c(-27, -23))
)

domit.par2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

#original has length 7573; filtered by deploy date has length 7488




##################
# Marshall-Qatar #
##################

marshall.qatar<- dat2 %>%
  filter(Source == 'Marshall')

ggplot() +
  geom_sf(data = mid.east) +
  geom_path(data = marshall.qatar, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(50, 53), ylim = c(24, 27))

# Filter by study extent
marshall.qatar2<- marshall.qatar %>%
  filter(Longitude < 53)


# Check if data have been sufficiently filtered
marshall.qatar2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)

ggplot() +
  geom_path(data = marshall.qatar2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Create interactive map
plotly::ggplotly(
  ggplot() +
    geom_sf(data = mid.east) +
    geom_path(data = marshall.qatar2, aes(Longitude, Latitude, color = factor(Ptt)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(50, 53), ylim = c(24, 27))
)




######################################################################
### Account for location error w/ continuous-time SSM based on CRW ###
######################################################################

# merge all filtered datasets back together
dat3 <- rbind(fuentes.gom3,
              lamont.gom3,
              sasso.gom2,
              fuentes.fdn,
              fuentes.par2,
              domit.par2,
              marshall.qatar2)


# Only keep IDs tracked w/ > 30 observations
dat3 %>%
  group_by(Ptt) %>%
  count() %>%
  data.frame()
#all IDs besides 161639 have >=30 obs

dat4 <- dat3 %>%
  group_split(Ptt) %>%
  discard( ~ nrow(.x) < 30) %>%
  bind_rows()




## Export processed data

# write.csv(dat4, "Processed_data/Prefiltered_Cm_Tracks.csv", row.names = FALSE)  #add info here to save your data


