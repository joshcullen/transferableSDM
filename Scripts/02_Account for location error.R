
##################################################################
### Account for location error in tracks via SSM in {foieGras} ###
##################################################################


library(tidyverse)
library(lubridate)
library(sf)
library(foieGras)
library(rnaturalearth)
library(tictoc)
library(plotly)
library(future)  #needed to properly run foieGras::osar() in parallel



## Load turtle data
dat<- read.csv("Raw_data/Master Sat Tag Dataset.csv") %>%
  mutate(Date = as_datetime(Date))


## Load land spatial layers
states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))
north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')

brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')

mid.east<- ne_countries(scale = 10, country = c("Qatar", "Iran", "Bahrain", "Saudi Arabia",
                                                "United Arab Emirates"), returnclass = 'sf')


## Assign 'Quality' equal to "G" for GPS locs
dat$Quality<- ifelse(dat$Type == 'FastGPS', 'G', dat$Quality)

## Generate time steps for data (along w/ SL, TA, and NSD)
dat2<- dat %>%
  rename(date = Date) %>%
  filter(Quality != "Z")  #remove Z LCs







## Remove observations before deployment; easiest to do this per region and source

###############
# Fuentes-GoM #
###############

fuentes.gom<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'GoM')

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

fuentes.gom2 %>%
  rename(id = Ptt, x = Longitude, y = Latitude) %>%
  dplyr::select(id, date, x, y) %>%
  bayesmove::shiny_tracks(epsg = 4326)


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




##############
# Lamont-GoM #
##############

lamont.gom<- dat2 %>%
  filter(Source == 'Lamont' & Region == 'GoM')

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
ggplot() +
  geom_path(data = lamont.gom2, aes(date, Longitude, color = factor(Ptt))) +
  theme_bw() +
  facet_wrap(~ Ptt, scales = "free")

# Check possible issues w/ 142658, 142659, 172677, and 175692
lamont.gom2 %>%
  filter(Ptt == 142658) %>%
  slice(1:50)  #remove obs before 2017-09-22

lamont.gom2 %>%
  filter(Ptt == 142659) %>%
  slice(1:50)  #remove obs before 2017-11-03

lamont.gom2 %>%
  filter(Ptt == 172677) %>%
  slice(1:10)  #remove obs before 2019-01-01

lamont.gom2 %>%
  filter(Ptt == 175692) %>%
  slice(1:50)  #remove obs before 2019-01-01


# Filter by date
lamont.gom3 <- lamont.gom2 %>%
  filter(!(Ptt %in% c(142658, 142659) & date < "2017-09-22")) %>%
  filter(!(Ptt %in% c(172677, 175692) & date < "2019-01-01"))

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

fuentes.fdn<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'Brazil' & Age == 'Adult')

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = fuentes.fdn, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -1))

# don't need to remove any (pre-deployment) locations




##################
# Fuentes-Parana #
##################

fuentes.par<- dat2 %>%
  filter(Source == 'Fuentes' & Region == 'Brazil' & Age == 'Juv')

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = fuentes.par, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-49, -38), ylim = c(-27, -15))

# don't need to remove any (pre-deployment) locations





################
# Domit-Parana #
################

domit.par<- dat2 %>%
  filter(Source == 'Domit' & Region == 'Brazil')

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = domit.par, aes(Longitude, Latitude, color = factor(Ptt))) +
  theme_bw() +
  coord_sf(xlim = c(-49, -46), ylim = c(-27, -23))

# Filter by study extent
domit.par2<- domit.par %>%
  filter(Longitude < -46)


# Check if data have been sufficiently filtered
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
              fuentes.par,
              domit.par2,
              marshall.qatar2)


# Only keep IDs tracked w/ > 30 observations
dat4 <- dat3 %>%
  group_split(Ptt) %>%
  discard( ~ nrow(.x) < 30) %>%
  bind_rows()

dat3 %>%
  group_by(Ptt) %>%
  count() %>%
  data.frame()
#all IDs have >=30 obs



## Reformat data for model from foieGras package
dat5<- dat4 %>%
  rename(id = Ptt, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)



## Run model (w/o predicting; no track regularization) for all IDs at once since not hierarchical model
tic()
fit.crw <- fit_ssm(dat5, vmax = 3, model = "crw", time.step = NA, control = ssm_control(verbose = 1))
toc()  #took 2 min where time.step = NA

# Check convergence and positive-definite Hessian matrices
oo <- 1
for (i in 1:9) {
  print(fit.crw[oo:(oo+9),])
  oo<- oo + 10
}
## looks good



# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit.crw, what = "fitted", type = 1, ask = TRUE)
plot(fit.crw, what = "fitted", type = 2, ask = TRUE)

fmap(x = fit.crw, y = NULL, what = "fitted")



## Grab results and plot

res.crw<- grab(fit.crw, what = "fitted", as_sf = FALSE)




# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat5 %>%
              filter(str_detect(id, "^205")),
            aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res.crw %>%
              filter(str_detect(id, "^205")),
            aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -1))


# Compare raw tracks vs fitted tracks (for juveniles tagged near Paranagua)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat5 %>%
              filter(!str_detect(id, "^205")),
            aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res.crw %>%
              filter(!str_detect(id, "^205")),
            aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-49, -38), ylim = c(-27, -15))


# Viz modeled tracks together
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = res.crw, aes(lon, lat, group = id, color = id)) +
  scale_color_viridis_d(option = "inferno") +
  theme_bw() +
  coord_sf(xlim = c(-49, -32), ylim = c(-27, -1))



# Interactively viz time-series of tracks
res.crw %>%
  dplyr::select(-c(x, y, s.se)) %>%
  filter(lon < 0 & lat > 0) %>%
  rename(x = lon, y = lat) %>%
  dplyr::select(id, date, x, y) %>%
  data.frame() %>%  #current version of shiny_tracks won't work w/ tibble format or when a column has all NAs
  bayesmove::shiny_tracks(., 4326)






## Export results

# write.csv(res.crw, "Processed_data/Processed_Cm_Tags_SSM.csv", row.names = FALSE)  #add info here to save your data







## Filter data only for GoM (Fuentes and Lamont[juv only]) for NOAA
id.noaa <- dat4 %>%
  filter(Region == 'GoM', Source %in% c('Lamont','Fuentes'), Age == 'Juv') %>%
  filter(!str_detect(Ptt, "^169")) %>%
  dplyr::select(Ptt) %>%
  distinct() %>%
  unlist()

dat.noaa <- res.crw %>%
  filter(id %in% id.noaa)

plotly::ggplotly(
  ggplot() +
    geom_sf(data = north.am) +
    geom_sf(data = states) +
    geom_path(data = dat.noaa, aes(lon, lat, color = factor(id)), alpha = 0.6) +
    theme_bw() +
    coord_sf(xlim = c(-96, -70), ylim = c(20, 31))
)

# Add metadata column for "source"
dat.noaa <- dat.noaa %>%
  mutate(source = case_when(str_detect(id, "^159") ~ "Fuentes",
                            TRUE ~ as.character("Lamont")
  ),
  age = "Juv",
  .before = id)

dat.noaa2 <- dat.noaa %>%
  dplyr::select(source, age, lon, lat)

#export
# write.csv(dat.noaa2, "Processed_data/Fuentes_Lamont_Cm_juv_tracks.csv", row.names = FALSE)
