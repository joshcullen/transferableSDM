
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
dat<- read.csv("Processed_data/Prefiltered_Cm_Tracks.csv") %>%
  mutate(date = as_datetime(date))


## Load land spatial layers
states<- ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(!name %in% c('Alaska', 'Hawaii'))
north.am<- ne_countries(scale = 10, continent = "North America", returnclass = 'sf')

brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')

mid.east<- ne_countries(scale = 10, country = c("Qatar", "Iran", "Bahrain", "Saudi Arabia",
                                                "United Arab Emirates"), returnclass = 'sf')



## Assign 'Quality' equal to "G" for GPS locs
dat$Quality<- ifelse(dat$Type == 'FastGPS', 'G', dat$Quality)



## Reformat data for model from foieGras package
dat2<- dat %>%
  rename(id = Ptt, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)


## Calculate avg speed per individual

speed <- dat2 %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326, remove = FALSE) %>%
  st_transform(crs = 3395) %>%
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  bayesmove::prep_data(coord.names = c('x','y'), id = "id") %>%
  mutate(speed = step/dt) %>%
  group_by(id) %>%
  summarize(mean.speed = median(speed, na.rm = TRUE))

# Mean speed across IDs
mean(speed$mean.speed) * 3600  #turtles typically can travel ~1 km in 1 hr



## Run model (w/o predicting; no track regularization) for all IDs at once since not hierarchical model
tic()
fit.crw <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = 4, control = ssm_control(verbose = 1))
toc()  #took 5 min where time.step = 4

# Check convergence and positive-definite Hessian matrices
oo <- 1
for (i in 1:9) {
  print(fit.crw[oo:(oo+9),])
  oo<- oo + 10
}
## looks good



# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit.crw, what = "predicted", type = 1, ask = TRUE)
plot(fit.crw, what = "predicted", type = 2, ask = TRUE)

map(x = fit.crw, y = NULL, what = "fitted")



## Grab results and plot

res.crw<- grab(fit.crw, what = "predicted", as_sf = FALSE)




# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat2 %>%
              filter(str_detect(id, "^205")),
            aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res.crw %>%
              filter(str_detect(id, "^205")),
            aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -1))


# Compare raw tracks vs fitted tracks (for juveniles tagged near Paranagua)
# ggplot() +
#   geom_sf(data = brazil) +
#   geom_path(data = dat2 %>%
#               filter(!str_detect(id, "^205")),
#             aes(lon, lat, group = id), color = 'black') +  #raw tracks
#   geom_path(data = res.crw %>%
#               filter(!str_detect(id, "^205")),
#             aes(lon, lat, group = id), color = "blue") +  #modeled tracks
#   theme_bw() +
#   facet_wrap(~id) +
#   coord_sf(xlim = c(-49, -38), ylim = c(-27, -15))


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

# write.csv(res.crw, "Processed_data/Processed_Cm_Tracks_SSM_4hr.csv", row.names = FALSE)







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
