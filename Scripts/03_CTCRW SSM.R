
### Fit CTCRW SSM to turtle data via {aniMotum} ###

library(tidyverse)
library(lubridate)
library(aniMotum)
library(sf)
library(sfarrow)
library(tictoc)
library(plotly)
library(trelliscopejs)
library(patchwork)
library(ptolemy)
library(fuzzyjoin)

source('Scripts/helper functions.R')


#################
### Load data ###
#################

dat <- read_csv("Processed_data/Prefiltered_Cm_Tracks.csv")

glimpse(dat)
summary(dat)


##################################################
### Wrangle data for analysis using {aniMotum} ###
##################################################

# Convert all 'Quality' values to "G" for FastGPS data and 'Date' to datetime format
dat <- dat %>%
  mutate(Quality = ifelse(Type == 'FastGPS', 'G', Quality))


# Rename columns for {foieGras}
dat2<- dat %>%
  rename(id = Ptt, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)  #reorders and subsets the columns

glimpse(dat2)




# Define bouts per PTT where gap > 3 days (for splitting fitted tracks)
int_tbl <- dat2 %>%
  rename(datetime = date) %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x, bout = cu_add_gaps(x = ., gap = 3, time_unit = "days") %>%
                          pull(bout_id))) %>%
  purrr::map(., function (x)
  {
    date <- NULL
    x <- x %>%
      group_by(.data[["bout"]]) %>%
      summarize(start = min(datetime), end = max(datetime))
  })


##################################################################
### Inspect time steps of transmissions for making predictions ###
##################################################################

# Calculate time step by ID
tmp <- dat2 %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x,
                        dt = difftime(c(date[-1], NA),
                                      date,
                                      units = "hours") %>%
                          as.numeric())
  ) %>%
  bind_rows()


ggplot(tmp, aes(date, dt)) +
  geom_point() +
  theme_bw() +
  facet_trelliscope(~id, nrow = 5, ncol = 5, scales = "free_x")

# Determine primary time step by ID
tmp %>%
  group_by(id) %>%
  summarize(mean = mean(dt, na.rm = TRUE),
            median = median(dt, na.rm = TRUE)) %>%
  data.frame()
# Mean/median time step is ~1-2 hrs for most IDs




#################
#### Run SSM ####
#################

# Estimate 'true' locations regularized at 4 hr time step
set.seed(2022)
tic()
fit_crw <- fit_ssm(dat2, vmax = 3, model = "crw", time.step = 4,
                          control = ssm_control(verbose = 1))
toc()  #took 2.5 min

print(fit_crw, n = nrow(fit_crw))  #all indiv. models converged
summary(fit_crw)

# Viz the filtered outliers (black x's), raw observations (blue points), and estimated locations (orange points), along with associated uncertainty (orange shading)
plot(fit_crw, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)



#############################
### Process fitted tracks ###
#############################

# Grab results and create data.frame
res_crw <- grab(fit_crw, what = "predicted")



# Separate tracks by region
gom.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'GoM',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)

br.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'Brazil',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)

qa.tracks <- res_crw %>%
  filter(id %in% unique(dat[dat$Region == 'Qatar',]$Ptt)) %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395)


# Generate hi-res coastline spatial layers
tic()
gom.sf <- ptolemy::extract_gshhg(data = gom.tracks, resolution = 'f', buffer = 500000)
toc()  #takes 1.5 min to run

tic()
br.sf <- ptolemy::extract_gshhg(data = br.tracks, resolution = 'f', buffer = 200000)
toc()  #takes 1 min to run

tic()
qa.sf <- ptolemy::extract_gshhg(data = qa.tracks, resolution = 'f', buffer = 100000)
toc()  #takes 1 min to run



# Defined bout periods per observation
gom.tracks2 <- gom.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'GoM',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()

br.tracks2 <- br.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'Brazil',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()

qa.tracks2 <- qa.tracks %>%
  split(.$id) %>%
  map2(.,
       int_tbl[names(int_tbl) %in% unique(dat[dat$Region == 'Qatar',]$Ptt)],
       ~fuzzyjoin::fuzzy_left_join(x = .x, y = .y,
                                   by = c(date = "start", date = "end"),
                                   match_fun = list(`>=`, `<=`))) %>%
  bind_rows()






# Remove location estimates from large gaps (> 3 days) and convert to DF from sf object
gom.tracks3 <- gom.tracks2 %>%
  drop_na(bout) %>%  #remove predictions that fall w/in 3-day time gaps
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()

br.tracks3 <- br.tracks2 %>%
  drop_na(bout) %>%  #remove predictions that fall w/in 3-day time gaps
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()

qa.tracks3 <- qa.tracks2 %>%
  drop_na(bout) %>%  #remove predictions that fall w/in 3-day time gaps
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  split(.$id) %>%
  purrr::map(~add_track_gaps(.x, "bout")) %>%  #add single rows of NA coords to create gaps when mapping
  bind_rows()




#############################
### Inspect mapped tracks ###
#############################

# GoM
ggplotly(
  ggplot() +
    geom_sf(data = gom.sf) +
    geom_path(data = gom.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = gom.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")


# Brazil
ggplotly(
  ggplot() +
    geom_sf(data = br.sf) +
    geom_path(data = br.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = br.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")


# Qatar
ggplotly(
  ggplot() +
    geom_sf(data = qa.sf) +
    geom_path(data = qa.tracks3, aes(x, y, group = id, color = id)) +  #modeled tracks
    theme_bw()
)

ggplot(data = qa.tracks3, aes(x, y, color = factor(bout))) +
  geom_path() +
  theme_bw() +
  facet_trelliscope(~ id, nrow = 3, ncol = 3, scales = "free")



### Map tracks by age class ###

age.df <- dat %>%
  dplyr::select(Ptt, Age) %>%
  rename(id = Ptt) %>%
  distinct()

dat3 <- rbind(gom.tracks3, br.tracks3, qa.tracks3) %>%
  drop_na(x, y) %>%
  st_as_sf(., coords = c('x','y'), crs = 3395) %>%
  st_transform(4326) %>%
  left_join(., age.df, by = "id") %>%
  mutate(Age = factor(Age, levels = c("Juv","Adult")))



p.gom <- ggplot() +
  geom_sf(data = gom.sf |>
            st_transform(4326)) +
  geom_sf(data = dat3 |>
            st_transform(4326), aes(color = Age), size = 0.5, alpha = 0.5) +
  scale_color_brewer("Life Stage", palette = "Dark2") +
  annotate(geom = "text", label = "Gulf of\nMexico", fontface = "italic",
           size = 8, x = -92, y = 25) +
  geom_text(aes(x = -97.5, y = 30, label = "(a)"), size = 10, fontface = "bold") +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = "white")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_sf(xlim = c(-98,-77), ylim = c(20,30.5))


p.br <- ggplot() +
  geom_sf(data = br.sf |>
            st_transform(4326)) +
  geom_sf(data = dat3 |>
            st_transform(4326), aes(color = Age), size = 0.5, alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  annotate(geom = "text", label = "Brazil", fontface = "italic", size = 8, x = -43, y = -8) +
  geom_text(aes(x = -47, y = -2.5, label = "(b)"), size = 10, fontface = "bold") +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = "white")
  ) +
  coord_sf(xlim = c(-49,-31), ylim = c(-26,-2))


p.qa <- ggplot() +
  geom_sf(data = qa.sf |>
            st_transform(4326)) +
  geom_sf(data = dat3 |>
            st_transform(4326), aes(color = Age), size = 0.5, alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  annotate(geom = "text", label = "Qatar", fontface = "italic", size = 8,
           x = 51.0, y = 25.3) +
  geom_text(aes(x = 50.5, y = 27.1, label = "(c)"), size = 10, fontface = "bold") +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = "white")
  ) +
  coord_sf(xlim = c(50.25,53), ylim = c(23.8,27.2))


p.gom / (p.br + p.qa)

ggsave("Tables_Figs/Figure S1.png", width = 5, height = 7, units = "in", dpi = 400)




###############################################
### Export fitted tracks and spatial layers ###
###############################################

# tracks
write.csv(gom.tracks3, "Processed_data/Processed_GoM_Cm_Tracks_SSM_4hr_aniMotum.csv", row.names = FALSE)
write.csv(br.tracks3, "Processed_data/Processed_Brazil_Cm_Tracks_SSM_4hr_aniMotum.csv", row.names = FALSE)
write.csv(qa.tracks3, "Processed_data/Processed_Qatar_Cm_Tracks_SSM_4hr_aniMotum.csv", row.names = FALSE)


# spatial layers
st_write_parquet(gom.sf, 'Environ_data/GoM_land.parquet')
st_write_parquet(br.sf, 'Environ_data/Brazil_land.parquet')
st_write_parquet(qa.sf, 'Environ_data/Qatar_land.parquet')
