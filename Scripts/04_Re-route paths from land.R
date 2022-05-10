
### Re-route tracks around land using {pathroutr} ###

library(tidyverse)
library(lubridate)
library(sf)
library(sfnetworks)
library(pathroutr)
library(ptolemy)
library(furrr)
# library(future)


## Load data
dat <- read.csv("Processed_data/Processed_Cm_Tracks_SSM_4hr.csv")

dat.gom <- dat %>%
  filter(lat > 0 & lon < -80 & lon > -100)

dat.gom.sf <- dat.gom %>%
  st_as_sf(., coords = c('lon','lat'), crs = 4326) %>%
  st_transform(3395) #%>%
  # group_by(id) %>%
  # summarize(do_union = FALSE) %>%
  # st_cast("MULTILINESTRING")

## Extract high-res coastline/land layer
gom.sf <- extract_gshhg(data = dat.gom.sf, resolution = 'f')


# Viz data
ggplot() +
  geom_sf(data = gom.sf) +
  geom_sf(data = dat.gom.sf, aes(color = id)) +
  theme_bw() +
  coord_sf()



## Identify segments on land
# segs_tbl <- get_barrier_segments(dat.gom.sf, gom.sf)
# segs_tbl

## Create Visibility Graph
vis_graph <- prt_visgraph(gom.sf)
vis_graph_sf <- sfnetworks::activate(vis_graph, "edges") %>%
  st_as_sf()

ggplot() +
  geom_sf(data = gom.sf) +
  geom_sf(data = vis_graph_sf, size = 0.5) +
  theme_void()

## Re-route model predictions via shortest path method
# segs_tbl <- segs_tbl %>%
#   prt_shortpath(vis_graph)
#
# ggplot() +
#   geom_sf(data = gom.sf) +
#   geom_sf(data = segs_tbl$geometry, color = "deepskyblue3") +
#   theme_void()

plan("multisession", workers = 10)  ## to run path re-routing in parallel across IDs
track_pts_fix <- dat.gom.sf %>%
  split(.$id) %>%
  future_map(., ~{.x %>%
      prt_trim(gom.sf) %>%
      prt_reroute(gom.sf, vis_graph) %>%
      prt_update_points(.x)}) %>%
  bind_rows()

track_pts_trim <- prt_trim(dat.gom.sf, gom.sf)
track_pts_fix <- prt_reroute(track_pts_trim, gom.sf, vis_graph)
track_pts_fix2 <- prt_update_points(track_pts_fix, track_pts_trim)

track_line_fixed <- track_pts_fix2 %>%
  group_by(id) %>%
  summarise(do_union = FALSE) %>%
  st_cast('MULTILINESTRING')

ggplot() +
  geom_sf(data = gom.sf) +
  geom_sf(data = track_line_fixed %>%
            filter(id == 172678), aes(color = id)) +
  theme_void()
