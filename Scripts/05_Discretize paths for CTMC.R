
### Create gridded paths for turt tracks ###

library(tidyverse)
library(lubridate)
library(sf)
library(crawlUtils)
library(furrr)
library(future)
library(vroom)
library(terra)
library(cmocean)
library(tictoc)
library(progress)

source("Scripts/make_hex_grid2.R")
source("Scripts/make_disc_path2.R")


### Load turtle tracks ###
dat <- vroom("Processed_data/Imputed_Cm_Tracks_SSM_30min.csv", delim = ",")
dat.orig <- vroom("Processed_data/Prefiltered_Cm_Tracks.csv", delim = ",") %>%
  filter(Region == "GoM") %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(3395)
quad_times <- read.csv("Environ_data/quadrature_pts.csv") %>%
  dplyr::pull(date) %>%
  as_date()

glimpse(dat)
glimpse(dat.orig)

dat.sf.p <- dat %>%
  # filter(ptt == 181800) %>%
  # filter(ptt %in% c(159776, 175692, 181796, 181800, 181807)) %>%
  st_as_sf(., coords = c('mu.x', 'mu.y'), crs = 3395)



### Create hexgrids based on imputed tracks ###
cellsize <- hex_size(area = 4000^2)  #16 km^2

plan("multisession", workers = availableCores() - 2)
tic()
hex_grid <- dat.sf.p %>%
  split(.$ptt) %>%
  future_map(~hex_grid_sfsample(sf_tracks = ., cellsize = cellsize, buffer = 15000)) %>%
  set_names(unique(dat.sf.p$ptt))
toc()
plan("sequential")  #return to single core
# takes 5 min to run on desktop and laptop

# dat.sf.l <- dat.sf.p %>%
#   group_by(ptt, rep) %>%
#   summarize(do_union = FALSE) %>%
#   st_cast("MULTILINESTRING")


# Download and clip high-res land spatial layer
# gom.sf <- ptolemy::extract_gshhg(data = dat.orig, resolution = 'f', buffer = 50000)

# bb = st_bbox(hex_grid$`181800`$poly)
# plotly::ggplotly(
#   ggplot() +
#     geom_sf(data = gom.sf, size = 0.2) +
#     geom_sf(data = hex_grid$`181800`$poly, fill = "transparent", size = 0.25) +
#     geom_sf(data = dat.sf.l %>%
#               filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
#     scale_color_viridis_d("Imputation", option = "turbo") +
#     # geom_sf(data = dat.sf.p %>%
#     #           filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
#     geom_sf(data = dat.orig %>%
#                  filter(Ptt == 181800), size = 1) +
#     coord_sf(xlim = c(bb["xmin"]-50000, bb["xmax"]+50000), ylim = c(bb["ymin"]-50000, bb["ymax"]+50000)) +
#     theme_bw() +
#     theme(panel.grid = element_blank())
# )


if(!dir.exists("Plots")) dir.create("Plots")



### Export plots of hexgrids w/ imputed tracks and original points ###

# PTTs <- unique(dat.sf.l$ptt)
#
# for (i in 1:length(PTTs)) {
#   message("plotting ", PTTs[i], " ...")
#
#   lines <- dat.sf.l %>%
#     filter(ptt == PTTs[i])
#   pts <- dat.sf.p %>%
#     filter(ptt == PTTs[i])
#   orig.pts <- dat.orig %>%
#     filter(Ptt == PTTs[i])
#   hexes <- hex_grid[[i]]$poly
#   bb <- st_bbox(hexes)
#
#   p_data <- ggplot() +
#     geom_sf(data = gom.sf, size = 0.2) +
#     geom_sf(data = hexes, fill = "transparent", size = 0.25) +
#     geom_sf(data = lines, aes(color = rep), size = 0.25, alpha = 0.6) +
#     scale_color_viridis_d(option = "turbo", guide = "none") +
#     # geom_sf(data = pts, aes(color = rep), size = 0.25, alpha = 0.6) +
#     geom_sf(data = orig.pts, size = 1) +
#     coord_sf(xlim = c(bb["xmin"] - 10000, bb["xmax"] + 10000),
#              ylim = c(bb["ymin"] - 10000, bb["ymax"] + 10000)) +
#     ggtitle(paste("PTT", PTTs[i])) +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           title = element_text(size = 18, face = "bold"))
#
#   nm = paste0("Plots/p_", PTTs[i],"_data.pdf")
#   ggsave(nm, p_data, dpi = "retina", width = 7, height = 7, units = 'in')
#
# }



### Discretize imputation paths for each turtle ###
plan("multisession", workers = availableCores() - 2)
tic()
disc_path <- make_disc_path(sims = dat.sf.p, hex_grid = hex_grid, quad_times = quad_times)
toc()
plan("sequential")  #return to single core


### Export products ###

# save(dat, disc_path, hex_grid, gom.sf, file = "Data_products/discretized_tracks.RData")


