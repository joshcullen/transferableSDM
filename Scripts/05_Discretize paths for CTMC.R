
### Create gridded paths for turt tracks ###
# Original code modified from Supplement of Brennan et al 2018 (https://static-content.springer.com/esm/art%3A10.1007%2Fs10980-018-0642-z/MediaObjects/10980_2018_642_MOESM3_ESM.html)

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

source("Scripts/make_hex_grid2.R")


## Load turtle tracks
dat <- vroom("Processed_data/Imputed_Cm_Tracks_SSM_30min.csv", delim = ",")
dat.orig <- vroom("Processed_data/Prefiltered_Cm_Tracks.csv", delim = ",") %>%
  filter(Region == "GoM") %>%
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(3395)

glimpse(dat)
glimpse(dat.orig)

dat.sf.p <- dat %>%
  # filter(ptt == 181800) %>%
  # filter(ptt %in% c(159776, 175692, 181796, 181800, 181807)) %>%
  st_as_sf(., coords = c('mu.x', 'mu.y'), crs = 3395)

cellsize <- hex_size(area = 4000^2)

plan("multisession", workers = availableCores() - 2)
tic()
hex_grid <- dat.sf.p %>%
  split(.$ptt) %>%
  future_map(~hex_grid_sfsample(sf_tracks = ., cellsize = cellsize, buffer = 15000)) %>%
  set_names(unique(dat.sf.p$ptt))
toc()
# takes 5 min to run

dat.sf.l <- dat.sf.p %>%
  group_by(ptt, rep) %>%
  summarize(do_union = FALSE) %>%
  st_cast("MULTILINESTRING")


# Download and clip high-res land spatial layer
gom.sf <- ptolemy::extract_gshhg(data = dat.orig, resolution = 'f', buffer = 50000)

bb = st_bbox(hex_grid$`181800`$poly)
plotly::ggplotly(
  ggplot() +
    geom_sf(data = gom.sf, size = 0.2) +
    geom_sf(data = hex_grid$`181800`$poly, fill = "transparent", size = 0.25) +
    geom_sf(data = dat.sf.l %>%
              filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
    scale_color_viridis_d("Imputation", option = "turbo") +
    # geom_sf(data = dat.sf.p %>%
    #           filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
    geom_sf(data = dat.orig %>%
                 filter(Ptt == 181800), size = 1) +
    coord_sf(xlim = c(bb["xmin"]-50000, bb["xmax"]+50000), ylim = c(bb["ymin"]-50000, bb["ymax"]+50000)) +
    theme_bw() +
    theme(panel.grid = element_blank())
)

  if(!dir.exists("Plots")) dir.create("Plots")



### Export plots of hexgrids w/ imputed tracks and original points
  PTTs <- unique(dat.sf.l$ptt)

  for (i in 1:length(PTTs)) {
    message("plotting ", PTTs[i], " ...")

    lines <- dat.sf.l %>%
      filter(ptt == PTTs[i])
    pts <- dat.sf.p %>%
      filter(ptt == PTTs[i])
    orig.pts <- dat.orig %>%
      filter(Ptt == PTTs[i])
    hexes <- hex_grid[[i]]$poly
    bb <- st_bbox(hexes)

    p_data <- ggplot() +
      geom_sf(data = gom.sf, size = 0.2) +
      geom_sf(data = hexes, fill = "transparent", size = 0.25) +
      geom_sf(data = lines, aes(color = rep), size = 0.25, alpha = 0.6) +
      scale_color_viridis_d(option = "turbo", guide = "none") +
      # geom_sf(data = pts, aes(color = rep), size = 0.25, alpha = 0.6) +
      geom_sf(data = orig.pts, size = 1) +
      coord_sf(xlim = c(bb["xmin"] - 10000, bb["xmax"] + 10000), ylim = c(bb["ymin"] - 10000, bb["ymax"] + 10000)) +
      ggtitle(paste("PTT", PTTs[i])) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            title = element_text(size = 18, face = "bold"))

    nm = paste0("Plots/p_", PTTs[i],"_data.pdf")
    ggsave(nm, p_data, dpi = "retina", width = 7, height = 7, units = 'in')

  }


## Load environmental layers
sst <- rast("Environ_data/GoM SST example.tif")
sst.proj <- sst %>%
  project('EPSG:3395') #%>%
  # raster::stack()  #MUST be of class RasterStack in order for path2ctmc() to work

# sst.proj.df <- as.data.frame(sst.proj, xy = TRUE)
# names(sst.proj.df)[3:6] <- c("Aug", "Sep", "Oct", "Nov")
# sst.proj.df <- sst.proj.df %>%
#   pivot_longer(cols = -c('x', 'y'), names_to = "month", values_to = "sst")
#
#
# ggplot() +
#   geom_raster(data = sst.proj.df %>% filter(month == "Nov"), aes(x, y, fill = sst)) +
#   scale_fill_cmocean("SST") +
#   geom_sf(data = hex_grid$`181800`$poly, fill = "transparent", size = 0.25) +
#   geom_sf(data = dat.sf.l %>%
#             filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
#   scale_color_manual(values = rainbow(20), guide = "none") +
#   ggtitle("Nov. SST") +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   # geom_sf(data = dat.sf.p %>%
#   #           filter(ptt == 181800), aes(color = rep), size = 0.25, alpha = 0.6) +
#   theme_bw() +
#   theme(panel.grid = element_blank())




### Export products ###

# save(dat, hex_grid, gom.sf, file = "Data_products/discretized_tracks.RData")


