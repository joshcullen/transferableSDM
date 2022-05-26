hex_grid_sfsample <- function(sf_tracks, cellsize, buffer) {
  attr(st_geometry(sf_tracks), "bbox") <- st_bbox(sf_tracks) + c(rep(-buffer,2), rep(buffer,2))
  crs1 <- st_crs(sf_tracks)

  hex_poly <- st_make_grid(sf_tracks, cellsize = cellsize, offset = c(0.5,0.5),
                           square = FALSE)
  hex_cent <- st_centroid(hex_poly)
  hex_nb = dnearneigh(hex_cent, cellsize-0.01*cellsize, cellsize+0.01*cellsize)

  #Pick back up from here
  visit = unique(over(spdata, hex_poly))
  neighbor_df = data.frame(
    hex_to = unlist(hex_nb),
    hex_from = rep(1:length(hex_cent), sapply(hex_nb, length))
  ) %>% filter(hex_from %in% visit)

  idx = sort(unique(c(neighbor_df$hex_from, neighbor_df$hex_to)))

  hex_poly = hex_poly[idx]
  row.names(hex_poly) = as.character(1:length(hex_poly))

  hex_cent = hex_cent[idx]
  hex_nb = dnearneigh(hex_cent, cellsize-0.01*cellsize, cellsize+0.01*cellsize)

  visit = unique(over(spdata, hex_poly))
  neighbor_df = data.frame(
    hex_from = rep(1:length(hex_cent), sapply(hex_nb, length)),
    hex_to = unlist(hex_nb)
  )
  neighbor_df$visit = neighbor_df$hex_from %in% visit

  to_next =  coordinates(hex_cent[neighbor_df$hex_to,]) - coordinates(hex_cent[neighbor_df$hex_from,])
  neighbor_df$bearing_to_next = (atan2(to_next[,1], to_next[,2])*180/pi) %>% ifelse(.<0, 360+., .)
  unit_xy = cbind(cos(pi*(90-neighbor_df$bearing_to_next)/180), sin(pi*(90-neighbor_df$bearing_to_next)/180))
  neighbor_df$north = as.vector(round(unit_xy %*% c(0,1),3))
  neighbor_df$east = as.vector(round(unit_xy %*% c(1,0),3))

  sp::proj4string(hex_poly) = p4s
  sp::proj4string(hex_cent) = p4s

  output <- list(
    poly=hex_poly,
    neighbor_df = neighbor_df,
    hex_centroids = hex_cent
  )
}
