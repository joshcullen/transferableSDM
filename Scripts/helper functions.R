
# function to convert list object from rerddapXtracto::rxtracto() into a data.frame for ggplot
array2df <- function(lon, lat, var, var.name, time) {
  dims <- dim(var)
  rast.df <- expand.grid(x = lon, y = lat, time = time)
  rast.df$var <- array(apply(var, 3, rbind), dims[1] * dims[2] * dims[3])
  names(rast.df)[4] <- var.name
  rast.df <- data.frame(rast.df) %>%
    mutate(across(everything(), as.vector))

  return(rast.df)
}

#---------------------------

### Utility functions to download GEBCO_2019 bathymetry data from Github gist
# https://gist.github.com/mdsumner/dea21c65674108574bab22cf6f011f8d

.gebco_template <- function() {
  #ri <- vapour::vapour_raster_info(gebco)
  #ri <- list(extent = ri$extent, dimension = ri$dimXY, projection  = ri$projection)
  #dput(ri)
  list(extent = c(-180, 180, -90, 90), dimension = c(86400L, 43200L), projection = "OGC:CRS84")
}

#---------------------------

.gebco_template_terra <- function() {
  ri <- .gebco_template()
  terra::rast(terra::ext(ri$extent), ncols = ri$dimension[1L], nrows = ri$dimension[2L], crs = ri$projection)
}

#---------------------------

# input lon,lat (and optionally width,height in metres) - default for an equal width,height region (118km default)
##    gives 256x<256+more depending on latitude>
# OR set y to a SpatRaster (for a small region, but in any projection)
#' @examples
#' get_gebco(c(147, -42))
#' get_gebco(terra::rast())
get_gebco <- function(x, y = rep(1.18e5, 2L), ..., resample = c("near", "bilinear", "cubic", "cubicspline" ), snap = "out", debug = FALSE) {
  if (missing(x)) x<- geosphere::randomCoordinates(1L)
  gebco <- "/vsicurl/https://public.services.aad.gov.au/datasets/science/GEBCO_2019_GEOTIFF/GEBCO_2019.tif"
  template <- .gebco_template_terra()
  if (is.numeric(x) && length(x) == 2L) { ##we think x is a longlat point
    stopifnot(is.numeric(y) && length(y) == 2L)
    coslat <- cos(pi/180 * x[2L])
    ex <-    rep(x, each = 2L) + (c(-1/coslat, 1/coslat, -1, 1) *  rep(y, each = 2L)/2) / 111111.12

    x <- terra::crop(template, terra::ext(ex), snap = snap)

  } else {
    if (inherits(x, "SpatRaster")) {

    } else {
      stop("can't interpret 'x', input x = c(lon, lat), y = c(width, height) OR input x = \"SpatRaster\"")
    }
  }
  extent <- c(terra::xmin(x), terra::xmax(x), terra::ymin(x), terra::ymax(x))
  dimension <- dim(x)[2:1]
  projection <- terra::crs(x)

  if (debug) {
    return(x)
  }
  v <- vapour::vapour_warp_raster(gebco, bands = 1L, extent = extent, dimension = dimension,
                                  projection = projection, resample = resample)

  terra::setValues(x, v[[1L]])
}
