
### addImageQuery ############################################################
##############################################################################
#' Add image query functionality to leaflet/mapview map.
#'
#' @details
#' This function enables Raster*/stars objects added to leaflet/mapview maps to
#' be queried. Standard query is on 'mousmove', but can be changed to 'click'.
#' Note that for this to work, the \code{layerId} needs to be the same as the
#' one that was set in \code{\link[leaflet]{addRasterImage}} or
#' \code{\link{addStarsImage}}. Currently only works for
#' numeric values (i.e. numeric/integer and factor values are supported).
#'
#' @param map the map with the RasterLayer to be queried.
#' @param x the RasterLayer that is to be queried.
#' @param band for stars layers, the band number to be queried.
#' @param group the group of the RasterLayer to be queried.
#' @param layerId the layerId of the RasterLayer to be queried. Needs to be the
#'   same as supplied in \code{\link[leaflet]{addRasterImage}} or
#'   \code{\link{addStarsImage}}.
#' @param project whether to project the RasterLayer to conform with leaflets
#'   expected crs. Defaults to \code{TRUE} and things are likely to go haywire
#'   if set to \code{FALSE}.
#' @param type whether query should occur on 'mousemove' or 'click'. Defaults
#'   to 'mousemove'.
#' @param digits the number of digits to be shown in the display field.
#' @param position where to place the display field. Default is 'topright'.
#' @param prefix a character string to be shown as prefix for the layerId.
#' @param className a character string to append to the control legend.
#' @param ... currently not used.
#'
#' @return
#' A leaflet map object.
#'
#' @examples
#' if (interactive()) {
#'   if (requireNamespace("plainview")) {
#'     library(leaflet)
#'     library(plainview)
#'
#'     leaflet() %>%
#'       addProviderTiles("OpenStreetMap") %>%
#'       addRasterImage(poppendorf[[1]], project = TRUE, group = "poppendorf",
#'                      layerId = "poppendorf") %>%
#'       addImageQuery(poppendorf[[1]], project = TRUE,
#'                     layerId = "poppendorf") %>%
#'       addLayersControl(overlayGroups = "poppendorf")
#'   }
#' }
#'
#' @importFrom raster projectExtent projectRaster as.matrix
#'
#' @export addImageQuery
#' @name addImageQuery
#' @rdname addImageQuery
addImageQuery = function(map,
                         x,
                         band = 1,
                         group = NULL,
                         layerId = NULL,
                         project = TRUE,
                         type = c("mousemove", "click"),
                         digits,
                         position = 'topright',
                         prefix = 'Layer',
                         className = "",
                         ...) {
  
  if (inherits(map, "mapview")) map = mapview2leaflet(map)
  
  type = match.arg(type)
  if (missing(digits)) digits = "null"
  if (is.null(group)) group = "stars"
  if (is.null(layerId)) layerId = group
  
  jsgroup <- gsub(".", "", make.names(group), fixed = TRUE)
  
  tmp <- leafem:::makepathStars(as.character(jsgroup))
  pathDatFn <- tmp[[2]][1]
  # starspathDatFn <- tmp[[3]][1]
  # datFn <- tmp[[4]][1]
  
  if (project) {
    if (inherits(x, "stars")) {
      if (utils::packageVersion("stars") >= "0.4-1") {
        projected = stars::st_warp(x, crs = 4326)
      } else {
        projected <- sf::st_transform(x, crs = 4326)
      }
    }
    if (inherits(x, "SpatRaster")) {
      projected <- terra::project(x, y = "EPSG:4326")
      }
    if (inherits(x, "Raster")) {
      projected = raster::projectRaster(
        x
        , raster::projectExtent(x, crs = sf::st_crs(4326)$proj4string)
        , method = "ngb"
      )
    }
  } else {
    projected <- x
  }
  
  pre <- paste0('var data = data || {}; data["', layerId, '"] = ')
  writeLines(pre, pathDatFn)
  cat('[', image2Array(projected, band = band), '];',
      file = pathDatFn, sep = "", append = TRUE)
  
  ## check for existing layerpicker control
  ctrlid = leafem:::getCallEntryFromMap(map, "addControl")
  ctrl_nm = paste("imageValues", layerId, sep = "-")
  imctrl = unlist(sapply(ctrlid, function(i) {
    ctrl_nm %in% map$x$calls[[i]]$args
  }))
  ctrlid = ctrlid[imctrl]
  
  # map = leaflet::clearControls(map)
  
  if (length(ctrlid) == 0) {
    # must add empty character instead of NULL for html with addControl
    map = leaflet::addControl(
      map,
      html = "",
      layerId = ctrl_nm,
      position = position,
      className = paste("info legend", className)
    )
  }
  
  sm <- leafem:::createFileId() #sample(1:1000, 1)
  map$dependencies <- c(map$dependencies,
                        leafem:::starsDataDependency(jFn = pathDatFn,
                                            counter = 1,
                                            group = paste0(layerId,"_",sm)))
  map$dependencies = c(map$dependencies,
                       list(htmltools::htmlDependency(
                         version = "0.0.1",
                         name = "joda",
                         src = system.file("htmlwidgets/lib/joda",
                                           package = "leafem"),
                         script = c("joda.js",
                                    "addImageQuery-bindings.js"))
                       ))
  
  bounds <- as.numeric(sf::st_bbox(projected))
  
  leaflet::invokeMethod(
    map
    , NULL
    , "addImageQuery"
    , layerId
    , bounds
    , type
    , digits
    , prefix
  )
}

###################################


stars2Array = function(x, band = 1) {
  if(length(dim(x)) == 2) layer = x[[1]] else layer = x[[1]][, , band]
  paste(
    sapply(seq(nrow(x[[1]])), function(i) {
      paste0(
        '['
        , gsub(
          "NA"
          , "null"
          , paste(as.numeric(layer[i, ]), collapse = ",")
        )
        , ']'
      )
    }),
    collapse = ","
  )
}


###################################


SpatRaster2Array = function(x) {
  x = terra::as.matrix(x, wide = TRUE)
  paste(
    sapply(seq(ncol(x)), function(i) {
      paste0(
        '['
        , gsub(
          "NA"
          , "null"
          , paste(as.matrix(x)[, i], collapse = ",")
        )
        , ']'
      )
    }),
    collapse = ","
  )
}


###################################


rasterLayer2Array = function(x) {
  x = as.matrix(x)
  paste(
    sapply(seq(ncol(x)), function(i) {
      paste0(
        '['
        , gsub(
          "NA"
          , "null"
          , paste(as.matrix(x)[, i], collapse = ",")
        )
        , ']'
      )
    }),
    collapse = ","
  )
}


###################################


image2Array = function(x, band = 1) {
  switch(class(x)[1],
         "SpatRaster" = SpatRaster2Array(x),
         "stars" = stars2Array(x, band = band),
         "RasterLayer" = rasterLayer2Array(x))
}


###################################


# Function that modifies existing leaflet::addLegend by adding an option for decreasing order
addLegend_decreasing <- function (map,
                                  position = c("topright", "bottomright", "bottomleft",
                                               "topleft"),
                                  pal,
                                  values,
                                  na.label = "NA",
                                  bins = 7,
                                  colors,
                                  opacity = 0.5,
                                  labels = NULL,
                                  labFormat = labelFormat(),
                                  title = NULL, className = "info legend", layerId = NULL,
                                  group = NULL, data = getMapData(map), decreasing = FALSE) {
  position <- match.arg(position)
  type <- "unknown"
  na.color <- NULL
  extra <- NULL
  if (!missing(pal)) {
    if (!missing(colors))
      stop("You must provide either 'pal' or 'colors' (not both)")
    if (missing(title) && inherits(values, "formula"))
      title <- deparse(values[[2]])
    values <- evalFormula(values, data)
    type <- attr(pal, "colorType", exact = TRUE)
    args <- attr(pal, "colorArgs", exact = TRUE)
    na.color <- args$na.color
    if (!is.null(na.color) && col2rgb(na.color, alpha = TRUE)[[4]] ==
        0) {
      na.color <- NULL
    }
    if (type != "numeric" && !missing(bins))
      warning("'bins' is ignored because the palette type is not numeric")
    if (type == "numeric") {
      cuts <- if (length(bins) == 1)
        pretty(values, bins)
      else bins
      
      if (length(bins) > 2)
        if (!all(abs(diff(bins, differences = 2)) <=
                 sqrt(.Machine$double.eps)))
          stop("The vector of breaks 'bins' must be equally spaced")
      n <- length(cuts)
      r <- range(values, na.rm = TRUE)
      cuts <- cuts[cuts >= r[1] & cuts <= r[2]]
      n <- length(cuts)
      p <- (cuts - r[1])/(r[2] - r[1])
      extra <- list(p_1 = p[1], p_n = p[n])
      p <- c("", paste0(100 * p, "%"), "")
      if (decreasing == TRUE){
        colors <- pal(rev(c(r[1], cuts, r[2])))
        labels <- rev(labFormat(type = "numeric", cuts))
      }else{
        colors <- pal(c(r[1], cuts, r[2]))
        labels <- rev(labFormat(type = "numeric", cuts))
      }
      colors <- paste(colors, p, sep = " ", collapse = ", ")
      
    }
    else if (type == "bin") {
      cuts <- args$bins
      n <- length(cuts)
      mids <- (cuts[-1] + cuts[-n])/2
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "bin", cuts))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "bin", cuts)
      }
      
    }
    else if (type == "quantile") {
      p <- args$probs
      n <- length(p)
      cuts <- quantile(values, probs = p, na.rm = TRUE)
      mids <- quantile(values, probs = (p[-1] + p[-n])/2,
                       na.rm = TRUE)
      if (decreasing == TRUE){
        colors <- pal(rev(mids))
        labels <- rev(labFormat(type = "quantile", cuts, p))
      }else{
        colors <- pal(mids)
        labels <- labFormat(type = "quantile", cuts, p)
      }
    }
    else if (type == "factor") {
      v <- sort(unique(na.omit(values)))
      colors <- pal(v)
      labels <- labFormat(type = "factor", v)
      if (decreasing == TRUE){
        colors <- pal(rev(v))
        labels <- rev(labFormat(type = "factor", v))
      }else{
        colors <- pal(v)
        labels <- labFormat(type = "factor", v)
      }
    }
    else stop("Palette function not supported")
    if (!any(is.na(values)))
      na.color <- NULL
  }
  else {
    if (length(colors) != length(labels))
      stop("'colors' and 'labels' must be of the same length")
  }
  legend <- list(colors = I(unname(colors)), labels = I(unname(labels)),
                 na_color = na.color, na_label = na.label, opacity = opacity,
                 position = position, type = type, title = title, extra = extra,
                 layerId = layerId, className = className, group = group)
  invokeMethod(map, data, "addLegend", legend)
}


###################################


# Function to change CRS of coords, but leave as data.frame instead of converting to {sf} obj
st_changeCoords = function(data, coords, crs_in, crs_out) {
  # returns data.frame w/ cols named 'x' and 'y' for transformed coords
  
  dat.sf <- sf::st_as_sf(data, coords = coords, crs = crs_in) %>% 
    sf::st_transform(dat.sf, crs = crs_out) %>% 
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) %>% 
    st_drop_geometry()
  
  return(dat.sf)
}
