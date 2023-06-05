
# Functions to hexify rasters

hexify_df = function(bounds, raster_list, hex_list, quad_list){
  # bounds = a vector of start and end dates in POSIXct format for a single PTT
  # raster_list = a list containing a {terra} SpatRaster object; can include 1 or more layers that represent different times
  # hex_list = a list object that contains an element named 'poly', which holds an {sf} POLYGON object representing the hexgrid for a single PTT
  # quad_list = a list of quadrature (time) points for the layers of each raster


  #create empty list to store results
  hex_cov1 <- vector("list", length(raster_list))

  #hexify each covariate
  for (j in 1:length(raster_list)) {
    #create vector of dates for indexing raster layers
    if(length(quad_list) == 1) {
      grid <- as_date(quad_list[[1]])
    } else {
      grid <- as_date(quad_list[[j]])
    }

    #create index to subset raster based on tracking period
    if(length(grid) == 1) {
      idx <- 1  #for static raster layers
    } else {
      bounds.seq <- seq(bounds[1], bounds[2], by = "month")
      idx <- which(grid %in% bounds.seq)
    }

    hex_cov1[[j]] <- hexify(x = raster_list[[j]], hex_list = hex_list, idx = idx)

    message(names(raster_list)[j], " hexified ...")
  }

  # wrangle hexified data into data.frame for each covariate
  hex_cov_df1 <- hex_cov1[[1]]
  hex_cov_df1$Time <- bounds[1]  ##ONLY DO THIS IF BATHYMETRY OR OTHER STATIC VARIABLE IS FIRST; sets all 'Time' values to that of starting month for given PTT
  if(length(raster_list) > 1) {
    for (j in 2:length(hex_cov1)) {
      hex_cov_df1 <- hex_cov_df1 %>%
        full_join(hex_cov1[[j]], by = c("hex", "Time"))
    }
  }
  names(hex_cov_df1)[-c(1:2)] <- names(raster_list)
  hex_cov_df1 <- hex_cov_df1 %>%
    arrange(hex, Time)
  # hex_cov_df1 <- lapply(hex_cov_df1, zoo::na.locf, na.rm = FALSE) %>%
  #   as.data.frame()

  # add values for every quadrature point for static covariates
  static <- which(map(raster_list, terra::nlyr) == 1)  #index static covariates
  if (length(static) > 0) {
    hex_cov_df1 <- hex_cov_df1 %>%
      mutate(across(names(static), zoo::na.locf))  #ensure that static values are consistent (and not NA) across all time intervals
  }


  return(hex_cov_df1)
}





hexify = function(x, hex_list, idx){
  # x = a single {terra} SpatRaster object representing a single covariate; can include 1 or more layers that represent different times
  # hex_list = a list object that contains an element named 'poly', which holds an {sf} POLYGON object representing the hexgrid for a single PTT
  # idx =  a vector of one or more integers used to subset the SpatRaster object 'x' based on the quadrature points


  #subset raster layer by index from tracking period
  x <- x[[idx]]

  #define projections for hexgrid and raster layer
  proj_h <- st_crs(hex_list$poly)
  if(is.na(proj_h)) stop("hex_list poly element must have CRS")
  proj_x <- crs(x)
  if(str_count(proj_x, "\\+proj=longlat")==1) {
    if(extent(x)@xmax > 180){
      x <- rotate(x)  #changes longitude to -180, 180 from 0, 360 range
    }
  }

  # Change projection of hexgrid and extract weighted values for each cell
  poly <- st_transform(hex_list$poly, proj_x)
  poly.ext <- ext(vect(poly))
  x <- crop(x, extend(poly.ext, res(x)))
  out <- terra::extract(x, vect(poly), weights = TRUE)


  # Calculate weighted mean value per hex grid cell
  out1 <- out %>%
    split(.$ID) %>%
    map(wt.mean) %>%
    bind_rows(.id = "hex")

  # in case there's an issue where a 'hex' column isn't created (potentially due to many NAs)
  if (!any(names(out1) == 'hex') & nrow(out1) == 1) {
    out1 <- data.frame(hex = unique(out$ID),
                       t(out1)
    )
    names(out1)[2] <- names(out)[2]
  }

  # Wrangle data to long-format
  out2 <- out1 %>%
    pivot_longer(cols = -hex, names_to = 'Time', values_to = 'val') %>%
    mutate(across(hex, as.numeric)) %>%
    arrange(hex, Time) %>%
    mutate(Time = as_date(Time))

  return(out2)
}




## Weighted mean wrapper function
  # m = a data.frame where rows represent each square raster cell that overlaps a given hexgrid cell and columns represent different quadrature (time) points; the last column is the weight attributed to each of the square raster cells


wt.mean <- function(m){
  if(is.null(m)) stop("error")

  if(nrow(m) > 1) {
    wt.vals <- apply(m[,-c(1, ncol(m)), drop=FALSE], 2, weighted.mean, w = m[,ncol(m)])
  } else {
    wt.vals <- m[,-c(1, ncol(m))]
  }

  return(wt.vals)
}
