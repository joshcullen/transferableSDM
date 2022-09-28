
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

# function to convert list object from rerddapXtracto::rxtracto() into a {terra} SpatRaster object
array2rast <- function(lon, lat, var, time, extent) {
  #lon: a vector of longitude values
  #lat: a vector of latitude values
  #var: an array of values for the variable of interest where the 1st two dimensions denote the spatial grid and the third dimensions represents the number of datetimes
  #time: a vector of dates associated with each 2D array from `var`
  #extent: a SpatRaster extent on which to spatially define the raster; ordered as xmin, xmax, ymin, ymax

  dims <- dim(var)

  # check dims in case "altitude/depth" is included; this code will only select 1st slice of altitude
  if (length(dims) == 4) {
    var <- var[,,1,]
  } else {
    var <- var
  }

  dims1 <- dim(var)

  rast.list <- vector("list", dims1[3])

  for (i in 1:dims1[3]) {
    rast.list[[i]] <- terra::rast(t(var[,,i]), crs = 'EPSG:4326', extent = extent) %>%
      terra::flip(direction = "vertical")
  }

  rast1 <- terra::rast(rast.list)
  names(rast1) <- time

  return(rast1)
}

#---------------------------
# Create a LINESTRING per row given start and end points
#must be labeled x1, y1, x2, y2

make_line <- function(x1, y1, x2, y2) {
  st_linestring(matrix(c(x1, x2, y1, y2), 2, 2))
}

#---------------------------

extract.covars.internal = function(data, layers, state.col, which_cat, dyn_names, ind, imputed, p) {
  ## data = data frame containing at least the id, coordinates (x,y), date-time (date), and
  ##      step length (step); if imputations are also included per id, then the column labeling the separate imputations (rep) should also be included
  ## layers = a raster object (Raster, RasterStack, RasterBrick) object containing environ covars
  ## state.col = character. The name of the column that contains behavioral states w/in
  ##             data (if present)
  ## which_cat = vector of names or numeric positions of discrete raster layers; NULL by default
  ## dyn_names = vector of names dynamic raster layers (in same order as layers); NULL by default
  ## ind = character/integer. The name or column position of the indicator column of data to be
  ##       matched w/ names of a dynamic raster
  ## imputed = logical. If TRUE, the functions calculate time intervals while also accounting for separate imputations
  ## p = a stored 'progressr' object for creating progress bar



  #Subset and prep data
  if (is.null(data$dt) & imputed == FALSE) {

    tmp<- data %>%
      # dplyr::filter(id == ind[i]) %>%
      dplyr::mutate(dt = difftime(date, dplyr::lag(date, 1), units = "secs")) %>%
      dplyr::mutate_at("dt", {. %>%
          as.numeric() %>%
          round()})
    tmp$dt<- c(purrr::discard(tmp$dt, is.na), NA)

  } else if (is.null(data$dt) & imputed == TRUE) {

    tmp<- data %>%
      split(.$rep) %>%
      map(., {. %>%
          dplyr::mutate(dt = difftime(date, dplyr::lag(date, 1), units = "secs")) %>%
          dplyr::mutate(dt = dt %>%
                          as.numeric() %>%
                          round())
      }) %>%
      bind_rows()
    tmp$dt<- c(purrr::discard(tmp$dt, is.na), NA)

  } else {
    tmp <- data
  }


  # if (!is.null(dyn_names) & !is.factor(tmp[,ind])) stop("The `ind` column must be a factor.")



  #Identify levels of categorical layer (if available)
  if (!is.null(which_cat)) {
    lev<- layers[[which_cat]]@data@attributes[[1]][,1]
  }

  #Turn dataset into rows of individual LINESTRINGs (in EPSG:3395)
  segment <- tmp %>%
    select(x1, y1, x2, y2) %>%
    pmap(make_line) %>%
    st_as_sfc(crs = 3395) %>%
    {tibble(month.year = tmp$month.year,
            strata = tmp$strata,
            obs = tmp$obs,
            geometry = .)} %>%
    st_sf()

  #Extract values from each line segment
  for (j in 1:n_distinct(segment$month.year)) {
    # print(j)

    #create subsetted data.frame of original given selected month.year
    tmp.sub <- tmp[tmp[[ind]] == unique(segment[[ind]])[j],]


    # Create time-matched raster stack
    time.ind <- map(layers[dyn_names], ~which(names(.x) == unique(segment[[ind]])[j]))
    layers.tmp <- layers

    # Replace missing raster for time interval w/ NA-filled raster
    if (length(time.ind) != length(unlist(time.ind))) {
      cond <- which(map(time.ind, length) == 0)
      time.ind[cond] <- 1
      layers.tmp[[names(cond)]] <- layers.tmp[[names(cond)]][[1]]
      terra::values(layers.tmp[[names(cond)]]) <- NA
    }

    layers.tmp[dyn_names] <- map2(layers.tmp[dyn_names], time.ind, ~{.x[[.y]]})
    layers.tmp <- rast(layers.tmp)


    tmp1<- terra::extract(layers.tmp, terra::vect(segment[segment[[ind]] == unique(segment[[ind]])[j],]),
                          along = TRUE, cells = FALSE)




    # aggregate data per step and calculate mean values
    tmp2 <- tmp1 %>%
      group_by(ID) %>%
      summarize(across(names(layers.tmp), mean, na.rm = TRUE)) %>%
      dplyr::mutate(n = as.vector(table(tmp1$ID)),
                    dist = tmp.sub$step,
                    .before = everything()) %>%
      dplyr::mutate(dt = as.numeric(tmp.sub$dt),
                    id = unique(data$id),
                    date = tmp.sub$date,
                    state = ifelse(!is.null(state.col), tmp.sub[,state.col], NA),
                    .after = dist) %>%
      dplyr::select(-ID)


    if (j == 1) {
      extr.covar <- tmp2
    } else {
      extr.covar<- rbind(extr.covar, tmp2)
    }

  }

  # extr.covar <- extr.covar %>%
  #   mutate(date = as_datetime(date))

  p()  #plot progress bar
  extr.covar
}

#----------------------------
extract.covars = function(data, layers, state.col = NULL, which_cat = NULL, dyn_names = NULL,
                          ind, imputed = TRUE) {
  ## data must be a data frame with "id" column, coords labeled "x" and "y" and datetime as POSIXct labeled "date"; optionally can have column that specifies behavioral state; if imputed = TRUE, a column named "rep" must be included to distinguish among imputations

  message("Prepping data for extraction...")
  tictoc::tic()

  dat.list <- bayesmove::df_to_list(data, "id")
  dat.list.rep <- map(dat.list, ~bayesmove::df_to_list(., ind = "rep"))

  ## Make raster data (stored in `layers`) usable in parallel
  .layers <- map(layers, terra::wrap)

  ## Create empty list to store results
  path.list <- vector("list", length = length(dat.list)) %>%
    purrr::set_names(names(dat.list))


  ## Analyze across IDs using for-loop
  for (i in 1:length(dat.list.rep)) {

    message("Extracting environmental values for PTT ", names(dat.list.rep)[i], "...")

    progressr::with_progress({
      #set up progress bar
      p<- progressr::progressor(steps = length(dat.list.rep[[i]]))

      # tictoc::tic()
      path <- furrr::future_map(dat.list.rep[[i]],
                                ~extract.covars.internal(data = .x, layers = map(.layers, terra::rast),
                                                         state.col = state.col,
                                                         which_cat = which_cat,
                                                         dyn_names = dyn_names, ind = ind,
                                                         imputed = imputed, p = p),
                                .options = furrr_options(seed = TRUE))
      # tictoc::toc()
    })

    path <- dplyr::bind_rows(path, .id = "rep")
    path.list[[i]] <- path
  }

  message("Exporting extracted values...")

  path.out <- dplyr::bind_rows(path.list)
  tictoc::toc()


  return(path.out)
}


#----------------------------

add_avail_steps <- function(data) {
  ## data = data.frame containing at least columns `id`, `x`, `y`, and `step`
  ## this function creates a set of 4 available steps per observed step in the 4 cardinal directions
  ## observed and available steps are grouped together by added column `strata`

  data2 <- data %>%
    rename(x1 = x, y1 = y) %>%
    split(.$rep) %>%
    map(., ~{.x %>%
        mutate(x2 = c(x1[-1], NA),
               y2 = c(y1[-1], NA),
               .after = y1)
      }) %>%
    bind_rows()

  data2$strata <- 1:nrow(data2)
  data2$obs <- 1


  ## define available steps
  avail.list <- vector("list", 4)
  tmp <- data2

  for (i in 1:4) {

    if (i == 1) {
      #step north
      tmp$x2 <- tmp$x1
      tmp$y2 <- tmp$y1 + tmp$step
    } else if (i == 2) {
      #step east
      tmp$x2 <- tmp$x1 + tmp$step
      tmp$y2 <- tmp$y1
    } else if (i == 3) {
      #step south
      tmp$x2 <- tmp$x1
      tmp$y2 <- tmp$y1 - tmp$step
    } else {
      #step west
      tmp$x2 <- tmp$x1 - tmp$step
      tmp$y2 <- tmp$y1
    }

    tmp$obs <- 0
    avail.list[[i]] <- tmp
  }

  avail <- bind_rows(avail.list)

  #add available steps to observed steps
  data3 <- rbind(data2, avail) %>%
    arrange(strata)

  return(data3)

}
