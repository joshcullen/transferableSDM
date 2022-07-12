
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

  # extr.covar<- data.frame()
  extr.covar<- matrix(NA, nrow = nrow(tmp) - 1, ncol = 6 + length(layers))
  colnames(extr.covar) <- c("n", "dist", "dt", "id", "date", "state", names(layers))

  #Extract values from each line segment
  for (j in 2:nrow(tmp)) {
    # print(j)
    segment<- tmp[(j-1):j, c("x","y")] %>%
      as.matrix() %>%
      st_linestring() %>%
      st_sfc(crs = 'EPSG:3395') %>%
      st_transform(terra::crs(layers[[1]]))


    # Create time-matched raster stack
    time.ind <- map(layers[dyn_names], ~which(names(.x) == tmp[[ind]][j]))
    layers.tmp <- layers

    # Replace missing raster for time interval w/ NA-filled raster
    if (length(time.ind) != length(unlist(time.ind))) {
      cond <- which(map(time.ind, length) == 0)
      time.ind[cond] <- 1
      layers.tmp[[names(cond)]] <- layers.tmp[[names(cond)]][[1]]
      terra::values(layers.tmp[[names(cond)]]) <- NA
    } else {
      dyn_names1 <- dyn_names
    }
    layers.tmp[dyn_names] <- map2(layers.tmp[dyn_names], time.ind, ~{.x[[.y]]})
    layers.tmp <- rast(layers.tmp)


    # .layers.tmp <- terra::wrap(layers.tmp)
    # .segment <- terra::wrap(terra::vect(segment))
    tmp1<- terra::extract(layers.tmp, terra::vect(segment), along = TRUE,
                          cellnumbers = FALSE) #%>%
      # purrr::map(., ~matrix(., ncol = nlyr(layers.tmp))) %>%
      # purrr::pluck(1) %>%
      # data.frame() %>%
      # purrr::set_names(names(layers))

    #subset to only include time-matched vars (by some indicator variable)
    # if (!is.null(ind)) {
    #
    #   cond<- tmp[j-1, ind]
    #   cond2<- levels(cond)[which(cond != levels(cond))]
    #   tmp1<- tmp1[,!stringr::str_detect(names(tmp1), paste(cond2, collapse="|")), drop=F]
    #
    #   ind1<- stringr::str_which(names(tmp1), as.character(cond))
    #   names(tmp1)[ind1]<- dyn_names
    #
    # }


    #calculate segment means if continuous and proportions spent in each class if categorical
    if (is.null(which_cat)) {
      covar.means<- data.frame(t(colMeans(tmp1)))
    } else {
      covar.means<- data.frame(t(colMeans(tmp1[,names(tmp1) != which_cat, drop = FALSE])))
      cat<- factor(tmp1[,which_cat], levels = lev)
      cat<- data.frame(unclass(t(table(cat)/length(cat))))

      covar.means<- cbind(covar.means, cat)  #Merge continuous and categorical vars
    }


    tmp2<- cbind(n = nrow(tmp1), dist = tmp$step[j-1], covar.means) %>%
      dplyr::mutate(dt = as.numeric(tmp$dt[j-1]),
                    id = unique(data$id),
                    date = tmp$date[j-1],
                    state = ifelse(!is.null(state.col), tmp[j-1,state.col], NA),
                    .after = dist) %>%
      dplyr::select(-ID)
    # tmp2<- tmp2[,!apply(is.na(tmp2), 2, any)]

    # extr.covar<- rbind(extr.covar, tmp2)
    extr.covar[j-1,] <- unlist(tmp2)
    # rm(.segment)
    # rm(.layers.tmp)
  }

  extr.covar <- as.data.frame(extr.covar) %>%
    mutate(date = as_datetime(date))

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


