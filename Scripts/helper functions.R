
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

extract.covars.internal = function(data, layers, state.col, which_cat, dyn_names, ind, along, imputed, p) {
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
  # if (is.null(data$dt) & imputed == FALSE) {
  #
  #   tmp<- data %>%
  #     # dplyr::filter(id == ind[i]) %>%
  #     dplyr::mutate(dt = difftime(date, dplyr::lag(date, 1), units = "secs")) %>%
  #     dplyr::mutate_at("dt", {. %>%
  #         as.numeric() %>%
  #         round()})
  #   tmp$dt<- c(purrr::discard(tmp$dt, is.na), NA)
  #
  # } else if (is.null(data$dt) & imputed == TRUE) {
  #
  #   tmp<- data %>%
  #     split(.$rep) %>%
  #     map(., {. %>%
  #         dplyr::mutate(dt = difftime(date, dplyr::lag(date, 1), units = "secs")) %>%
  #         dplyr::mutate(dt = dt %>%
  #                         as.numeric() %>%
  #                         round())
  #     }) %>%
  #     bind_rows()
  #   tmp$dt<- c(purrr::discard(tmp$dt, is.na), NA)
  #
  # } else {
  #   tmp <- data
  # }

  tmp <- data

  # if (!is.null(dyn_names) & !is.factor(tmp[,ind])) stop("The `ind` column must be a factor.")



  #Identify levels of categorical layer (if available)
  if (!is.null(which_cat)) {
    lev<- layers[[which_cat]]@data@attributes[[1]][,1]
  }

  #Turn dataset into rows of individual LINESTRINGs (in EPSG:3395)

  if (along) {

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
                          along = along, cells = FALSE)




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
  } else if (!along) {

    # Extract environ covariate values by points (in EPSG:3395 proj)
    pts <- tmp %>%
      # select(x1, y1, x2, y2) %>%
      # pmap(make_line) %>%
      # st_as_sfc(crs = 3395) %>%
      # {tibble(month.year = tmp$month.year,
      #         strata = tmp$strata,
      #         obs = tmp$obs,
      #         geometry = .)} %>%
      # st_sf()
      sf::st_as_sf(., coords = c('x2','y2'), crs = 3395)

    #Extract values from each line segment
    for (j in 1:n_distinct(pts[[ind]])) {
      # print(j)

      #create subsetted data.frame of original given selected month.year
      tmp.sub <- tmp[tmp[[ind]] == unique(pts[[ind]])[j],]


      # Create time-matched raster stack
      time.ind <- map(layers[dyn_names], ~which(names(.x) == unique(pts[[ind]])[j]))
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


      tmp1<- terra::extract(layers.tmp, terra::vect(pts[pts[[ind]] == unique(pts[[ind]])[j],]),
                            along = along, cells = FALSE)




      # add extracted covariates to original dataset
      tmp2 <- tmp1 %>%
        # group_by(ID) %>%
        # summarize(across(names(layers.tmp), mean, na.rm = TRUE)) %>%
        # dplyr::mutate(n = as.vector(table(tmp1$ID)),
        #               dist = tmp.sub$step,
        #               .before = everything()) %>%
        # dplyr::mutate(dt = as.numeric(tmp.sub$dt),
        #               id = unique(data$id),
        #               date = tmp.sub$date,
        #               state = ifelse(!is.null(state.col), tmp.sub[,state.col], NA),
        #               .after = dist) %>%
        dplyr::select(-ID) %>%
        cbind(tmp.sub, .)



      if (j == 1) {
        extr.covar <- tmp2
      } else {
        extr.covar <- rbind(extr.covar, tmp2)
      }
    }
  }

  # extr.covar <- extr.covar %>%
  #   mutate(date = as_datetime(date))

  p()  #plot progress bar
  extr.covar
}

#----------------------------
extract.covars = function(data, layers, state.col = NULL, which_cat = NULL, dyn_names = NULL, along = TRUE,
                          ind, imputed = FALSE) {
  ## data must be a data frame with "id" column, coords labeled "x" and "y" and datetime as POSIXct labeled "date"; optionally can have column that specifies behavioral state; if imputed = TRUE, a column named "rep" must be included to distinguish among imputations

  message("Prepping data for extraction...")
  tictoc::tic()

  dat.list <- bayesmove::df_to_list(data, "id")
  # dat.list.rep <- map(dat.list, ~bayesmove::df_to_list(., ind = "rep"))

  ## Make raster data (stored in `layers`) usable in parallel
  .layers <- map(layers, terra::wrap)

  ## Create empty list to store results
  # path.list <- vector("list", length = length(dat.list)) %>%
  #   purrr::set_names(names(dat.list))


  ## Analyze across IDs using for-loop
  # for (i in 1:length(dat.list)) {

    message("Extracting environmental values for PTTs...")

    progressr::with_progress({
      #set up progress bar
      p<- progressr::progressor(steps = length(dat.list))

      # tictoc::tic()
      path <- furrr::future_map(dat.list,
                                ~extract.covars.internal(data = .x, layers = map(.layers, terra::rast),
                                                         state.col = state.col,
                                                         which_cat = which_cat,
                                                         dyn_names = dyn_names, ind = ind, along = along,
                                                         imputed = imputed, p = p),
                                .options = furrr_options(seed = TRUE))
      # tictoc::toc()
    })

    path <- dplyr::bind_rows(path, .id = "id")
    # path.list[[i]] <- path
  # }

  # message("Exporting extracted values...")

  # path.out <- dplyr::bind_rows(path.list)
  tictoc::toc()


  return(path)
}


#----------------------------

# Function to generate new coordinates for different forms of SSF


avail_steps <- function(data, nsteps) {

  # create list to store results for each strata
  avail.list <- vector("list", nrow(data))


  for (i in 1:nrow(data)) {

    ## generate random steps and angles
    rand.sl <- sample(data$step, size = nsteps)
    rand.bearing <- runif(nsteps, min = 0, max = 2*pi)


    if (i == 1) {  #add NA for turning angles during first step since not possible to calculate

      tmp <- data %>%
        slice(rep(i, each = nsteps)) %>%
        mutate(x2 = x1 + (rand.sl * cos(rand.bearing)),
               y2 = y1 + (rand.sl * sin(rand.bearing)),
               step = rand.sl,
               angle = NA,
               obs = 0
        )

    } else {

      coords_tmin1 <- data[i-1, c('x1','y1')]  #coordinates from previous used step

      tmp <- data %>%
        slice(rep(i, each = nsteps)) %>%
        mutate(x2 = x1 + (rand.sl * cos(rand.bearing)),
               y2 = y1 + (rand.sl * sin(rand.bearing)),
               step = rand.sl,
               obs = 0
        )

      # calculate turning angles using 3 sets of consecutive coordinates
      ang <- atan2(tmp$y1 - coords_tmin1$y1, tmp$x1 - coords_tmin1$x1) -
        atan2(tmp$y2 - tmp$y1, tmp$x2 - tmp$x1)
      ang <- ifelse(ang > pi, ang - 2 * pi,
                    ifelse(ang < -pi, 2 * pi + ang,
                           ang))

      tmp$angle <- ang
    }

    avail.list[[i]] <- tmp
  }

  avail <- bind_rows(avail.list)

  return(avail)
}


#----------------------------

add_avail_steps <- function(data, nsteps) {
  ## data = data.frame containing at least columns `id`, `x`, `y`, and `step`
  ## nsteps = integer. The number of random steps to characterize the available habitat. Step lengths drawn from empirical distrib and turning angles from uniform distrib
  ## observed and available steps are grouped together by added column `strata`

  tictoc::tic()

  data2 <- data %>%
    rename(x1 = x, y1 = y) %>%
    split(.$id) %>%
    map(., ~{.x %>%
        mutate(x2 = c(x1[-1], NA),
               y2 = c(y1[-1], NA),
               .after = y1)
      }) %>%
    bind_rows()

  data2$strata <- 1:nrow(data2)
  data2$obs <- 1


  ## generate coordinates for available steps; also adds new values for steps and angles
  dat.list <- split(data2, data2$id)
  avail <- dat.list %>%
    furrr::future_map(~avail_steps(data = ., nsteps = nsteps),
                      .options = furrr_options(seed = TRUE)) %>%
    bind_rows()



  #add available steps to observed steps
  data3 <- rbind(data2, avail) %>%
    arrange(strata)

  tictoc::toc()

  return(data3)

}

#----------------------------

# Center and scale across multiple layers (representing different time periods) for a given variable
scale_across_covar <- function(layer) {

  mean1 <- layer %>%
    terra::values() %>%
    mean(na.rm = TRUE)

  sd1 <- layer %>%
    terra::values() %>%
    sd(na.rm = TRUE)


  (layer - mean1) / sd1
}

#----------------------------

### Calculate probabilities of used and available steps based on time model results ###
calc_time_probs <- function(mod, dat, covar.names, p) {
  # mod = the rstan model fit for the time model
  # dat = the input for the time model that includes cols for id1, step, dt, strata, covars, and obs
  # covar names = a vector of names (in proper order) for the covars on which to make predictions
  # p = a stored 'progressr' object for creating progress bar

  # message("Calculating probabilities from time model for PTT ", names(dat$id)[1], "...")

  dat <- data.frame(dat)

  #calculate probabilities
  betas <- rstan::extract(mod, pars = 'betas')$betas
  b <- rstan::extract(mod, pars = 'b')$b
  strat <- unique(dat$strata)
  prob.all <- vector("list", length(strat))

  for (i in 1:length(strat)) {
    # print(i)

    #get parameters
    id1 <- dat[dat$strata == strat[i], "id1"][1]
    dist1 <- dat[dat$strata == strat[i], "step"][1]
    betas1 <- betas[,id1,]
    b1 <- b[,id1]

    #define design vector
    xmat <- as.matrix(cbind(1, dat[dat$strata == strat[i], covar.names]))
    #calculate mean using linear algebra
    mean1 <- dist1 * exp(xmat %*% t(betas1))
    #calculate a and b parameters of gamma distribution
    a1 <- b1 * t(mean1)

    #calculate probabilities for the realized and potential steps
    prob <- rep(NA, ncol(a1))
    for (j in 1:ncol(a1)) {
      #calculate posterior median of gamma density
      prob[j] <- median(dgamma(dat[dat$strata == strat[i], "dt"][j], a1[,j], b1))
    }

    prob.all[[i]] <- prob
  }

  p()  #plot progress bar

  prob.all <- unlist(prob.all)
  return(prob.all)
}

#----------------------------

### Function to test different sets of initial values in iterative manner
run.HMMs.internal = function(data, K, Par0, state.names, p, seeds) {
  #data = A data.frame that includes columns returned from momentuHMM::prepData()
  #K = integer. The number of states to be estimated by the HMM
  #Par0 = A named list that includes the initial parameters for each of the data streams as required for momentuHMM::fitHMM
  #state.names = A vector of names (i.e., character strings) to be assigned to each state
  #p = A function generated by progressr::progressor() to update progress bar
  #seeds = integer. A random seed to be applied when running separate HMMs in parallel

  set.seed(seeds)

  # Step length
  stepMean0 <- runif(K,
                     min = Par0$step[1:K] / 2,
                     max = Par0$step[1:K] * 2)
  stepSD0 <- runif(K,
                   min = Par0$step[(K+1):(K*2)] / 2,
                   max = Par0$step[1:K] * 2)
  whichzero_sl <- which(data$step == 0)
  propzero_sl <- length(whichzero_sl)/nrow(data)
  zeromass0_sl <- c(propzero_sl, rep(0, K-1))        #for zero distances by state



  # Fit model
  if(propzero_sl > 0) {  #don't include zero mass if no 0s present
    stepPar0 <- c(stepMean0, stepSD0, zeromass0_sl)
  } else {
    stepPar0 <- c(stepMean0, stepSD0)
  }

  anglePar0 <- Par0$angle


  hmm.res <- fitHMM(data = data, nbStates = K,
                    Par0 = list(step = stepPar0, angle = anglePar0),
                    dist = list(step = "gamma", angle = "wrpcauchy"),
                    formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                    estAngleMean = list(angle=TRUE),
                    stateNames = state.names
  )

  # Update progress bar
  p()

  return(hmm.res)
}


#-----------------------------------


run.HMMs = function(data, K, Par0, state.names, niter) {
  #data = A data.frame that includes columns returned from momentuHMM::prepData()
  #K = integer. The number of states to be estimated by the HMM
  #Par0 = A named list that includes the initial parameters for each of the data streams as required for momentuHMM::fitHMM
  #state.names = A vector of names (i.e., character strings) to be assigned to each state
  #niter = integer. The number of randomly perturbed initial values to use to determine the model w/ the global maximum likelihood

  # Convert data.frame into list of identical elements to fit w/ different initial values
  hmm.list <- lapply(1:niter, function(x) data)

  # Fit model
  plan(multisession, workers = availableCores() - 2)
  seeds <- future.apply::future_lapply(seq_along(hmm.list), FUN = function(x) sample(1:1e6, 1),
                                       future.chunk.size = Inf, future.seed = TRUE)  #set seed per list element

  handlers(handler_progress(incomplete=".", complete="*", current="o", clear = FALSE))
  progressr::with_progress({
    #set up progress bar
    p<- progressr::progressor(steps = length(hmm.list))

    hmm.res <- furrr::future_map2(hmm.list, seeds,
                                  ~run.HMMs.internal(data = .x,
                                                     K = K,
                                                     Par0 = Par0,
                                                     state.names = state.names,
                                                     p = p,
                                                     seeds = .y),
                                  .options = furrr_options(seed = TRUE))
  })

  plan(sequential)

  # Extract likelihoods of fitted models
  allnllk <- unlist(lapply(hmm.res, function(m) m$mod$minimum))

  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)

  # Best fitting model
  res <- hmm.res[[whichbest]]

  return(res)
}

#-----------------------------------

#' A function to obtain posterior means of marginals from an inla object
#'
#' This function allows you to get the posterior means of variances from a fitted inla object. This is useful because the output of inla is stored as precisions.
#' @param r.out Fitted inla() object
#' @keywords inla
#' @export
#' @examples
#' inla_emarginal()

inla_emarginal <- function(r.out){
  results <- sapply(r.out$marginals.hyperpar,
                    function(y)
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))

  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}


#-----------------------------------

#' A function to obtain posterior modes of marginals from an inla object
#'
#' This function allows you to get the posterior mode of variances from a fitted inla object. This is useful because the output of inla is stored as precisions.
#' @param r.out Fitted inla() object
#' @keywords inla
#' @export
#' @examples
#' inla_mmarginal()

inla_mmarginal <- function(r.out){
  results <- sapply(r.out$marginals.hyperpar,
                    function(y)
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))

  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}

#-----------------------------------


# Function that adds single row of NA coords to create gaps in tracks (when mapping w/ ggplot2)

add_track_gaps <- function(data, ind) {

  # find rows where gaps should be placed
  row.ind <- which(diff(data[[ind]]) > 0) + 1

  # insert NA for coords at each gap
  tmp <- data[row.ind,] %>%
    mutate(date = date - 30*60) %>%  #subtract 30 min from end of gap
    mutate(x = NA,
           y = NA)

  data.out <- rbind(data, tmp) %>%
    dplyr::arrange(date)

  return(data.out)
}
