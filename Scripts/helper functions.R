
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
      sf::st_as_sf(., coords = c('x','y'), crs = 3395)

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


avail_steps <- function(data, pop, nsteps) {

  # create list to store results for each strata
  avail.list <- vector("list", nrow(data))


  for (i in 1:nrow(data)) {

    ## generate random steps and angles
    rand.sl <- sample(pop$step, size = nsteps)
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
    furrr::future_map(~avail_steps(data = ., pop = data2, nsteps = nsteps),
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

#-----------------------------------

# Function to mask out land (or some polygon) from in-water UDs (or some other polygon)
# Per StackOverflow answer at https://stackoverflow.com/questions/71289669/intersection-keeping-non-intersecting-polygons-in-sf-r

st_mask <- function(layer, mask) {

  # First check if an intersection occurs
  tmp <- unlist(st_intersects(layer, mask))

  if (length(tmp > 0) & unique(st_geometry_type(layer)) %in% c("POLYGON","MULTIPOLYGON")) {
    intersection <- st_intersection(layer, mask)

    if (!is.null(names(intersection))) {  #when there's a df w/ column named 'geometry'
      diff1 <- st_difference(layer, st_union(intersection$geometry))
    } else {  #when object is simply an sf layer, not as df
      diff1 <- st_difference(layer, st_union(intersection))
    }


  } else if (length(tmp > 0) & unique(st_geometry_type(layer)) == "POINT") {
    diff1 <- layer[!st_intersects(layer, mask) %>% lengths > 0,]

  } else {  #for tracks/polygons where no intersection w/ land

    diff1 <- layer

  }

  return(diff1)
}

#-----------------------------------

# Updated version of the same function in {bayesmove}; fixes issue w/ specifying wrong N for ID/tseg combos
expand_behavior=function(dat, theta.estim, obs, nbehav, behav.names, behav.order) {

  #Assign behaviors (via theta) to each track segment
  theta.estim1<- apply(theta.estim[,1:nbehav], 1, function(x) x/sum(x)) %>%
    t()  #normalize probs for only first 3 behaviors being used
  theta.estim1<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim1)
  theta.estim1$id<- as.character(theta.estim1$id)
  names(theta.estim1)<- c("id", "tseg", behav.names)  #define behaviors

  nobs<- data.frame(id = obs$id, tseg = obs$tseg)
  nobs <- nobs %>%
    left_join(dat %>%
                dplyr::group_by(.data$id, .data$tseg) %>%
                dplyr::tally() %>%
                dplyr::ungroup(),
              by = c('id', 'tseg'))


  for (i in 1:nrow(theta.estim1)) {
    ind<- which(dat$id == theta.estim1$id[i] & dat$tseg == theta.estim1$tseg[i])

    if (i == 1) {
      theta.estim2<- rep(theta.estim1[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim1), byrow = TRUE)
    } else {
      tmp<- rep(theta.estim1[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim1), byrow = TRUE)

      theta.estim2<- rbind(theta.estim2, tmp)
    }
  }

  colnames(theta.estim2)<- names(theta.estim1)
  theta.estim2<- data.frame(theta.estim2, time1 = dat$time1, date = dat$date,
                            stringsAsFactors = FALSE)

  #Change col classes to numeric besides ID
  ind1<- which(names(theta.estim2) != "id")
  theta.estim2<- theta.estim2 %>%
    dplyr::mutate_at(names(theta.estim2)[ind1], as.numeric) %>%
    dplyr::select(.data$id, .data$tseg, .data$time1, .data$date, dplyr::everything())

  #Change into long format
  theta.estim.long<- tidyr::pivot_longer(theta.estim2, cols = -c(1:4),
                                         names_to = "behavior", values_to = "prop")
  theta.estim.long$behavior<- factor(theta.estim.long$behavior,
                                     levels = behav.names[behav.order])
  theta.estim.long<- theta.estim.long %>%
    dplyr::arrange(.data$behavior) %>%
    dplyr::mutate_at("date", lubridate::as_datetime)

  theta.estim.long
}

#-----------------------------------

# Updated version of the same function in {bayesmove}; fixes issue w/ function speed and efficiency
get_MAP = function(dat, nburn) {

  tmp <- dat[,(nburn + 2):ncol(dat)]  #subset only columns after burn-in period
  MAP.est <- apply(tmp, 1, function(x) which.max(x)) + nburn  #find max LML per ID

  return(MAP.est)
}

#-----------------------------------

#### internal function calculating predicted-to-expected ratio for each class-interval (from ecospat::ecospat.boyce())
boyce.internal <- function(interval, obs, fit) {
  pi <- sum(as.numeric(obs >= interval[1] & obs < interval[2])) / length(obs)
  ei <- sum(as.numeric(fit >= interval[1] & fit < interval[2])) / length(fit)

  tmp <- data.frame(f = round(pi/ei,10),  #F-ratio for predicted over expected obs
                    perc.use = pi)  #Percentage of total observations w/in a given bin
  return(tmp)
}


#-----------------------------------

# A modified version of ecospat::ecospat.boyce() that accounts for NAs from extracted habitat suitability values
# Also modified to accept vector of bin breaks
boyce <- function(fit, obs, nbins = 10, bin.method = c("seq","quantile"),
                PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman') {



  if (inherits(fit,"SpatRaster")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- extract(fit, obs, ID = FALSE) %>%
        pull(1)
      obs <- obs[!is.na(obs)]
    }
    fit <- terra::values(fit)
    fit <- fit[!is.na(fit)]
  }

  mini <- min(fit,obs)
  maxi <- max(fit,obs)

  if (bin.method == "quantile"){

    bins <- quantile(fit, seq(0, 1, length.out = nbins + 1))
      interval <- cbind(bins[-length(bins)], bins[-1])
      interval[length(bins) - 1,2] <- interval[length(bins) - 1,2] + 0.1

  } else if (bin.method == "seq") {

    bins <- seq(mini, maxi, length.out = nbins + 1)
    interval <- cbind(bins[-length(bins)], bins[-1])
    interval[length(bins) - 1,2] <- interval[length(bins) - 1,2] + 0.1
  } else {

  stop("Need to specify bin method as either 'seq' or 'quantile'")
}

  boyce.index.res <- apply(interval, 1, boyce.internal, obs, fit) %>%
    bind_rows()
  f <- boyce.index.res$f
  # f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], bins[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  # if(length(nclass)==1 & nclass == 0) {
  #   HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  # }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS, perc.use = boyce.index.res$perc.use))
}

#--------------------------------

# Slightly modified version of gbm::plot.gbm() to create partial dependency plots


gbm.pdp.vals <- function (x, i.var = 1, n.trees = x$n.trees, continuous.resolution = 100,
                          return.grid = TRUE, type = c("link", "response"), level.plot = TRUE,
                          contour = FALSE, number = 4, overlap = 0.1, col.regions = viridis::viridis,
                          ...)
{
  type <- match.arg(type)
  if (all(is.character(i.var))) {
    i <- match(i.var, x$var.names)
    if (any(is.na(i))) {
      stop("Requested variables not found in ", deparse(substitute(x)),
           ": ", i.var[is.na(i)])
    }
    else {
      i.var <- i
    }
  }
  if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
    warning("i.var must be between 1 and ", length(x$var.names))
  }
  if (n.trees > x$n.trees) {
    warning(paste("n.trees exceeds the number of tree(s) in the model: ",
                  x$n.trees, ". Using ", x$n.trees, " tree(s) instead.",
                  sep = ""))
    n.trees <- x$n.trees
  }
  if (length(i.var) > 3) {
    warning("plot.gbm() will only create up to (and including) 3-way ",
            "interaction plots.\nBeyond that, plot.gbm() will only return ",
            "the plotting data structure.")
    return.grid <- TRUE
  }
  grid.levels <- vector("list", length(i.var))
  for (i in 1:length(i.var)) {
    if (is.numeric(x$var.levels[[i.var[i]]])) {
      grid.levels[[i]] <- seq(from = min(x$var.levels[[i.var[i]]]),
                              to = max(x$var.levels[[i.var[i]]]), length = continuous.resolution)
    }
    else {
      grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]],
                                            levels = x$var.levels[[i.var[i]]])) - 1
    }
  }
  X <- expand.grid(grid.levels)
  names(X) <- paste("X", 1:length(i.var), sep = "")
  if (is.null(x$num.classes)) {
    x$num.classes <- 1
  }
  y <- .Call("gbm_plot", X = as.double(data.matrix(X)), cRows = as.integer(nrow(X)),
             cCols = as.integer(ncol(X)), n.class = as.integer(x$num.classes),
             i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees),
             initF = as.double(x$initF), trees = x$trees, c.splits = x$c.splits,
             var.type = as.integer(x$var.type), PACKAGE = "gbm")
  if (x$distribution$name == "multinomial") {
    X$y <- matrix(y, ncol = x$num.classes)
    colnames(X$y) <- x$classes
    if (type == "response") {
      X$y <- exp(X$y)
      X$y <- X$y/matrix(rowSums(X$y), ncol = ncol(X$y),
                        nrow = nrow(X$y))
    }
  }
  else if (is.element(x$distribution$name, c("bernoulli", "pairwise")) &&
           type == "response") {
    X$y <- 1/(1 + exp(-y))
  }
  else if ((x$distribution$name == "poisson") && (type == "response")) {
    X$y <- exp(y)
  }
  else if (type == "response") {
    warning("`type = \"response\"` only implemented for \"bernoulli\", ",
            "\"poisson\", \"multinomial\", and \"pairwise\" distributions. ",
            "Ignoring.")
  }
  else {
    X$y <- y
  }
  f.factor <- rep(FALSE, length(i.var))
  for (i in 1:length(i.var)) {
    if (!is.numeric(x$var.levels[[i.var[i]]])) {
      X[, i] <- factor(x$var.levels[[i.var[i]]][X[, i] +
                                                  1], levels = x$var.levels[[i.var[i]]])
      f.factor[i] <- TRUE
    }
  }
  names(X)[1:length(i.var)] <- x$var.names[i.var]
  if (return.grid) {
    return(X)
  }
  # nx <- length(i.var)
  # if (nx == 1L) {
  #   gbm:::plotOnePredictorPDP(X, ...)
  # }
  # else if (nx == 2) {
  #   gbm:::plotTwoPredictorPDP(X, level.plot = level.plot, contour = contour,
  #                       col.regions = col.regions, ...)
  # }
  # else {
  #   gbm:::plotThreePredictorPDP(X, nx = nx, level.plot = level.plot,
  #                         contour = contour, col.regions = col.regions, number = number,
  #                         overlap = overlap, ...)
  # }
}

#---------------------------

# Wrapper around gbm.pdp.vals() to collate predicted values for each covar

gbm.pdp <- function(fit, ...) {

  covars <- fit$var.names

  # Store predictions in list per covar
  tmp <- vector("list", length(covars)) %>%
    purrr::set_names(covars)

  for (i in 1:length(covars)) {
    tmp[[i]] <- gbm.pdp.vals(fit, i.var = i, ...)
  }

  tmp.df <- bind_rows(tmp) %>%
    pivot_longer(cols = -y, names_to = "covar", values_to = "x") %>%
    drop_na()

  return(tmp.df)
}

#----------------------------

# Function to fit hierarchical Gaussian Process regression

fit_hgpr <- function(data, covars, pcprior, mesh.seq, nbasis, degree, alpha, age.class,
                     int.strategy, method) {
  #data = data.frame containing all properly formatted columns for fitting the HGPR model
  #covars = vector of names for bathym, npp, and sst based on how used in formula expression
  #pcprior = stores rho_0 and sigma_0, respectively
  #mesh.seq = list of min and max values per covar for generating a sequence to create 1D mesh
  #nbasis = number of basis functions for approximating GP
  #degree = degree for defining 1D mesh of GP
  #alpha = for calculating Matern covariance matrix
  #age.class = logical; TRUE or FALSE whether the model should account for age class differences
  #int.strategy  = character for integration strategy to use from INLA (see ?control.inla); use 'auto' in most cases and 'eb' when model is failing
  #method = "corr" or "hybrid" for specifying whether model should be run as correlative or hybrid form


  # Define weighted response variable
  y <- data$obs / data$wts

  # Define groups by id1
  ngroup <- n_distinct(data$id1)

  # Define 1D mesh per covar
  mesh.list <- vector("list", length(covars))
  for (i in 1:length(covars)) {
    mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                       length.out = nbasis),
                                   degree = degree,
                                   boundary = 'free')
  }

  # Calculate Matern covariance matrix for SPDE (using PC priors)
  spde.list <- vector("list", length(covars))
  for (i in 1:length(covars)) {
    spde.list[[i]] <-  inla.spde2.pcmatern(mesh.list[[i]],
                                           alpha = alpha,
                                           prior.range = c(pcprior[[i]][1], 0.05),
                                           prior.sigma = c(pcprior[[i]][2], 0.05))
  }


  ################################
  ### Include random GP slopes ###
  A.list.id <- vector("list", length(covars))
  index.list.id <- vector("list", length(covars))

  for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.list.id[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                       loc = data[[covars[[i]]]],
                                       group = data$id1,
                                       n.group = ngroup
    )
    index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
                                                n.spde = spde.list[[i]]$n.spde,
                                                n.group = ngroup)
  }
  ################################




  # Structure remaining objects based on inclusion of age.class group (or not)
  if (age.class == FALSE) {

    # Define A matrix and index
    A.list.pop <- vector("list", length(covars))
    index.list.pop <- vector("list", length(covars))

    for (i in 1:length(covars)) {
      A.list.pop[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                          loc = data[[covars[[i]]]]
      )
      index.list.pop[[i]] <-  inla.spde.make.index(paste(covars[i], "pop", sep = "."),
                                                   n.spde = spde.list[[i]]$n.spde)
    }

    # Create INLA stack of A matrices and other data
    st.est <- inla.stack(data = list(y = y),
                         A = c(A.list.pop, A.list.id,
                               1, 1, 1),
                         effects = c(index.list.pop, index.list.id,
                                     list(Intercept = rep(1, nrow(data)), id1 = data$id1,
                                          log.sst = data$log.sst)))

    # Define formula for HGPR RSF model
    if (method == 'hybrid') {

      formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +  #fixed terms
        # pop-level terms
        f(log.bathym.pop, model=spde.list[[1]]) +
        f(log.npp.pop, model=spde.list[[2]]) +
        f(log.sst.pop, model=spde.list[[3]]) +
        # id-level terms
        f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group,
          control.group = list(model = 'iid')) +
        f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
        f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
        f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))

    } else if (method == 'corr') {

      formula <-  y ~ -1 + Intercept +  #fixed terms
        # pop-level terms
        f(log.bathym.pop, model=spde.list[[1]]) +
        f(log.npp.pop, model=spde.list[[2]]) +
        f(log.sst.pop, model=spde.list[[3]]) +
        # id-level terms
        f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group,
          control.group = list(model = 'iid')) +
        f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
        f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
        f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))

    }



  } else if (age.class == TRUE) {

    # Define groups by Age1
    ngroup.age <- n_distinct(data$Age1)

    # Define A matrix and index
    A.list.age <- vector("list", length(covars))
    index.list.age <- vector("list", length(covars))

    for (i in 1:length(covars)) {
      A.list.age[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                          loc = data[[covars[[i]]]],
                                          group = data$Age1,
                                          n.group = ngroup.age
      )
      index.list.age[[i]] <-  inla.spde.make.index(paste(covars[i], "age", sep = "."),
                                                   n.spde = spde.list[[i]]$n.spde,
                                                   n.group = ngroup.age)
    }


    # Create INLA stack of A matrices and other data
    st.est <- inla.stack(data = list(y = y),
                         A = c(A.list.age, A.list.id,
                               1, 1, 1),
                         effects = c(index.list.age, index.list.id,
                                     list(Intercept = rep(1, nrow(data)), id1 = data$id1,
                                          log.sst = data$log.sst)))


    # Define formula for HGPR RSF model
    if (method == 'hybrid') {

    formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +  #fixed terms
      # life stage-level terms
      f(log.bathym.age, model=spde.list[[1]], group = log.bathym.age.group,
        control.group = list(model = 'iid')) +
      f(log.npp.age, model=spde.list[[2]], group = log.npp.age.group, control.group = list(model = 'iid')) +
      f(log.sst.age, model=spde.list[[3]], group = log.sst.age.group, control.group = list(model = 'iid')) +
      # id-level terms
      f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
      f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
      f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
      f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))

    } else if (method == 'corr') {

      formula <-  y ~ -1 + Intercept +  #fixed terms
        # life stage-level terms
        f(log.bathym.age, model=spde.list[[1]], group = log.bathym.age.group,
          control.group = list(model = 'iid')) +
        f(log.npp.age, model=spde.list[[2]], group = log.npp.age.group, control.group = list(model = 'iid')) +
        f(log.sst.age, model=spde.list[[3]], group = log.sst.age.group, control.group = list(model = 'iid')) +
        # id-level terms
        f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
        f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
        f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
        f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))
    }

  }


  ## Run model
  stack.data <-  inla.stack.data(st.est)
  set.seed(2023)

if (method == "hybrid") {
  tic()
  hgpr.fit <- inla(formula, data=stack.data, family="Poisson", #Ntrials = 1,
                   control.predictor = list(A = inla.stack.A(st.est), compute = FALSE),
                   control.fixed = list(  #physiologically-informed SST component
                     mean = list(log.sst = 30*6.592, `I(log.sst^2)` = 30*-1),
                     prec = list(log.sst = 0.005, `I(log.sst^2)` = 0.005)),
                   weights = data$wts,
                   # num.threads = 1:1,
                   # inla.mode="experimental",
                   control.inla = list(int.strategy = int.strategy)
  )  #for greater reproducibility
  toc()

} else if (method == "corr") {
  tic()
  hgpr.fit <- inla(formula, data=stack.data, family="Poisson", #Ntrials = 1,
                   control.predictor = list(A = inla.stack.A(st.est), compute = FALSE),
                   weights = data$wts,
                   # num.threads = 1:1,
                   # inla.mode="experimental",
                   control.inla = list(int.strategy = int.strategy)
  )  #for greater reproducibility
  toc()

} else {
  stop("Must select either 'corr' or 'hybrid' for method")
}




  return(hgpr.fit)
}

#-----------------------------------

predict.hgpr <- function(cov_list, model_fit, covars, mesh.seq, nbasis, degree, age.class, method) {
  #cov_list = list of different environ covars as SpatRasters
  #model_fit = a fitted HGPR INLA model object
  #covars = character vector of covar names as used in model_fit
  #mesh.seq = list of covar ranges as used for model_fit
  #nbasis = number of basis functions used for model_fit
  #degree = degree value used for model_fit
  #alpha = alpha value used for model_fit
  #age.class = logical; TRUE or FALSE whether the model accounts for age class differences
  #method = "corr" or "hybrid" for specifying whether model should be run as correlative or hybrid form


  if (!method %in% c('hybrid','corr')) {
    stop("Must select either 'corr' or 'hybrid' for method")
  }

  # Define 1D mesh per covar
  mesh.list <- vector("list", length(covars))
  for (i in 1:length(covars)) {
    mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                       length.out = nbasis),
                                   degree = degree,
                                   boundary = 'free')
  }

  # Define vector of month.years for indexing
  my.ind <- names(cov_list$npp)

  # Set up rasters to store predictions
  if (age.class == FALSE) {
    rast.hgpr <- rep(cov_list$bathym, nlyr(cov_list$npp))
    names(rast.hgpr) <- my.ind
  } else if (age.class == TRUE) {
    rast.pred.juv <- rast.pred.adult <- rep(cov_list$bathym, nlyr(cov_list$npp))
    names(rast.pred.juv) <- names(rast.pred.adult) <- my.ind
  }





  # Make spatial predictions per month.year
  for (i in 1:nlyr(cov_list$npp)) {

    # Subset covars by month.year
    gp_vars <- data.frame(log.bathym = as.vector(terra::values(cov_list$bathym)) %>%
                            abs() %>%
                            log(),
                          log.npp = as.vector(terra::values(cov_list$npp[[my.ind[i]]])) %>%
                            log(),
                          log.sst = as.vector(terra::values(cov_list$sst[[my.ind[i]]])) %>%
                            log()) %>%
      mutate(row_id = 1:nrow(.)) %>%
      drop_na(log.bathym, log.npp, log.sst)

    if (method == "hybrid") {
      fixed_vars <- data.frame(Intercept = 1,
                               log.sst = as.vector(terra::values(cov_list$sst[[my.ind[i]]])) %>%
                                 log(),
                               log.sst2 = as.vector(terra::values(cov_list$sst[[my.ind[i]]])) %>%
                                 log() %>%
                                 . ^ 2) %>%
        mutate(row_id = 1:nrow(.)) %>%
        filter(row_id %in% gp_vars$row_id)
    }




    # Generate matrices for covariate raster data (for prediction)
    A.mat <- vector("list", length(covars))
    for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
      A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = gp_vars[[covars[[j]]]])
    }


    if (age.class == FALSE) {

      # Define coeff values from HGPR
      coeff1 <- model_fit$summary.random[1:3] %>%
        map(., ~pull(.x, mean))

      if (method == "hybrid") {
        # Define coeff values of fixed terms from HGPR
        coeff2 <- model_fit$summary.fixed$mean
      }


      # Make predictions on intensity of use from model for GP terms
      hgpr.pred <- A.mat %>%
        map2(.x = ., .y = coeff1,
             ~{.x %*% .y %>%
                 as.vector()}
        ) %>%
        bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
        rowSums()  #sum up all predictions across covars

      if (method == "hybrid") {
        # Make predictions using linear terms
        hgpr.pred2 <- as.matrix(vars2[,1:3]) %*% coeff2
      }

      # Store results in raster stack
      terra::values(rast.hgpr[[i]]) <- NA  # initially store all NAs for locs w/o predictions

      if (method == "hybrid") {
        terra::values(rast.hgpr[[i]])[gp_vars$row_id] <- hgpr.pred + hgpr.pred2[,1]
      } else {
        terra::values(rast.hgpr[[i]])[gp_vars$row_id] <- hgpr.pred
      }


    } else if (age.class == TRUE) {

      # Replicate list elements for mapping over lists
      A.mat2 <- rep(A.mat, each = 2)

      # Define coeff values from HGPR
      coeff1 <- model_fit$summary.random[1:3] %>%
        map(., ~dplyr::select(.x, mean)) %>%
        map(., ~mutate(.x, age = rep(1:2, each = 6))) %>%
        map(., ~split(.x, .x$age)) %>%
        flatten() %>%
        map(., pull, mean) %>%
        set_names(paste(rep(covars, each = 2), rep(1:2, length(covars)), sep = "_"))

      if (method == "hybrid") {
        # Define coeff values of fixed terms from HGPR
        coeff2 <- model_fit$summary.fixed$mean
      }


      # Make predictions via linear algebra
      hgpr.pred <- A.mat2 %>%
        map2(.x = ., .y = coeff1,
             ~{.x %*% .y %>%
                 as.vector()}
        ) %>%
        set_names(names(coeff1)) %>%
        bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
        split.default(., str_extract(names(.), "[0-9]$")) %>%  #split preds into list by age.class
        map(., rowSums) %>%  #sum up all predictions across covars
        set_names(c("Juv","Adult"))

      if (method == "hybrid") {
        # Make predictions using linear terms
        hgpr.pred2 <- as.matrix(vars2[,1:3]) %*% coeff2
      }


      # Store results in raster stack
      terra::values(rast.pred.juv[[i]]) <- terra::values(rast.pred.adult[[i]]) <- NA  # initially store all NAs for locs w/o predictions

      if (method == "hybrid") {
        terra::values(rast.pred.juv[[i]])[gp_vars$row_id] <- hgpr.pred$Juv + hgpr.pred2[,1]
        terra::values(rast.pred.adult[[i]])[gp_vars$row_id] <- hgpr.pred$Adult + hgpr.pred2[,1]
      } else {
        terra::values(rast.pred.juv[[i]])[gp_vars$row_id] <- hgpr.pred$Juv
        terra::values(rast.pred.adult[[i]])[gp_vars$row_id] <- hgpr.pred$Adult
      }

    }

  }


  if (age.class == FALSE) {
    return(rast.hgpr)
  } else if (age.class == TRUE) {
    return(list(Juv = rast.pred.juv,
                Adult = rast.pred.adult))
  }

}

#------------------------

#Function to normalize data (e.g., predicted raster values) via min-max scaling
# In this case, min and max values will be determined from entire set of predictions (provided by user)
# This function is intended to be used for values from SpatRaster stack
normalize <- function(x) {

  range <- values(x) %>%
    as.vector() %>%
    range(na.rm = TRUE)

  (x - range[1]) / (range[2] - range[1])
}


#----------------------------

# Function to fit hierarchical generalized linear model

fit_hglm <- function(data) {
  #data = data.frame containing all properly formatted columns for fitting the HGPR model
  #age.class = logical; TRUE or FALSE whether the model accounts for age class differences


  # create vector of ID values
  id.vals <- unique(data$id1)

  # Define formula expression
  RSF.formula <- obs ~ log.bathym + I(log.bathym ^ 2) + log.npp + I(log.npp ^ 2) + log.sst + I(log.sst ^ 2) +
    f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
    f(id2, log.bathym, values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
    f(id3, I(log.bathym ^ 2), values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
    f(id4, log.npp, values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
    f(id5, I(log.npp ^ 2), values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
    f(id6, log.sst, values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
    f(id7, I(log.sst ^ 2), values = id.vals, model = "iid",
      hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))


  # Fit model
  set.seed(2023)
  hglm.fit <- inla(RSF.formula, family = "Binomial", Ntrials = 1, data = data, #weights = data$wts,
              control.fixed = list(
                mean = 0,
                prec = list(default = 1e-3)),
              control.compute = list(waic = TRUE,
                                     dic = TRUE)#, num.threads = 1:1
  )

  return(hglm.fit)

}

#-----------------------------------

predict.hglm <- function(cov_list, model_fit, age.class) {
  #cov_list = list of different environ covars as SpatRasters
  #model_fit = a fitted HGPR INLA model object
  #age.class = logical; TRUE or FALSE whether the model accounts for age class differences


  # Define vector of month.years for indexing
  my.ind <- names(cov_list$npp)

  # Set up rasters to store predictions
  if (age.class == FALSE) {
    rast.pred <- rep(cov_list$bathym, nlyr(cov_list$npp))
    names(rast.pred) <- my.ind
  } else if (age.class == TRUE) {
    rast.pred.juv <- rast.pred.adult <- rep(cov_list$bathym, nlyr(cov_list$npp))
    names(rast.pred.juv) <- names(rast.pred.adult) <- my.ind
  }





  # Make spatial predictions per month.year
  for (i in 1:nlyr(cov_list$npp)) {

    # Subset covars by month.year
    vars <- data.frame(intercept = 1,
                       log.bathym = as.vector(terra::values(cov_list$bathym)) %>%
                         abs() %>%
                         log(),
                       log.bathym2 = as.vector(terra::values(cov_list$bathym)) %>%
                         abs() %>%
                         log() %>%
                         sapply(., function(x) x^2),
                       log.npp = as.vector(terra::values(cov_list$npp[[my.ind[i]]])) %>%
                         log(),
                       log.npp2 = as.vector(terra::values(cov_list$npp[[my.ind[i]]])) %>%
                         log() %>%
                         sapply(.,function(x) x^2),
                       log.sst = as.vector(terra::values(cov_list$sst[[my.ind[i]]])) %>%
                         log(),
                       log.sst2 = as.vector(terra::values(cov_list$sst[[my.ind[i]]])) %>%
                         log() %>%
                         sapply(., function(x) x^2))



    if (age.class == FALSE) {

      # Define coeff values
      coeff1 <- model_fit$summary.fixed$mean

      # Make predictions on intensity of use from model
      hglm.pred <- as.matrix(vars) %*% coeff1  #make predictions

      # Store results in raster stack
      terra::values(rast.pred[[i]]) <- hglm.pred


    } else if (age.class == TRUE) {

      stop("This doesn't exist yet")
    }

  }


  if (age.class == FALSE) {
    return(rast.pred)
  } else if (age.class == TRUE) {
    return(list(Juv = rast.pred.juv,
                Adult = rast.pred.adult))
  }

}

#-----------------------------------

# Function to preferentially sample locations at shallower depths to better resolve these functional responses
sample_rast_gradient <- function(rast, n_pts) {
  ## rast = SpatRaster; a layer where gradient will be rescaled from 0 to 1
  ## n_pts = integer; number of points to probabilistically sample

  # Rescale rast from 0 to 1 (to use as prob surface)
  probs <- scales::rescale(values(rast), to = c(0,1))
  probs <- ifelse(is.na(probs), 0, probs)  #convert NAs to 0

  # Sample cells from new probability values
  pts <- sample(x = 1:ncell(rast),
                size = n_pts,
                replace = TRUE,
                prob = probs
  )

  # Create DF of available point coordinates
  avail <- xyFromCell(rast, pts) %>%
    data.frame() #%>%
  # mutate(month.year = sample(names(cov_list$sst), size = n(), replace = T))

  return(avail)
}
