
# Function to discretize continuous animal trajectories for CT(S)MC model

make_disc_path_internal = function(sims1, hex_poly, quad_times, p) {
  # sims: a {sf} POINT object with columns for 'ptt', 'rep', and 'datetime' for a single indiv.
  # hex_poly: a {sf} POLYGON object representing the hexgrid for a particular ID
  # quad_times: a vector of POSIXct dates corresponding w/ the environmental rasters

  # Change 'rep' from character string to integer
  sims1$rep <- str_replace(sims1$rep, pattern = "\\d+_*", replacement = "") %>%  #convert rep to int
    as.integer()

  # Change datetime to match format for environ raster quadrature points; set to 1st of each month
  # sims1 <- sims1 %>%
  #   mutate(date = as_date(datetime)) %>%
  #   mutate(date = str_replace(date, pattern = "-..$", replacement = "-01"))

  bound = list(min(sims1$date), max(sims1$date))
  cont_path_time <- sims1 %>%
    # filter(rep == 1) %>%
    dplyr::select(datetime, date)

  ## Add hex ids
  disc_path <- data.frame(cont_path_time,
                          move_hex_to = unlist(st_intersects(sims1, hex_poly))
                          ) %>%
    mutate(move_hex_from = c(NA, move_hex_to[-n()])) %>%
    left_join(.,
              data.frame(date = quad_times,
                         cov_trans = 1),
              by = "date") %>%
    dplyr::arrange(datetime)

  # Make quadrature times closed on the right (not left)
  disc_path1 <- disc_path %>%
    mutate(grid_time = ifelse(!is.na(cov_trans), date, NA) %>%
             zoo::na.locf()) %>%
    # mutate(grid_time = ifelse(!is.na(cov_trans), NA, grid_time) %>% zoo::na.locf(na.rm=F) %>% crawl::intToPOSIX()) %>%
    mutate(across(c("date", "grid_time"), as_date)) %>%
    filter(date >= bound[[1]])

  # Add delta and time since movement (tsm)
  disc_path2 <- disc_path1 %>%
    mutate(elapsed_time = c(0, cumsum(diff(as.numeric(datetime)))) / 3600,  #in hours
           move = case_when(move_hex_from != move_hex_to ~ 1,
                            TRUE ~ 0,
                            is.na(move_hex_from) ~ 0) %>%
             cumsum() %>%
             ifelse(move_hex_from == move_hex_to, NA, .),
           tempA = zoo::na.locf(move, fromLast = TRUE, na.rm = FALSE),
           tempB = c(0, diff(elapsed_time))) %>%
    group_by(tempA) %>%
    mutate(tsm = cumsum(tempB) %>%  #in hours
             ifelse(. == 0, NA, .)) %>%
    ungroup() %>%
    dplyr::select(-contains("temp"))


  disc_path3 <- disc_path2 %>%
    mutate(move_hex_to = zoo::na.locf(move_hex_to),
           move_hex_from = ifelse(is.na(move_hex_from), move_hex_to, move_hex_from),
           delta = c(NA, diff(elapsed_time))) %>%
    mutate(rep = unique(sims1$rep), .before = "datetime")

  p()  #print progress bar

  return(disc_path3)
}

#---------------------------------
make_disc_path = function(sims, hex_grid, quad_times) {
  # sims: a {sf} POINT object with columns for 'ptt', 'rep', and 'datetime' for all indivs
  # hex_grid: a list of {sf} POLYGON objects representing the hexgrid for each ID, where each element is a list w/ element named 'poly'
  # quad_times: a vector of POSIXct dates corresponding w/ the environmental rasters


  message("Filtering locations outside of quad_times...")
  sims <- sims %>%
    mutate(date = as_date(datetime) %>%
             str_replace(pattern = "-..$", replacement = "-01") %>%
             as_date()) %>%
    filter(date >= min(quad_times) & date <= max(quad_times))
  PTTs <- unique(sims$ptt)

  hex_grid1 <- hex_grid[as.character(PTTs)]

  num_reps <- str_replace(sims$rep, pattern = "\\d+_*", replacement = "") %>%  #convert rep to int
    as.integer() %>%
    max()

  # hex_cov_df = readRDS(hex_df)
  # num_reps = max(sims$reps)

  # create storage data frame
  # fit_frame <- expand.grid(ptt = PTTs, rep = 1:num_reps)
  # quad_times <- unique(hex_cov_df$Time)

  # perform fitting in parallel over imputations
  fit_list <- vector("list", length(PTTs)) %>%
    set_names(PTTs)

  # to generate progress bar
  # pb <- progress_bar$new(
    # format = "  downloading [:bar] :percent eta: :eta",
    # total = length(PTTs), clear = FALSE, width = 60)



  # spatially discretize path
  for (i in 1:length(PTTs)) {

    message("Starting ", PTTs[i],"...")

    # define aesthetics of progress bar
    progressr::handlers(progressr::handler_progress(incomplete = ".", complete = "*",
                                                    current = "o", clear = FALSE))

    sims_list <- sims[sims$ptt == PTTs[i],] %>%
      split(.$rep)

    progressr::with_progress({
      #set up progress bar
      p<- progressr::progressor(steps = length(sims_list))

      disc_path_list <- future_map(sims_list,
                                   ~make_disc_path_internal(., hex_grid1[[i]]$poly, quad_times, p),
                                   .options = furrr_options(seed = 2022)
      )
    })

    fit_list[[i]] <- bind_rows(disc_path_list) %>%
      mutate(ptt = PTTs[i], .before = rep)

    # pb$tick()  #update progress bar

  }


  message("Merging list into data.frame...")
  fit_df <- bind_rows(fit_list)

  return(fit_df)
}
