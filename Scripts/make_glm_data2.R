
# Transform data to be analyzed by Poisson GLM

make_glm_data = function(disc_path, neighbor_df, hex_cov_df, p){
  # disc_path = a data.frame containing the discretized path info for a single imputation of a PTT
  # neighbor_df = a data.frame defining the possible movements among neighboring hexgrid cells for a PTT
  # hex_cov_df = a data.frame containing the hexified environ covariate values for the quadrature time points of a given PTT


  # if leftover 'geometry' column from previous {sf} object, remove since this intereferes with use  of other functions
  disc_path <- disc_path %>%
    dplyr::select(-geometry)

  # join discretized paths, neighboring cells, and environ covariates together
  glm_data <- left_join(disc_path, neighbor_df, by = c("move_hex_from" = "hex_from")) %>%
    mutate(z = ifelse(!is.na(move) & hex_to == move_hex_to, 1, 0)) %>%
    left_join(hex_cov_df, by = c("move_hex_from" = "hex", "grid_time" = "Time"))

  if (n_distinct(na.omit(glm_data$move)) != sum(glm_data$z)) {
    stop("Not all moves are to adjacent cells. Need to adjust time step and/or grid size to prevent this from happening.")
  }


  # add a 'move' to end of track if one not present
  zorig <- glm_data$z[nrow(glm_data)]
  if(!zorig){
    glm_data$z[nrow(glm_data)] <- 1
    glm_data$move[(nrow(glm_data)-6+1):nrow(glm_data)] <- max(glm_data$move, na.rm=T)+1
  }

  # define previous bearing for each move
  prev_move_data <- glm_data %>%
    filter(z == 1) %>%
    dplyr::select(move, bearing_to_next) %>%
    mutate(prev_bearing = c(NA, bearing_to_next[-n()])) %>%
    dplyr::select(-bearing_to_next)

  # add previous bearing to main data.frame
  glm_data1 <- glm_data %>%
    left_join(prev_move_data, by = c("move")) %>%
    mutate(prev_bearing = zoo::na.locf(prev_bearing, fromLast = TRUE, na.rm = FALSE))

  if(!zorig){
    glm_data1$z[nrow(glm_data1)] <- 0
    glm_data1$move[(nrow(glm_data1)-6+1):nrow(glm_data1)] <- NA
  }

  # define distance of previous move on unit vector scale [-1,1]
  t1 <- glm_data1$datetime[glm_data1$z == 1 & glm_data1$move == 1]
  glm_data2 <- glm_data1 %>%
    mutate(prev_bearing = ifelse(datetime <= t1, NA, prev_bearing),
           prev_move_u = cos(pi * (90 - prev_bearing) / 180),
           prev_move_v = sin(pi * (90 - prev_bearing) / 180),
           prev_move = vector2scalar(prev_move_u, prev_move_v, bearing_to_next))

  glm_data2$prev_move = ifelse(is.na(glm_data2$prev_move), 0, glm_data2$prev_move)

  p()  #print progress bar

  return(glm_data2)
}


# Convert two vectors into a single scalar value where bearing is in degrees
vector2scalar = function(u, v, bearing, rotation_deg = NULL){
  unit_xy <- cbind(cos(pi * (90 - bearing) / 180), sin(pi * (90 - bearing) / 180))
  if (is.null(rotation_deg)) {
    uv <- cbind(u, v)
  } else {
    uv <- rotate_vector(u, v, rotation_deg)
  }
  scalar <- rowSums(unit_xy * uv)

  return(scalar)
}


# Function to rotate vectors for a given angle
rotate_vector = function(u, v, deg){
  theta <- deg * pi / 180
  R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  out <- t(t(R) %*% t(cbind(u,v)))

  return(out)
}
