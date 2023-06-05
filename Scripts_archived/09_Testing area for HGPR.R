
dat <- rsf.pts_10s %>%
  # filter(id %in% c(181800, 181796, 159776)) %>%
  mutate(id1 = as.integer(factor(id)),
         id2 = id1,
         id3 = id1,
         id4 = id1)
covars <- c('log.bathym','log.npp','log.sst')
# x <- tmp$sst
y <- dat$obs / dat$wts2

# sigma0=sd(y)
pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))
ngroup <- n_distinct(dat$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)


bri.gpr <- function(x, y, pcprior, nbasis=5, degree=2, alpha=2, xout=x,
                    sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y))))
    stop("missing or infinite values in inputs are not allowed")
  # mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis), degree = degree, boundary = 'free')
  mesh.list <- vector("list", length(covars))

  for (i in 1:length(covars)) {  #buffer min and max values by exp(2)
    mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                       length.out = nbasis),
                                   degree = degree,
                                   boundary = 'free')
  }

  nu <-  alpha - 1/2
  kappa0 <- sqrt(8 * nu)/rho0
  tau0 <- 1 / (4 * kappa0^3 * sigma0^2)^0.5
  if(missing(pcprior)){

    spde.list <- vector("list", length(covars))
    for (i in 1:length(covars)) {
      spde.list[[i]] <- inla.spde2.matern(mesh.list[[i]], alpha=alpha, constr = FALSE,
                                          B.tau = cbind(log(tau0), 1, 0),
                                          B.kappa = cbind(log(kappa0), 0, 1),
                                          theta.prior.prec = 1e-4)
    }

  }else{

    spde.list <- vector("list", length(covars))
    for (i in 1:length(covars)) {
      spde.list[[i]] <-  inla.spde2.pcmatern(mesh.list[[i]],
                                             alpha=alpha,
                                             prior.range=c(pcprior[[i]][1],0.05),
                                             prior.sigma=c(pcprior[[i]][2],0.05))
    }
  }




  ############################
  ### Population-level GPs ###
  A.list.pop <- vector("list", length(covars))
  index.list.pop <- vector("list", length(covars))

  for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.list.pop[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                        loc=dat[[covars[[i]]]]
                                        )
    index.list.pop[[i]] <-  inla.spde.make.index(paste(covars[i], "pop", sep = "."),
                                                 n.spde = spde.list[[i]]$n.spde)
  }
  ############################


  ################################
  ### Include random GP slopes ###
  A.list.id <- vector("list", length(covars))
  index.list.id <- vector("list", length(covars))

  for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
    A.list.id[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                    loc=dat[[covars[[i]]]],
                                    group = dat$id1,
                                    n.group = ngroup,
                                    )
    index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
                                                n.spde = spde.list[[i]]$n.spde,
                                                n.group = ngroup)

    # mapper.A <- inlabru::bru_mapper_multi(list(
    #   main = inlabru::bru_mapper(mesh.list[[i]], indexed = TRUE),
    #   group = inlabru::bru_mapper_index(n = ngroup)
    # ))
    #
    # mapper.index <- inlabru::bru_mapper_multi(list(
    #   field.main = inlabru::bru_mapper(mesh.list[[i]], indexed = TRUE),
    #   field.group = inlabru::bru_mapper_index(n = ngroup)
    # ))
  }
  ################################


  ################################
  ### Include random GP slopes ###
  # A.list.id <- vector("list", length(covars))
  # index.list.id <- vector("list", length(covars))
  #
  # for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  #   A.list.id[[i]] <- inla.spde.make.A(mesh.list[[i]],
  #                                      loc=dat[[covars[[i]]]],
  #                                      repl = dat$id1
  #   )
  #   index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
  #                                               n.spde = spde.list[[i]]$n.spde,
  #                                               n.repl = max(dat$id1))
  # }
  ################################


  # index.list <- vector("list", length(covars))
  # for (i in 1:length(covars)) {
  #   index.list[[i]] <-  inla.spde.make.index(covars[i],
  #                                            n.spde = spde.list[[i]]$n.spde,
  #                                            n.group = ngroup)
  # }

  # Function to replicate pop A matrix to match ID-level
  # mat_expand <- function(A, n) {
  #   tmp <- do.call("cbind", rep(list(A), n))
  #   return(tmp)
  # }
  #
  # A.list.pop2 <- A.list.pop %>%
  #   map(., mat_expand, n = ngroup)
  # index.list.pop2 <- index.list.pop %>%
  #   map_depth(., 2, rep, 3)


  st.est <- inla.stack(data=list(y=y),
                       A=c(A.list.pop, A.list.id,
                                1, 1, 1),
                       effects=c(index.list.pop, index.list.id,
                                      list(Intercept = rep(1, nrow(dat)), id1 = dat$id1,
                                           log.sst = dat$log.sst)))
  # st.pred.bath <- inla.stack(data=list(y=NA),
  #                            A=A.list[which(1:length(A.list) %% 2 == 0)[1]],
  #                            effects=index.list[1],
  #                            tag="pred.bath")
  # st.pred.npp <- inla.stack(data=list(y=NA),
  #                           A=A.list[which(1:length(A.list) %% 2 == 0)[2]],
  #                           effects=index.list[2],
  #                           tag="pred.npp")
  # st.pred.sst <- inla.stack(data=list(y=NA),
  #                           A=A.list[which(1:length(A.list) %% 2 == 0)[3]],
  #                           effects=index.list[3],
  #                           tag="pred.sst")
  # st.all <- inla.stack(st.pop, st.id)
  formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +
    f(log.bathym.pop, model=spde.list[[1]]) +
    f(log.npp.pop, model=spde.list[[2]]) +
    f(log.sst.pop, model=spde.list[[3]]) +
    f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
    f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
    f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
    # f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid'), hyper = list(beta = gaus.prior)) +
    # f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid'), hyper = list(beta = gaus.prior)) +
    # f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid'), hyper = list(beta = gaus.prior)) +
    f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))



  formula <-  y ~ -1 + Intercept +
    f(log.bathym.pop, model=spde.list[[1]]) +
    f(log.npp.pop, model=spde.list[[2]]) +
    f(log.sst.pop, model=spde.list[[3]]) +
    f(log.bathym.id, model=spde.list[[1]], replicate = log.bathym.id.repl) +
    f(log.npp.id, model=spde.list[[2]], replicate = log.npp.id.repl) +
    f(log.sst.id, model=spde.list[[3]], replicate = log.sst.id.repl) +
    f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))





  data <-  inla.stack.data(st.est)
  set.seed(2023)
  tic()
  result <-  inla(formula, data=data,  family="poisson", weights = dat$wts2,
                  control.predictor = list(A=inla.stack.A(st.est),compute=TRUE),
                  control.fixed = list(
                    mean = c(6.592, -1),
                    prec = c(100, 100)),
                  num.threads = 1:1
  )
  toc()  #took 11 min to run single-threaded; took 7.5 min on laptop
  # pred.bath.ind <- inla.stack.index(sestpred, tag='pred.bath')$data
  # pred.npp.ind <- inla.stack.index(sestpred, tag='pred.npp')$data
  # pred.sst.ind <- inla.stack.index(sestpred, tag='pred.sst')$data
  # list(xout=xout,
  #      mean=result$summary.fitted.values$mean[ii],
  #      lcb=result$summary.fitted.values$"0.025quant"[ii],
  #      ucb=result$summary.fitted.values$"0.975quant"[ii],
  #      inlaobj=result)
}

# ind.list <- list(pred.bath.ind, pred.npp.ind, pred.sst.ind)

summary(result)



##############################
### Marginal effects plots ###
##############################

# Generate matrices for covariate raster data (for prediction)
A.test <- vector("list", length(covars))
newdat.list <- list(log.bathym = seq(min(dat$log.bathym), max(dat$log.bathym), length.out = 500),
                    log.npp = seq(min(dat$log.npp), max(dat$log.npp), length.out = 500),
                    log.sst = seq(min(dat$log.sst), max(dat$log.sst), length.out = 500))
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.test[[i]] <- inla.spde.make.A(mesh.list[[i]], loc=newdat.list[[covars[[i]]]])
}

# Replicate list elements for mapping over lists
A.test <- rep(A.test, each = max(dat$id1))

# Store resulting GP coeffs per covar into a list
pred.coeffs <- result$summary.random[4:6] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list.id[[1]]$log.bathym.id.group)) %>%
  map(., ~split(.x, .x$id)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(dat$id1)), rep(1:max(dat$id1), length(covars)), sep = "_"))

# Make predictions via linear algebra
pred.test <- A.test %>%
  map2(.x = ., .y = pred.coeffs,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  # bind_cols() %>%
  # rowSums()  #sum up all predictions across covars
  set_names(paste(rep(covars, each = max(dat$id1)), rep(1:max(dat$id1), length(covars)), sep = "_")) %>%
  bind_rows() %>%
  mutate(bathym = newdat.list$log.bathym,
         npp = newdat.list$log.npp,
         sst = newdat.list$log.sst) %>%
  mutate(across(everything(), exp)) %>%
  pivot_longer(cols = -c(bathym, npp, sst), names_to = "label", values_to = "mean") %>%
  separate(col = label, into = c('covar','id'), sep = "_")


# Wrangle data so that x and y values match up in long format
pred.test2 <- pred.test %>%
  mutate(x = case_when(str_detect(covar, "log.bathym") ~ bathym,
                       str_detect(covar, "log.npp") ~ npp,
                       str_detect(covar, "log.sst") ~ sst)) %>%
  dplyr::select(-c(bathym, npp, sst)) %>%
  mutate(covar = gsub(pattern = "log.", replacement = "", x = covar))


# Facet of all IDs and covars
ggplot() +
  geom_line(data = pred.test2, aes(x = x, y = log(mean), color = factor(id)), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~ covar, scales = "free")


# Depth
ggplot() +
  geom_line(data = pred.test2 %>%
              filter(covar == 'bathym'), aes(x = x, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(x = c(0,800), y = c(0,1.5e6)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


# NPP
ggplot() +
  geom_line(data = pred.test2 %>%
              filter(covar == 'npp'), aes(x = x / 1000, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(y = c(0,5000)) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


# SST
ggplot() +
  geom_line(data = pred.test2 %>%
              filter(covar == 'sst'), aes(x = x, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  lims(y = c(0,2.5e5)) +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


#--------------------------------------



# Make spatial predictions

# For ID 181796 in 2020-08
rsf.pts_10s %>%
  filter(id == 181796, obs == 1) %>%
  group_by(month.year) %>%
  count()  # most obs are in August

newdat.181796 <- data.frame(log.bathym = as.vector(terra::values(cov_list$bathym)) %>%
                       abs() %>%
                       log(),
                     log.npp = as.vector(terra::values(cov_list$npp$`2020-08-01`)) %>%
                       log(),
                     log.sst = as.vector(terra::values(cov_list$sst$`2020-08-01`)) %>%
                       log()) %>%
  mutate(row_id = 1:nrow(.)) %>%
  drop_na(log.bathym, log.npp, log.sst)
# names(newdat) <- c('log.bathym','log.bathym2','log.npp','log.npp2','log.sst','log.sst2')
summary(newdat.181796)

A.181796 <- vector("list", length(covars))
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.181796[[i]] <- inla.spde.make.A(mesh.list[[i]], loc = newdat.181796[[covars[[i]]]])
}


# Store resulting GP coeffs per covar into a list
ind <- dat[dat$id == 181796,]$id1[1]
pred.coeffs <- result$summary.random[4:6] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list.id[[1]]$log.bathym.id.group)) %>%
  map(., ~dplyr::filter(.x, id == ind)) %>%
  map(., pull, mean)


# Make predictions via linear algebra
x181796.pred <- A.181796 %>%
  map2(.x = ., .y = pred.coeffs,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  bind_cols() %>%
  rowSums()  #sum up all predictions across covars


# Define dummy raster for storing predictions
rast.pred.181796 <- cov_list$bathym
terra::values(rast.pred.181796) <- NA  # initially store all NAs for locs w/o predictions
terra::values(rast.pred.181796)[newdat.181796$row_id] <- x181796.pred


rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
names(rast.pred.181796.df)[3] <- "pred"
bbox <- ext(rast.pred.181796)

ggplot() +
  geom_raster(data = rast.pred.181796.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_path(data = rsf.pts_10s %>%
              filter(id == 181796, month.year == '2020-08-01', obs == 1), aes(x, y),
            color = "chartreuse", linewidth = 0.5) +
  labs(x="",y="", title = "ID 181796: August 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 20))
