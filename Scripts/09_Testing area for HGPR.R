
dat <- rsf.pts_10s %>%
  filter(id %in% c(181800, 181796, 159776)) %>%
  mutate(id1 = as.numeric(factor(id)))
covars <- c('log.bathym','log.npp','log.sst')
# x <- tmp$sst
y <- dat$obs / dat$wts2

# sigma0=sd(y)
pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))
ngroup <- n_distinct(dat$id1)


bri.gpr <- function(x, y, pcprior, nbasis=5, degree=2, alpha=2, xout=x,
                    sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y))))
    stop("missing or infinite values in inputs are not allowed")
  # mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis), degree = degree, boundary = 'free')
  mesh.list <- vector("list", length(covars))

  for (i in 1:length(covars)) {  #buffer min and max values by exp(2)
    mesh.list[[i]] <- inla.mesh.1d(seq(min(dat[[covars[i]]]) - 2, max(dat[[covars[i]]]) + 2, length.out = nbasis),
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

  A.list <- vector("list", length(covars) * 2)
  ind <- rep(1:length(covars), each = 2)

  for (i in 1:(length(covars) * 2)) { #one matrix for model estimation and another for generating predictions for plotting
    A.list[[i]] <- inla.spde.make.A(mesh.list[[ind[i]]],
                                    loc=dat[[covars[[ind[i]]]]],
                                    group = dat$id1,
                                    n.group = ngroup
                                    )
  }



  ################################
  ### Include random GP slopes ###
  A.rand.bathym <- inla.spde.make.A(mesh.list[[ind[1]]], loc = dat[[covars[[ind[1]]]]],
                                    # index = seq_len(nrow(dat)),
                                    group = dat$id1, n.group = max(id.vals))
  index.rand.bathym <- inla.spde.make.index(paste("rand", covars[1], sep = "."), n.spde = spde.list[[1]]$n.spde,
                                            n.group = max(id.vals))
  ################################


  index.list <- vector("list", length(covars))
  for (i in 1:length(covars)) {
    index.list[[i]] <-  inla.spde.make.index(covars[i],
                                             n.spde = spde.list[[i]]$n.spde,
                                             n.group = ngroup)
  }
  st.est <- inla.stack(data=list(y=y),
                       A=append(A.list[1:length(A.list) %% 2 == 1],
                                list(1, 1)),
                       effects=append(index.list,
                                      list(Intercept = rep(1, nrow(dat)), id1 = dat$id1)),
                       tag="est")
  st.pred.bath <- inla.stack(data=list(y=NA),
                             A=A.list[which(1:length(A.list) %% 2 == 0)[1]],
                             effects=index.list[1],
                             tag="pred.bath")
  st.pred.npp <- inla.stack(data=list(y=NA),
                            A=A.list[which(1:length(A.list) %% 2 == 0)[2]],
                            effects=index.list[2],
                            tag="pred.npp")
  st.pred.sst <- inla.stack(data=list(y=NA),
                            A=A.list[which(1:length(A.list) %% 2 == 0)[3]],
                            effects=index.list[3],
                            tag="pred.sst")
  sestpred <- inla.stack(st.est, st.pred.bath, st.pred.npp, st.pred.sst)
  formula <-  y ~ -1 + Intercept +
    f(log.bathym, model=spde.list[[1]], group = log.bathym.group, control.group = list(model = 'iid')) +
    f(log.npp, model=spde.list[[2]], group = log.npp.group, control.group = list(model = 'iid')) +
    f(log.sst, model=spde.list[[3]], group = log.sst.group, control.group = list(model = 'iid')) +
    f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) #+
  # f(rand.log.bathym, model=spde.list[[1]])
  data <-  inla.stack.data(sestpred)
  tic()
  result <-  inla(formula, data=data,  family="poisson", weights = rep(dat$wts2, 4),
                  control.predictor= list(A=inla.stack.A(sestpred),compute=TRUE)#,
                  # control.fixed = list(
                  #   mean = 0,
                  #   prec = list(default = 1e-3))
  )
  toc()  #took 1 min to run
  pred.bath.ind <- inla.stack.index(sestpred, tag='pred.bath')$data
  pred.npp.ind <- inla.stack.index(sestpred, tag='pred.npp')$data
  pred.sst.ind <- inla.stack.index(sestpred, tag='pred.sst')$data
  # list(xout=xout,
  #      mean=result$summary.fitted.values$mean[ii],
  #      lcb=result$summary.fitted.values$"0.025quant"[ii],
  #      ucb=result$summary.fitted.values$"0.975quant"[ii],
  #      inlaobj=result)
}

ind.list <- list(pred.bath.ind, pred.npp.ind, pred.sst.ind)




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
pred.coeffs <- result$summary.random[-4] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list[[1]]$log.bathym.group)) %>%
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
  lims(x = c(0,800)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~id, scales = "free")


# NPP
ggplot() +
  geom_line(data = pred.test2 %>%
              filter(covar == 'npp'), aes(x = x / 1000, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~id, scales = "free")


# SST
ggplot() +
  geom_line(data = pred.test2 %>%
              filter(covar == 'sst'), aes(x = x, y = mean, color = factor(id)), linewidth = 1) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_wrap(~id, scales = "free")


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

#### PICK BACK UP FROM HERE ####


# Replicate list elements for mapping over lists
A.test <- rep(A.test, each = max(dat$id1))

# Store resulting GP coeffs per covar into a list
pred.coeffs <- result$summary.random[-4] %>%  #remove random intercept term
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, id = index.list[[1]]$log.bathym.group)) %>%
  map(., ~split(.x, .x$id)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = max(dat$id1)), rep(1:max(dat$id1), length(covars)), sep = "_"))


# Specify terms for ID 181796
ind <- which(unique(rsf.pts_10s$id) == 181796)
coeff.181796 <- random.coeffs2 %>%
  filter(ID == ind) %>%
  pull(mean)
x181796.pred <- as.matrix(newdat) %*% coeff.181796  #make predictions


rast.pred.181796 <- cov_list$bathym
terra::values(rast.pred.181796) <- exp(x181796.pred)


rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
names(rast.pred.181796.df)[3] <- "pred"

ggplot() +
  geom_raster(data = rast.pred.181796.df, aes(x, y, fill = log(pred))) +
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
