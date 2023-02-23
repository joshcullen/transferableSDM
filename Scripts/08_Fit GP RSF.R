
### Fit RSF as GLMM w/ Gaussian Process Prior on SST ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_30 <- read_csv("Processed_data/GoM_Cm_RSFprep_30x.csv")
rsf.pts_50 <- read_csv("Processed_data/GoM_Cm_RSFprep_50x.csv")

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



####################
### Fit GLMM RSF ###
####################

# Center and scale covars; remove rows w/ incomplete observations
rsf.pts_10s <- rsf.pts_10 %>%
  drop_na(bathym, k490, npp, sst)
# obs.ind_10 <- which(rsf.pts_10s$obs == 1)
# rsf.pts_10s <- rsf.pts_10s %>%
#   mutate(bathym.s = (bathym - mean(bathym[obs.ind_10])) / sd(bathym),
#          k490.s = (k490 - mean(k490[obs.ind_10])) / sd(k490),
#          npp.s = (npp - mean(npp[obs.ind_10])) / sd(npp),
#          sst.s = (sst - mean(sst[obs.ind_10])) / sd(sst))
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(bathym.s = as.numeric(scale(bathym)),
         k490.s = as.numeric(scale(k490)),
         npp.s = as.numeric(scale(npp)),
         sst.s = as.numeric(scale(sst))
  )

rsf.pts_30s <- rsf.pts_30 %>%
  drop_na(bathym, k490, npp, sst)
obs.ind_30 <- which(rsf.pts_30s$obs == 1)
rsf.pts_30s <- rsf.pts_30s %>%
  mutate(bathym.s = (bathym - mean(bathym[obs.ind_30])) / sd(bathym),
         k490.s = (k490 - mean(k490[obs.ind_30])) / sd(k490),
         npp.s = (npp - mean(npp[obs.ind_30])) / sd(npp),
         sst.s = (sst - mean(sst[obs.ind_30])) / sd(sst))

rsf.pts_50s <- rsf.pts_50 %>%
  drop_na(bathym, k490, npp, sst)
obs.ind_50 <- which(rsf.pts_50s$obs == 1)
rsf.pts_50s <- rsf.pts_50s %>%
  mutate(bathym.s = (bathym - mean(bathym[obs.ind_50])) / sd(bathym),
         k490.s = (k490 - mean(k490[obs.ind_50])) / sd(k490),
         npp.s = (npp - mean(npp[obs.ind_50])) / sd(npp),
         sst.s = (sst - mean(sst[obs.ind_50])) / sd(sst))



# Infinitely-weighted logistic regression
rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, 5000, 1)
rsf.pts_30s$wts <- ifelse(rsf.pts_30s$obs == 0, 5000, 1)
rsf.pts_50s$wts <- ifelse(rsf.pts_50s$obs == 0, 5000, 1)


# Explore used vs available habitat values
rsf.pts_10s %>%
  # mutate(across(c(k490, npp), log)) %>%
  pivot_longer(cols = c(bathym, k490, npp, sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")

# Log-transform skewed covars to allow model fitting
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_30s <- rsf.pts_30s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_50s <- rsf.pts_50s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))


# Check Pearson corrs
cor(rsf.pts_10s[,c('bathym','k490','npp','sst')])  #NPP and K490 highly corr (0.78); remove k490


# Now explore transformed distributions
rsf.pts_10s %>%
  pivot_longer(cols = c(log.bathym, log.npp, log.sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)
rsf.pts_30s$wts2 <- ifelse(rsf.pts_30s$obs == 0, A / sum(rsf.pts_30s$obs == 0), 1e-6)
rsf.pts_50s$wts2 <- ifelse(rsf.pts_50s$obs == 0, A / sum(rsf.pts_50s$obs == 0), 1e-6)




## Mixed RSF via INLA
rsf.pts_10s$id1 <- as.numeric(factor(rsf.pts_10s$id))
rsf.pts_10s$id2 <- rsf.pts_10s$id1
rsf.pts_10s$id3 <- rsf.pts_10s$id1
rsf.pts_10s$id4 <- rsf.pts_10s$id1
rsf.pts_10s$id5 <- rsf.pts_10s$id1
rsf.pts_10s$id6 <- rsf.pts_10s$id1
rsf.pts_10s$id7 <- rsf.pts_10s$id1

rsf.pts_10s <- arrange(rsf.pts_10s, id1)


rsf.pts_30s$id1 <- as.numeric(factor(rsf.pts_30s$id))
rsf.pts_30s$id2 <- rsf.pts_30s$id1
rsf.pts_30s$id3 <- rsf.pts_30s$id1
rsf.pts_30s$id4 <- rsf.pts_30s$id1
rsf.pts_30s$id5 <- rsf.pts_30s$id1
rsf.pts_30s$id6 <- rsf.pts_30s$id1
rsf.pts_30s$id7 <- rsf.pts_30s$id1

rsf.pts_30s <- arrange(rsf.pts_30s, id1)


rsf.pts_50s$id1 <- as.numeric(factor(rsf.pts_50s$id))
rsf.pts_50s$id2 <- rsf.pts_50s$id1
rsf.pts_50s$id3 <- rsf.pts_50s$id1
rsf.pts_50s$id4 <- rsf.pts_50s$id1
rsf.pts_50s$id5 <- rsf.pts_50s$id1
rsf.pts_50s$id6 <- rsf.pts_50s$id1
rsf.pts_50s$id7 <- rsf.pts_50s$id1

rsf.pts_50s <- arrange(rsf.pts_50s, id1)

# create vector of ID values
id.vals <- unique(rsf.pts_10s$id1)



dat <- rsf.pts_10s #%>%
  # filter(id %in% c(181800, 181796, 159776))
covars <- c('log.bathym','log.npp','log.sst')
# x <- tmp$sst
y <- dat$obs / dat$wts2

sigma0=sd(y)
pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))



bri.gpr <- function(x, y, pcprior, nbasis=5, degree=2, alpha=2, xout=x,
                    sigma0=sd(y), rho0 = 0.25*(max(x) - min(x))){
  if (!all(is.finite(c(x, y))))
    stop("missing or infinite values in inputs are not allowed")
  # mesh <- inla.mesh.1d(seq(min(xout),max(xout),length.out = nbasis), degree = degree, boundary = 'free')
  mesh.list <- vector("list", length(covars))

  for (i in 1:length(covars)) {
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
    A.list[[i]] <- inla.spde.make.A(mesh.list[[ind[i]]], loc=dat[[covars[[ind[i]]]]])
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
    index.list[[i]] <-  inla.spde.make.index(covars[i], n.spde = spde.list[[i]]$n.spde)
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
  formula <-  y ~ -1 + Intercept + f(log.bathym, model=spde.list[[1]]) + f(log.npp, model=spde.list[[2]]) + f(log.sst, model=spde.list[[3]]) +
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

pred.vals <- vector("list", length(covars))
for (i in 1:length(covars)) {
  pred.vals[[i]] <- list(x=dat[[covars[i]]],
                         mean=result$summary.fitted.values$mean[ind.list[[i]]],
                         lcb=result$summary.fitted.values$"0.025quant"[ind.list[[i]]],
                         ucb=result$summary.fitted.values$"0.975quant"[ind.list[[i]]]) %>%
    bind_cols() %>%
    mutate(across(x:ucb, exp))
}


# Depth
ggplot() +
  # geom_ribbon(data = pred.vals[[1]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals[[1]], aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

# NPP
ggplot() +
  # geom_ribbon(data = pred.vals[[2]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals[[2]], aes(x = x / 1000, y = mean), linewidth = 1.5) +
  theme_bw() +
  # lims(x = c(0,300)) +
  labs(x = expression(paste("NPP (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

# SST
ggplot() +
  # geom_ribbon(data = pred.vals[[3]], aes(x = x, ymin = lcb, ymax = ucb), alpha = 0.5) +
  geom_line(data = pred.vals[[3]], aes(x = x, y = mean), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))









#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'k490', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('k490', 'npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list[["bathym"]][cov_list[["bathym"]] == 0.0000] <- NA


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')




####################################
### Generate predictive surfaces ###
####################################

rast.sep.20 <- data.frame(log.bathym = log(abs(terra::values(cov_list$bathym))) %>%
                            as.vector(),
                          log.npp = log(terra::values(cov_list$npp$`2020-09-01`)) %>%
                            as.vector(),
                          log.sst = log(terra::values(cov_list$sst$`2020-09-01`)) %>%
                            as.vector()) %>%
  mutate(row_id = 1:nrow(.)) %>%
  drop_na(log.bathym, log.npp, log.sst)



# Generate matrices for covariate raster data (for prediction)
A.pop.sep.20 <- vector("list", length(covars))
for (i in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.pop.sep.20[[i]] <- inla.spde.make.A(mesh.list[[i]], loc=rast.sep.20[[covars[[i]]]])
}

# Store resulting GP coeffs per covar into a list
pred.coeffs <- result$summary.random[-4] %>%  #remove random intercept term
  map(., ~pull(.x, mean))

# Make predictions via linear algebra
pred.sep.20 <- A.pop.sep.20 %>%
  map2(.x = ., .y = pred.coeffs,
       ~{.x %*% .y %>%
           as.vector()}
      ) %>%
  bind_cols() %>%
  rowSums()  #sum up all predictions across covars


# Define dummy raster for storing predictions
rast.pred <- cov_list$bathym
terra::values(rast.pred) <- NA  # initially store all NAs for locs w/o predictions
terra::values(rast.pred)[rast.sep.20$row_id] <- pred.sep.20
# plot(rast.pred)


rast.pred.df <- as.data.frame(rast.pred, xy = TRUE)
names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)


# Generate predictive surface for GP at pop-level
ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "Population Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20))
