
### Fit HGPR RSF at different spatial scales ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10_5km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_10_10km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_10km.csv")
rsf.pts_10_20km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_20km.csv")
rsf.pts_10_40km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_40km.csv")
rsf.pts_10_80km <- read_csv("Processed_data/GoM_Cm_RSFprep_10x_80km.csv")

rsf.list <- list(
  sc.5 = rsf.pts_10_5km,
  sc.10 = rsf.pts_10_10km,
  sc.20 = rsf.pts_10_20km,
  sc.40 = rsf.pts_10_40km,
  sc.80 = rsf.pts_10_80km
)

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10_5km)
summary(rsf.pts_10_5km)



#####################
### Fit HGPR RSFs ###
#####################

# Remove rows w/ incomplete observations; log-transform covars
rsf.list2 <- rsf.list %>%
  map(., ~{.x %>%
      drop_na(bathym, npp, sst) %>%
      mutate(log.bathym = log(abs(bathym)),
             log.npp = log(npp),
             log.sst = log(sst))}
      )

#### PICK BACK UP FROM HERE ####


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)




# Add ID as integer for INLA
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(id1 = as.integer(factor(id)))

# Specify model params, predictors, and data
covars <- c('log.bathym','log.npp','log.sst')
y <- rsf.pts_10s$obs / rsf.pts_10s$wts2

pcprior <- list(bathym = c(1,10), npp = c(1,10), sst = c(1,10))  #stores rho_0 and sigma_0, respectively
ngroup <- n_distinct(rsf.pts_10s$id1)
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)

nbasis <- 5  #number of basis functions for approximating GP
degree <- 2  #degree for defining 1D mesh of GP
alpha <- 2  #for calculating Matern covariance matrix


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



############################
### Population-level GPs ###
A.list.pop <- vector("list", length(covars))
index.list.pop <- vector("list", length(covars))

for (i in 1:length(covars)) {
  A.list.pop[[i]] <- inla.spde.make.A(mesh.list[[i]],
                                      loc = rsf.pts_10s[[covars[[i]]]]
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
                                     loc = rsf.pts_10s[[covars[[i]]]],
                                     group = rsf.pts_10s$id1,
                                     n.group = ngroup,
  )
  index.list.id[[i]] <-  inla.spde.make.index(paste(covars[i], "id", sep = "."),
                                              n.spde = spde.list[[i]]$n.spde,
                                              n.group = ngroup)
}
################################


# Create INLA stack of A matrices and other data
st.est <- inla.stack(data = list(y = y),
                     A = c(A.list.pop, A.list.id,
                           1, 1, 1),
                     effects = c(index.list.pop, index.list.id,
                                 list(Intercept = rep(1, nrow(rsf.pts_10s)), id1 = rsf.pts_10s$id1,
                                      log.sst = rsf.pts_10s$log.sst)))

# Define formula for HGPR RSF model
formula <-  y ~ -1 + Intercept + log.sst + I(log.sst^2) +  #fixed terms
  # pop-level terms
  f(log.bathym.pop, model=spde.list[[1]]) +
  f(log.npp.pop, model=spde.list[[2]]) +
  f(log.sst.pop, model=spde.list[[3]]) +
  # id-level terms
  f(log.bathym.id, model=spde.list[[1]], group = log.bathym.id.group, control.group = list(model = 'iid')) +
  f(log.npp.id, model=spde.list[[2]], group = log.npp.id.group, control.group = list(model = 'iid')) +
  f(log.sst.id, model=spde.list[[3]], group = log.sst.id.group, control.group = list(model = 'iid')) +
  f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE)))


## Run model
data <-  inla.stack.data(st.est)
set.seed(2023)
tic()
hgpr.fit <- inla(formula, data=data, family="poisson", weights = rsf.pts_10s$wts2,
                 control.predictor = list(A = inla.stack.A(st.est), compute = TRUE),
                 control.fixed = list(  #physiologically-informed SST component
                   mean = list(log.sst = 10*6.592, `I(log.sst^2)` = 10*-1),
                   prec = list(log.sst = 0.1, `I(log.sst^2)` = 0.1)),
                 num.threads = 1:1)  #for greater reproducibility
toc()  #took 12.5 min to run single-threaded

summary(hgpr.fit)
