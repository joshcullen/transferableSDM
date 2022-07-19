
### Run time model to determine effect of covariates on movement rate

# library(bayesmove)
library(tidyverse)
# library(sf)
# library(terra)
library(lubridate)
library(furrr)
# library(future)
library(vroom)
library(tictoc)
# library(R2jags)
library(rstan)
library(MCMCvis)


dat <- vroom("Processed_data/Input for time model.csv", delim = ",")

# Change time step from secs to mins
dat$dt <- dat$dt/60

# Retain only observed steps
dat <- dat %>%
  filter(obs == 1)

# Remove all observations w/ missing environ covar values
dat2 <- dat %>%
  drop_na(bathym, Chla, Kd490, SST)


# Center and scale covariates
dat2 <- dat2 %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         chla.s = scale(Chla) %>%
           as.vector(),
         kd490.s = scale(Kd490) %>%
           as.vector(),
         sst.s = scale(SST) %>%
           as.vector())


# Check correlation among covars
PerformanceAnalytics::chart.Correlation(dat2[,c('bathym.s','chla.s','kd490.s','sst.s')])
#strong corr (1.00) between Chla and Kd490; omit Kd490 from subsequent analyses


### Run time model w/ JAGS ###

model <- function(){

  #likelihood
  for (i in 1:nobs){
    #mean of gamma distribution
    mu[i] <- dist[i] * exp(b0 + b1*bathym[i] + b2*chla[i] + b3*sst[i])
    # b[i] <- exp(g0+g1*Miss[i])
    #calculate the corresponding a[i] parameter
    a[i] <- mu[i] * b
    #likelihood
    dt[i] ~ dgamma(a[i], b)
  }


  #priors
  b0 ~ dnorm(0,0.01)
  b1 ~ dnorm(0,1)
  b2 ~ dnorm(0,1)
  b3 ~ dnorm(0,1)

  b ~ dexp(1)
  # g0 ~ dnorm(0,0.01)
  # g1 ~ dnorm(0,0.01)
}

# data
nobs <- nrow(dat2)
dt <- dat2$dt
dist <- dat2$dist
bathym <- dat2$bathym.s
chla <- dat2$chla.s
sst <- dat2$sst.s
dat.list <- list(nobs=nobs, dt=dt, dist=dist, bathym=bathym, chla=chla, sst=sst)


#set parameters to track
params <- c('b','b0','b1','b2','b3')



## run model

n.iter <- 5000  #number of iterations per chain
n.thin <- 10  #how to thin MCMC results
n.burnin <- n.iter / 2  #number of iterations to discard as burn-in
n.chains <- 3  #number of MCMC chains


tic()
res <- jags.parallel(model.file = model, parameters.to.save = params, data = dat.list,
                    n.chains = 3, n.burnin = 2500, n.iter = 5000,
                    n.thin = 10, DIC = TRUE, jags.seed = 123)
toc()
# takes 45 min to run 5000 iterations


res

MCMCsummary(res)
MCMCtrace(res, ind = TRUE, iter = 750, pdf = FALSE)
par(mfrow=c(1,1))
MCMCplot(res, excl = "deviance")

res.summ<- res$BUGSoutput$summary



### Run time model w/ Stan ###

# data
N <- nrow(dat2)
dt <- dat2$dt
dist <- dat2$dist
bathym <- dat2$bathym.s
chla <- dat2$chla.s
sst <- dat2$sst.s
dat.list <- list(nobs=nobs, dt=dt, dist=dist, bathym=bathym, chla=chla, sst=sst)

stan.model <- '
data {
  int<lower=1> N;                         // sample size
  int<lower=1> ID;                        // number of individuals
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst                           // Sea surface temperature (C)
}


transformed data {
  real bathym.s;
  real chla.s
  real sst.s

  // Center and scale covariates
  bathym.s = mean(bathym) / sd(bathym);
  chla.s = mean(chla) / sd(chla);
  sst.s = mean(sst) / sd(sst)
}
'
