
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
library(R2jags)
library(rstan)
library(MCMCvis)
library(arrow)


# dat <- vroom("Processed_data/Input for time model.csv", delim = ",")
dat <- read_parquet("Processed_data/Input for time model.parquet")

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
cor(dat2[,c('bathym.s','chla.s','kd490.s','sst.s')])
# PerformanceAnalytics::chart.Correlation(dat2[,c('bathym.s','chla.s','kd490.s','sst.s')])
#strong corr (0.99) between Chla and Kd490; omit Kd490 from subsequent analyses


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
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


dat3 <- dat2 %>%
  slice_sample(n = 10000)

# data
dat.list <- list(
  N = nrow(dat3),
  ID = n_distinct(dat3$id),
  dt = dat3$dt,
  dist = dat3$dist,
  bathym = dat3$bathym,
  chla = dat3$chla,
  sst = dat3$sst
  )

stan.model <- '
data {
  int N;                                  // sample size
  int ID;                                 // number of individuals
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst;                          // Sea surface temperature (C)
}


transformed data {
  vector[N] bathym_s;
  vector[N] chla_s;
  vector[N] sst_s;

  // Center and scale covariates
  bathym_s = (bathym - mean(bathym)) / sd(bathym);
  chla_s = (chla - mean(chla)) / sd(chla);
  sst_s = (sst - mean(sst)) / sd(sst);
}


parameters {
    real<lower=0> b;
    real b0;
    real b1;
    real b2;
    real b3;
}


transformed parameters {
  vector[N] mu;
  vector[N] a;

  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0 + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
    // calculate the corresponding a[i] parameter
    a[i] = mu[i] * b;
  }
}


model {
  // likelihood
  for (i in 1:N){
    // b[i] = exp(g0+g1*Miss[i]);

    // likelihood
    dt[i] ~ gamma(a[i], b);
  }


  // priors
  b0 ~ normal(0,0.01);
  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);

  b ~ exponential(1);
}


generated quantities {
  vector[N] y_hat;

  for (i in 1:N){
    // mean of gamma distribution
    y_hat[i] = dist[i] * exp(b0 + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
  }
}
'

mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 2022)

params <- c('b','b0','b1','b2','b3')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod1, params = params)

bayesplot::mcmc_acf(mod1, regex_pars = "^b")
bayesplot::mcmc_intervals(mod1, regex_pars = "^b")
bayesplot::mcmc_pairs(mod1, regex_pars = "^b")
bayesplot::mcmc_neff(neff_ratio(mod1, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)

#posterior prediction
y_hat <- vector("double", 4000)

mu_hat <- extract(mod1, par = 'mu')$mu
a_hat2 <- rowMeans(a_hat)
b_hat <- extract(mod1, par = 'b')$b %>%
  as.vector()

for (i in 1:4000) {
  y_hat[i] <- dgamma(1, shape = a_hat2[i], rate = b_hat[i])
}
ppc_dens_overlay(dat.list$dt,
                 extract(mod1, par = 'y_hat')$y_hat[1:200,])
