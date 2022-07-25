
### Run time model to determine effect of covariates on movement rate

library(tidyverse)
library(lubridate)
library(furrr)
# library(vroom)
library(tictoc)
# library(R2jags)
library(rstan)
library(MCMCvis)
library(bayesplot)
library(arrow)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# dat <- vroom("Processed_data/Input for time model.csv", delim = ",")
dat <- read_parquet("Processed_data/Input for time model.parquet")

# Change time step from secs to mins
dat$dt <- dat$dt/60

# Retain only observed steps
dat <- dat %>%
  filter(obs == 1)

# Remove all observations w/ missing bathym values (since NA values were assigned to land)
# dat2 <- dat %>%
#   drop_na(bathym)
dat2 <- dat


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
cor(dat2[,c('bathym.s','chla.s','kd490.s','sst.s')] %>%
      drop_na())
# PerformanceAnalytics::chart.Correlation(dat2[,c('bathym.s','chla.s','kd490.s','sst.s')])
#strong corr (0.99) between Chla and Kd490; omit Kd490 from subsequent analyses



## Explore relationship between dt and each of the covars

ggplot(dat, aes(bathym, dt)) +
  geom_point() +
  theme_bw()
# will benefit from inclusion of a quadratic term

ggplot(dat, aes(Chla, dt)) +
  geom_point() +
  theme_bw()
# will benefit from inclusion of quadratic term

ggplot(dat, aes(SST, dt, color = id)) +
  geom_point() +
  theme_bw()
# will benefit from inclusion of quadratic term (especially varying by ID)


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


##############################
### Run time model w/ Stan ###
##############################

## Prior predictive simulation ##

set.seed(2022)
iter <- 50
dat2$id <- factor(dat2$id) %>%  #need to convert to integer for Stan
  as.numeric()
ID <- unique(dat2$id)
nID <- length(ID)

pp_list <- vector("list", iter)

for (i in 1:iter){
  print(i)

  # define priors
  mean1 <- rnorm(1, 0, 1)
  sd1 <- abs(rnorm(1, 0, 1))

  b0 <- vector("double", nID)
  b0 <- rnorm(nID, mean1, sd1)

  b1_p <- rnorm(1, 0, 1)
  b2_p <- rnorm(1, 0, 1)
  b3_p <- rnorm(1, 0, 1)

  b_p <- abs(rnorm(1, 0, 1))


  # mean of gamma distribution
  mu = dat2$dist * exp(b0[dat2$id] + b1_p*dat2$bathym.s + b2_p*dat2$chla.s + b3_p*dat2$sst.s)

  # calculate the corresponding a[i] parameter
  a = mu * b_p

  dt_p <- vector("double", length(a))  #store estimated 'dt'
  for (j in 1:length(a)) {
    dt_p[j] <- rgamma(1, a[j], b_p)
  }


  pp_list[[i]] <- dt_p
}
names(pp_list) <- 1:length(pp_list)



prior_pred <- bind_rows(pp_list, .id = "iter") %>%
  t() %>%
  data.frame() %>%
  mutate(iter = rownames(.)) %>%
  pivot_longer(cols = -c("iter"), names_to = "obs", values_to = "dt")



ggplot() +
  geom_line(data = prior_pred, aes(dt, group = iter), stat = "density", alpha = 0.1, color = "lightblue") +
  geom_density(data = dat2, aes(dt)) +
  xlim(0,1000) +
  theme_bw()
# there's an oddly high number of obs w/ dt approx 125 min

tmp <- dat2 %>%
  filter(dt > 100 & dt < 130)
table(tmp$id)
# common across IDs; this is likely an artefact of setting the imputed time interval to 2 hrs (120 min)





## Run model

# Do this step if not running the prior predictive simulation
dat2$id <- factor(dat2$id) %>%  #need to convert to integer for Stan
  as.numeric()

dat3 <- dat2 %>%
  slice_sample(n = 10000)

# Define objects pertaining to imputation of missing covars
bathym_missidx <- which(is.na(dat3$bathym))
n_bathym_miss <- length(bathym_missidx)
chla_missidx <- which(is.na(dat3$chla))
n_chla_miss <- length(chla_missidx)
sst_missidx <- which(is.na(dat3$sst))
n_sst_miss <- length(sst_missidx)


# data
dat.list <- list(
  N = nrow(dat3),
  ID = dat3$id,
  nID = n_distinct(dat3$id),
  dt = dat3$dt,
  dist = dat3$dist,
  bathym = scale(dat3$bathym)[,1],
  chla = scale(dat3$chla)[,1],
  sst = scale(dat3$sst)[,1],
  bathym2 = scale(dat3$bathym ^ 2)[,1],
  chla2 = scale(dat3$chla ^ 2)[,1],
  sst2 = scale(dat3$sst ^ 2)[,1],
  n_bathym_miss = n_bathym_miss,
  bathym_missidx = bathym_missidx,
  n_chla_miss = n_chla_miss,
  chla_missidx = chla_missidx,
  n_sst_miss = n_sst_miss,
  sst_missidx = sst_missidx
  )

# Replace missing covar values w/ Inf
dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
dat.list$chla <- ifelse(is.na(dat.list$chla), Inf, dat.list$chla)
dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)
dat.list$bathym2 <- ifelse(is.na(dat.list$bathym2), Inf, dat.list$bathym2)
dat.list$chla2 <- ifelse(is.na(dat.list$chla2), Inf, dat.list$chla2)
dat.list$sst2 <- ifelse(is.na(dat.list$sst2), Inf, dat.list$sst2)


stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
  int n_chla_miss;                        // number of missing Chla vals
  int chla_missidx[n_chla_miss];          // index for missing Chla vals
  int n_bathym_miss;                      // number of missing bathymetry vals
  int bathym_missidx[n_bathym_miss];      // index for missing bathymetry vals
  int n_sst_miss;                         // number of missing SST vals
  int sst_missidx[n_sst_miss];            // index for missing SST vals
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst;                          // Sea surface temperature (C)
  vector[N] bathym2;                      // Squared bathymetric depth (m)
  vector[N] chla2;                        // Squared chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst2;                         // Squared sea surface temperature (C)

}


parameters {
    real<lower=0> b;                      // Rate param for gamm distrib of likelihood
    vector[nID] b0_id;                    // Random intercepts by ID
    real bBathym;                         // Linear reg. coeff for bathymetry
    real bBathym2;                        // Quadratic reg. coeff for bathymetry
    real bChla;                           // Linear reg. coeff for chlorophyll-a
    real bChla2;                          // Quadratic reg. coeff for chlorophyll-a
    real bSST;                            // Linear reg. coeff for sea surface temp
    real bSST2;                           // Quadratic reg. coeff for sea surface temp

    real mean1;                           // Population mean of intercept
    real<lower=0> sd1;                    // Population sd of intercept

    vector[n_bathym_miss] bathym_miss;    // Index of missing bathymetry values
    vector[n_chla_miss] chla_miss;        // Index of missing chlorophyll-a values
    vector[n_sst_miss] sst_miss;          // Index of missing sea surface temp values
    vector[n_bathym_miss] bathym_miss2;   // Index of squared missing bathymetry values
    vector[n_chla_miss] chla_miss2;       // Index of squared missing chlorophyll-a values
    vector[n_sst_miss] sst_miss2;         // Index of squared missing sea surface temp values

    real nu_bathym;                       // Mean value of standardized bathymetry
    real nu_chla;                         // Mean value of standardized chlorophyll-a
    real nu_sst;                          // Mean value of standardized sea surface temp
    real<lower=0> sigma_bathym;           // SD of bathymetry values
    real<lower=0> sigma_chla;             // SD of chlorophyll-a values
    real<lower=0> sigma_sst;              // SD of sea surface temp values
}


transformed parameters {
  vector[N] mu;
  real<lower=0> a[N];
  vector[N] bathym_merge;
  vector[N] chla_merge;
  vector[N] sst_merge;
  vector[N] bathym_merge2;
  vector[N] chla_merge2;
  vector[N] sst_merge2;


// Merge imputed w/ observed covariate vals
  bathym_merge = bathym;
    for (i in 1:n_bathym_miss) {
        bathym_merge[bathym_missidx[i]] = bathym_miss[i];
    }

  chla_merge = chla;
    for (i in 1:n_chla_miss) {
        chla_merge[chla_missidx[i]] = chla_miss[i];
    }

  sst_merge = sst;
    for (i in 1:n_sst_miss) {
        sst_merge[sst_missidx[i]] = sst_miss[i];
    }

  bathym_merge2 = bathym2;
    for (i in 1:n_bathym_miss) {
        bathym_merge2[bathym_missidx[i]] = bathym_miss2[i];
    }

  chla_merge2 = chla2;
    for (i in 1:n_chla_miss) {
        chla_merge2[chla_missidx[i]] = chla_miss2[i];
    }

  sst_merge2 = sst2;
    for (i in 1:n_sst_miss) {
        sst_merge2[sst_missidx[i]] = sst_miss2[i];
    }


// Calculate params for gamma distrib
  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym_merge[i] + bBathym2*bathym_merge2[i] +
    bChla*chla_merge[i] + bChla2*chla_merge2[i] + bSST*sst_merge[i] + bSST2*sst_merge2[i]);

    // calculate the corresponding a[i] parameter
    a[i] = mu[i] * b;
  }
}


model {
  // likelihood
  for (i in 1:N){
    dt[i] ~ gamma(a[i], b);
  }


  // priors
  for (j in 1:nID){
    b0_id[j] ~ normal(mean1,sd1);
  }

  mean1 ~ normal(0,1);
  sd1 ~ normal(0,1);

  [sigma_chla, sigma_bathym, sigma_sst] ~ normal(0,2);
  [nu_chla, nu_bathym, nu_sst, b, bBathym, bBathym2, bChla, bChla2, bSST, bSST2] ~ normal(0, 1);
  chla_merge ~ normal(nu_chla, sigma_chla);
  bathym_merge ~ normal(nu_bathym, sigma_bathym);
  sst_merge ~ normal(nu_sst, sigma_sst);

  bathym_merge2 ~ normal(0,1);
  chla_merge2 ~ normal(0,1);
  sst_merge2 ~ normal(0,1);

}


generated quantities {
  vector[N] y_hat;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym_merge[i] + bBathym2*bathym_merge2[i] +
    bChla*chla_merge[i] + bChla2*chla_merge2[i] + bSST*sst_merge[i] + bSST2*sst_merge2[i]);
  }
}
'

mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 8675309)

params <- c('b','b0_id','bBathym','bBathym2','bChla','bChla2','bSST','bSST2','mean1','sd1','nu_chla','sigma_chla','nu_bathym',
            'sigma_bathym','nu_sst','sigma_sst')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod1, params = params)

# posterior predictive check
ppc_dens_overlay(dat3$dt,
                 rstan::extract(mod1, pars = 'y_hat')$y_hat[1:200,]) +
  xlim(0,1000)


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
