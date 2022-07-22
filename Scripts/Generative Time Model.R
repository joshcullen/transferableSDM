

## Simulate data (time interval) from generative model ##

# mu = dist * exp(b0_id + b1*bathym + b2*chla + b3*sst)
# a = mu * b

# dt ~ gamma(a, b)

library(rstan)
library(tidyverse)
library(MCMCvis)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



##############################
### Random intercept model ###
##############################

set.seed(2022)

# Define data
ID <- rep(1:5, each = 1000)
nID <- n_distinct(ID)
N <- length(ID)
dist <- rgamma(N, 15, 0.2)

bathym <- runif(N, min = -100, max = 0)
chla <- rnorm(N, 350, 10)
sst <- rnorm(N, 28, 2)

bathym.s <- scale(bathym)
chla.s <- scale(chla)
sst.s <- scale(sst)

# Define coefficients
mean1 <- -0.4
sd1 <- 1
b0_id <- rnorm(nID, mean = mean1, sd = sd1)
b1 <- -0.5
b2 <- 0
b3 <- 0.5
b <- 0.05

# Generate time interval (dt) using subset of real data for `dist` and environ covars
mu <- dist * exp(b0_id[ID] + b1*bathym.s + b2*chla.s + b3*sst.s)
a <- mu * b
dt <- rgamma(N, shape = a, rate = b)





## Prior predictive simulation ##

plot(density(dt), lwd = 2, col = "red")

iter <- 1000
n <- 5000

pp_list <- vector("list", iter)

for (i in 1:iter){
  # define priors
  mean2 <- rnorm(1, 0, 1)
  sd2 <- abs(rnorm(1, 0, 1))

  b0 <- vector("double", nID)
  b0 <- rnorm(nID, mean2, sd2)

  b1_p <- rnorm(n, 0, 1)
  b2_p <- rnorm(n, 0, 1)
  b3_p <- rnorm(n, 0, 1)

  b_p <- abs(rnorm(1, 0, 1))

  # mean of gamma distribution
  mu = dist * exp(b0[ID] + b1*bathym.s + b2*chla.s + b3*sst.s)
  # calculate the corresponding a[i] parameter
  a = mu * b
  dt_p <- rgamma(n, a, b)

  pp_list[[i]] <- dt_p
}

for (i in 1:length(pp_list)) {
  lines(density(pp_list[[i]]), lwd = 0.25, xlim = c(0,1000))
}

lines(density(dt), lwd = 2, col = "red")




dat.list <- list(
  N = N,
  ID = ID,
  nID = nID,
  dt = dt,
  dist = dist,
  bathym = bathym,
  chla = chla,
  sst = sst
)

stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
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
    vector[nID] b0_id;
    real b1;
    real b2;
    real b3;

    real mean1;
    real<lower=0> sd1;
}


transformed parameters {
  vector[N] mu;
  real<lower=0> a[N];

  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
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

  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);

  b ~ normal(0,1);

  mean1 ~ normal(0,1);
  sd1 ~ normal(0,1);
}


generated quantities {
  vector[N] y_hat;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_s[i] + b2*chla_s[i] + b3*sst_s[i]);
  }
}
'



mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
             seed = 8675309)

params <- c('b','b0_id','b1','b2','b3','mean1','sd1')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod1, params = params)

bayesplot::mcmc_acf(mod1, regex_pars = "^b", pars = c("mean1","sd1"))
bayesplot::mcmc_intervals(mod1, regex_pars = "^b", pars = c("mean1","sd1"))
bayesplot::mcmc_pairs(mod1, regex_pars = "^b", pars = c("mean1","sd1"))
bayesplot::mcmc_neff(neff_ratio(mod1, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)


# posterior predictive check
ppc_dens_overlay(dt,
                 rstan::extract(mod1, pars = 'y_hat')$y_hat[1:200,])









### Missing data model; completely at random (i.e., Bayesian imputation) ###

set.seed(2022)

covars <- data.frame(bathym = bathym,
                     chla = chla,
                     sst = sst)

rowNA <- sample(1:N, 0.25*N, replace = FALSE)  # determine which rows have missing values (25% of data)
colNA <- sample(1:ncol(covars), 0.25*N, replace = TRUE)  # determine which covars have missing values

for (i in 1:length(rowNA)) {
  covars[rowNA[i],colNA[i]] <- NA
}

head(covars, n=50)

bathym_missidx <- which(is.na(covars$bathym))
chla_missidx <- which(is.na(covars$chla))
sst_missidx <- which(is.na(covars$sst))


dat.list <- list(
  N = N,
  ID = ID,
  nID = nID,
  dt = dt,
  dist = dist,
  bathym = covars$bathym,
  chla = covars$chla,
  sst = covars$sst,
  bathym_missidx = bathym_missidx,
  n_bathym_miss = length(bathym_missidx),
  chla_missidx = chla_missidx,
  n_chla_miss = length(chla_missidx),
  sst_missidx = sst_missidx,
  n_sst_miss = length(sst_missidx)
)

stan.model <- '
functions {
  vector merge_missing(array[] int miss_indexes, vector x_obs, vector x_miss) {
    int N_obs = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];
    vector[N_obs] merged;
    merged = x_obs;
    for (i in 1:N_miss){
      merged[miss_indexes[i]] = x_miss[i];
      return merged;
    }
  }
}


data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst;                          // Sea surface temperature (C)
  int n_bathym_miss;
  int n_chla_miss;
  int n_sst_miss;
  int bathym_missidx[n_bathym_miss];
  int chla_missidx[n_chla_miss];
  int sst_missidx[n_sst_miss];
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
    vector[nID] b0_id;
    real b1;
    real b2;
    real b3;

    real mean1;
    real<lower=0> sd1;

    real nu1;
    real nu2;
    real nu3;
    real<lower=0> sigma1;
    real<lower=0> sigma2;
    real<lower=0> sigma3;
}


transformed parameters {
  vector[N] mu;
  real<lower=0> a[N];
  real bathym_merge[N];
  real chla_merge[N];
  real sst_merge[N];
  vector[n_bathym_miss] bathym_impute;
  vector[n_chla_miss] chla_impute;
  vector[n_sst_miss] sst_impute;

  bathym_merge = merge_missing(bathym_missidx, to_vector(bathym_s), bathym_impute);
  chla_merge = merge_missing(chla_missidx, to_vector(chla_s), chla_impute);
  sst_merge = merge_missing(sst_missidx, to_vector(sst_s), sst_impute);

  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_merge[i] + b2*chla_merge[i] + b3*sst_merge[i]);
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

  b1 ~ normal(0,1);
  b2 ~ normal(0,1);
  b3 ~ normal(0,1);

  b ~ normal(0,1);

  mean1 ~ normal(0,1);
  sd1 ~ normal(0,1);

  bathym_merge ~ normal(nu1,sigma1);
  chla_merge ~ normal(nu2,sigma2);
  sst_merge ~ normal(nu3,sigma3);

  nu1 ~ normal(0,0.5);
  nu2 ~ normal(0,0.5);
  nu3 ~ normal(0,0.5);

  sigma1 ~ exponential(1);
  sigma2 ~ exponential(1);
  sigma3 ~ exponential(1);
}


generated quantities {
  vector[N] y_hat;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_merged[i] + b2*chla_merged[i] + b3*sst_merged[i]);
  }
}
'



mod2 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
             seed = 8675309)

params <- c('b','b0_id','b1','b2','b3','mean1','sd1')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
