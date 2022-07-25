

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








#####################################################################
### Random intercept model w/ missing data (completely at random) ###
#####################################################################

# This modeling approach incorporates full Bayesian imputation

set.seed(2022)

covars <- data.frame(bathym = bathym,
                     chla = chla,
                     sst = sst)

rowNA <- sample(1:N, 0.25*N, replace = FALSE)  # determine which rows have missing values (25% of data)
colNA <- sample(1:ncol(covars), 0.25*N, replace = TRUE)  # determine which covars have missing values

for (i in 1:length(rowNA)) {
  covars[rowNA[i], colNA[i]] <- NA
}

head(covars, n=50)

chla_missidx <- which(is.na(covars$chla))
n_chla_miss <- length(chla_missidx)
bathym_missidx <- which(is.na(covars$bathym))
n_bathym_miss <- length(bathym_missidx)
sst_missidx <- which(is.na(covars$sst))
n_sst_miss <- length(sst_missidx)



dat.list <- list(
  N = N,
  ID = ID,
  nID = nID,
  dt = dt,
  dist = dist,
  bathym = scale(covars$bathym)[,1],
  chla = scale(covars$chla)[,1],
  sst = scale(covars$sst)[,1],
  n_chla_miss = n_chla_miss,
  chla_missidx = chla_missidx,
  n_bathym_miss = n_bathym_miss,
  bathym_missidx = bathym_missidx,
  n_sst_miss = n_sst_miss,
  sst_missidx = sst_missidx
)
dat.list$chla <- ifelse(is.na(dat.list$chla), Inf, dat.list$chla)
dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)

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

}


parameters {
    real<lower=0> b;
    vector[nID] b0_id;
    real b1;
    real b2;
    real b3;

    real mean1;
    real<lower=0> sd1;

    vector[n_chla_miss] chla_miss;
    vector[n_bathym_miss] bathym_miss;
    vector[n_sst_miss] sst_miss;

    real nu_chla;
    real nu_bathym;
    real nu_sst;
    real<lower=0> sigma_chla;
    real<lower=0> sigma_bathym;
    real<lower=0> sigma_sst;
}


transformed parameters {
  vector[N] mu;
  real<lower=0> a[N];
  vector[N] chla_merge;
  vector[N] bathym_merge;
  vector[N] sst_merge;

  // merge observed and imputed Chla vals
  chla_merge = chla;
    for (i in 1:n_chla_miss) {
        chla_merge[chla_missidx[i]] = chla_miss[i];
    }

  // merge observed and imputed bathymetry vals
  bathym_merge = bathym;
    for (i in 1:n_bathym_miss) {
        bathym_merge[bathym_missidx[i]] = bathym_miss[i];
    }

  // merge observed and imputed SST vals
  sst_merge = sst;
    for (i in 1:n_sst_miss) {
        sst_merge[sst_missidx[i]] = sst_miss[i];
    }

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

  mean1 ~ normal(0,1);
  sd1 ~ normal(0,1);

  [sigma_chla, sigma_bathym, sigma_sst] ~ normal(0,2);
  [nu_chla, nu_bathym, nu_sst, b, b1, b2, b3] ~ normal(0, 1);
  chla_merge ~ normal(nu_chla, sigma_chla);
  bathym_merge ~ normal(nu_bathym, sigma_bathym);
  sst_merge ~ normal(nu_sst, sigma_sst);

}


generated quantities {
  vector[N] y_hat;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(b0_id[ID[i]] + b1*bathym_merge[i] + b2*chla_merge[i] + b3*sst_merge[i]);
  }
}
'



mod2 <- rstan::stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
             seed = 8675309)

params <- c('b','b0_id','b1','b2','b3','mean1','sd1','nu_chla','sigma_chla','nu_bathym',
            'sigma_bathym','nu_sst','sigma_sst')
print(mod2, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod2, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod2, params = params)

# posterior predictive check
ppc_dens_overlay(dt,
                 rstan::extract(mod2, pars = 'y_hat')$y_hat[1:200,])












#####################################################################
### Random intercept model w/ missing data (completely at random) ###
#####################################################################

set.seed(2022)

# Define coefficients
mu_b0 <- -0.4
mu_b1 <- -0.5
mu_b2 <- 0
mu_b3 <- 0.5
Mu <- c(mu_b0, mu_b1, mu_b2, mu_b3)

sigma_b0 <- 1
sigma_b1 <- 0.6
sigma_b2 <- 1.3
sigma_b3 <- 1.7
sigma <- c(sigma_b0, sigma_b1, sigma_b2, sigma_b3)

rho <- 0.65
Rho <- matrix(rho, nrow = 4, ncol = 4)
diag(Rho) <- 1

Sigma <- diag(sigma) %*% Rho %*% diag(sigma)

betas <- mvtnorm::rmvnorm(nID, Mu, Sigma, method = "chol")

b <- 0.05

# Generate time interval (dt) using subset of real data for `dist` and environ covars
mu <- dist * exp(betas[ID,1] + betas[ID,2]*bathym.s[,1] + betas[ID,3]*chla.s[,1] +
                   betas[ID,4]*sst.s[,1])
a <- mu * b
dt <- rgamma(N, shape = a, rate = b)



dat.list <- list(
  N = N,
  ID = ID,
  nID = nID,
  dt = dt,
  dist = dist,
  bathym = scale(covars$bathym)[,1],
  chla = scale(covars$chla)[,1],
  sst = scale(covars$sst)[,1],
  n_chla_miss = n_chla_miss,
  chla_missidx = chla_missidx,
  n_bathym_miss = n_bathym_miss,
  bathym_missidx = bathym_missidx,
  n_sst_miss = n_sst_miss,
  sst_missidx = sst_missidx
)
dat.list$chla <- ifelse(is.na(dat.list$chla), Inf, dat.list$chla)
dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)

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

}


parameters {
    real<lower=0> b;
    matrix[nID,4] betas;
    vector[4] mu_betas;
    cholesky_factor_corr[4] L_Rho;
    vector<lower=0>[4] sigma;

    vector[n_chla_miss] chla_miss;
    vector[n_bathym_miss] bathym_miss;
    vector[n_sst_miss] sst_miss;

    real nu_chla;
    real nu_bathym;
    real nu_sst;
    real<lower=0> sigma_chla;
    real<lower=0> sigma_bathym;
    real<lower=0> sigma_sst;
}


transformed parameters {
  vector[N] mu;
  real<lower=0> a[N];
  vector[N] chla_merge;
  vector[N] bathym_merge;
  vector[N] sst_merge;
  matrix[4,4] L_Sigma;

  // generate coefficients for regression from MVN w/ Cholesky decomp
  L_Sigma = diag_pre_multiply(sigma, L_Rho);

  // merge observed and imputed Chla vals
  chla_merge = chla;
    for (i in 1:n_chla_miss) {
        chla_merge[chla_missidx[i]] = chla_miss[i];
    }

  // merge observed and imputed bathymetry vals
  bathym_merge = bathym;
    for (i in 1:n_bathym_miss) {
        bathym_merge[bathym_missidx[i]] = bathym_miss[i];
    }

  // merge observed and imputed SST vals
  sst_merge = sst;
    for (i in 1:n_sst_miss) {
        sst_merge[sst_missidx[i]] = sst_miss[i];
    }

  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(betas[ID[i],1] + betas[ID[i],2]*bathym_merge[i] + betas[ID[i],3]*chla_merge[i] + betas[ID[i],4]*sst_merge[i]);

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
  to_vector(mu_betas) ~ normal(0,1);
  sigma ~ normal(0,1);
  L_Rho ~ lkj_corr_cholesky(2);

  for (j in 1:nID) {
    betas[j,] ~ multi_normal_cholesky(mu_betas, L_Sigma);
  }


  [sigma_chla, sigma_bathym, sigma_sst] ~ normal(0,2);
  [nu_chla, nu_bathym, nu_sst, b] ~ normal(0, 1);
  chla_merge ~ normal(nu_chla, sigma_chla);
  bathym_merge ~ normal(nu_bathym, sigma_bathym);
  sst_merge ~ normal(nu_sst, sigma_sst);

}


generated quantities {
  vector[N] y_hat;
  matrix[4,4] Rho;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(betas[ID[i],1] + betas[ID[i],2]*bathym_merge[i] + betas[ID[i],3]*chla_merge[i] + betas[ID[i],4]*sst_merge[i]);
  }

  // Correlations from variance-covariance matrix
  Rho = multiply_lower_tri_self_transpose(L_Rho);
}
'


mod3 <- rstan::stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
                    seed = 8675309)

params <- c('b','betas','mu_betas','sigma','Rho','nu_chla','sigma_chla','nu_bathym',
            'sigma_bathym','nu_sst','sigma_sst')
print(mod3, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod3, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod3, params = params)

# posterior predictive check
ppc_dens_overlay(dt,
                 rstan::extract(mod3, pars = 'y_hat')$y_hat[1:200,]) +
  xlim(0,1000)

