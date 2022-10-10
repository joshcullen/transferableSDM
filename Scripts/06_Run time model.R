
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
dat2 <- dat  #all points moved off land; NAs represent missing depths or interpolation close to land


# Center and scale covariates
dat2 <- dat2 %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         kd490.s = scale(Kd490) %>%
           as.vector(),
         npp.s = scale(NPP) %>%
           as.vector(),
         sst.s = scale(SST) %>%
           as.vector())


# Check correlation among covars
cor(dat2[,c('bathym.s','kd490.s','npp.s','sst.s')] %>%
      drop_na())
dat2 %>%
  drop_na(bathym, Kd490, NPP, SST) %>%
  slice_sample(n = 50000) %>%
  dplyr::select(bathym.s, kd490.s, npp.s, sst.s) %>%
  pairs()
#strongest corr is between NPP and Kd490 (0.59); leaving all vars in for subsequent analyses



## Explore relationship between dt and each of the covars

ggplot(dat2 %>%
         slice_sample(n = 500000), aes(bathym, dt)) +
  geom_point() +
  theme_bw()

ggplot(dat2 %>%
         slice_sample(n = 500000), aes(Kd490, dt)) +
  geom_point() +
  theme_bw()

ggplot(dat2 %>%
         slice_sample(n = 500000), aes(NPP, dt)) +
  geom_point() +
  theme_bw()

ggplot(dat2 %>%
         slice_sample(n = 500000), aes(SST, dt, color = id)) +
  geom_point() +
  theme_bw()




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

  bBathym_p <- rnorm(1, 0, 1)
  bChla_p <- rnorm(1, 0, 1)
  bSST_p <- rnorm(1, 0, 1)

  b_p <- abs(rnorm(1, 0, 1))


  # mean of gamma distribution
  mu = dat2$dist * exp(b0[dat2$id] + bBathym_p*dat2$bathym.s + bChla_p*dat2$chla.s + bSST_p*dat2$sst.s)

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
# common across IDs; this is likely an artifact of setting the imputed time interval to 2 hrs (120 min)





## Run model

# Do this step if not running the prior predictive simulation
dat2$id <- factor(dat2$id) %>%  #need to convert to integer for Stan
  as.numeric()

table(dat2$id)


dat3 <- dat2 #%>%
  # slice_sample(n = 0.01*nrow(dat2)) %>%
  # arrange(id, date)

# Define objects pertaining to imputation of missing covars
bathym_missidx <- which(is.na(dat3$bathym))
n_bathym_miss <- length(bathym_missidx)
k490_missidx <- which(is.na(dat3$Kd490))
n_k490_miss <- length(k490_missidx)
npp_missidx <- which(is.na(dat3$NPP))
n_npp_miss <- length(npp_missidx)
sst_missidx <- which(is.na(dat3$SST))
n_sst_miss <- length(sst_missidx)


# data
dat.list <- list(
  N = nrow(dat3),
  ID = dat3$id,
  nID = n_distinct(dat3$id),
  dt = dat3$dt,
  dist = dat3$dist,
  bathym = scale(dat3$bathym)[,1],
  k490 = scale(dat3$Kd490)[,1],
  npp = scale(dat3$NPP)[,1],
  sst = scale(dat3$SST)[,1],
  n_bathym_miss = n_bathym_miss,
  bathym_missidx = bathym_missidx,
  n_k490_miss = n_k490_miss,
  k490_missidx = k490_missidx,
  n_npp_miss = n_npp_miss,
  npp_missidx = npp_missidx,
  n_sst_miss = n_sst_miss,
  sst_missidx = sst_missidx
  )

# Replace missing covar values w/ Inf
dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
dat.list$k490 <- ifelse(is.na(dat.list$k490), Inf, dat.list$k490)
dat.list$npp <- ifelse(is.na(dat.list$npp), Inf, dat.list$npp)
dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)


stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
  int n_bathym_miss;                      // number of missing bathymetry vals
  int bathym_missidx[n_bathym_miss];      // index for missing bathymetry vals
  int n_k490_miss;                        // number of missing K490 vals
  int k490_missidx[n_k490_miss];          // index for missing K490 vals
  int n_npp_miss;                         // number of missing NPP vals
  int npp_missidx[n_npp_miss];            // index for missing NPP vals
  int n_sst_miss;                         // number of missing SST vals
  int sst_missidx[n_sst_miss];            // index for missing SST vals
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] k490;                         // Light attenuation coefficient at 490 nm (units)
  vector[N] npp;                          // Net primary productivity (units)
  vector[N] sst;                          // Sea surface temperature (C)

}


parameters {
    real<lower=0> mu_b;
    vector<lower=0>[nID] b;
    vector[nID] b0_id;
    real bBathym;
    real bK490;
    real bNPP;
    real bSST;

    real b0_bar;
    real<lower=0> sd_b0;

    vector[n_bathym_miss] bathym_miss;
    vector[n_k490_miss] k490_miss;
    vector[n_npp_miss] npp_miss;
    vector[n_sst_miss] sst_miss;

    //real<lower=0> sigma_k490;
    //real<lower=0> sigma_npp;

    //real a_k490;
    //real a_npp;
    //real bBk490;
    //real bBk490_2;
    //real bKnpp;
    //real bBnpp;
}


model {
  vector[N] mu;
  vector[N] a;
  vector[N] bathym_merge;
  vector[N] k490_merge;
  vector[N] npp_merge;
  vector[N] sst_merge;
  //vector[N] nu_k490;
  //vector[N] nu_npp;

  bathym_merge = bathym;
  k490_merge = k490;
  npp_merge = npp;
  sst_merge = sst;

  bathym_merge[bathym_missidx] = bathym_miss;
  k490_merge[k490_missidx] = k490_miss;
  npp_merge[npp_missidx] = npp_miss;
  sst_merge[sst_missidx] = sst_miss;


  // priors
  b0_id ~ normal(b0_bar, sd_b0);

  b0_bar ~ normal(0,1);
  sd_b0 ~ normal(0,1);

  mu_b ~ normal(0, 0.5);
  b ~ normal(mu_b, 0.5);

  [bBathym, bK490, bNPP, bSST] ~ normal(0, 1);

  bathym_merge ~ normal(0, 1);
  //k490_merge ~ normal(nu_k490, sigma_k490);
  //npp_merge ~ normal(nu_npp, sigma_npp);
  k490_merge ~ normal(0, 1);
  npp_merge ~ normal(0, 1);
  sst_merge ~ normal(0, 1);

  //nu_k490 = a_k490 + bBk490*bathym_merge + bBk490_2*square(bathym_merge);  // mean of Kd490 is quadratic function of bathym
  //nu_npp = a_npp + bKnpp*k490_merge + bBnpp*bathym_merge;  //mean of NPP is function of bathym and Kd490
  //sigma_k490 ~ normal(0,1);
  //sigma_npp ~ normal(0,1);

  //[bBk490, bBk490_2, bKnpp, bBnpp] ~ normal(0, 0.5);
  //[a_k490, a_npp] ~ normal(0, 1);



  for (i in 1:N) {
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym_merge[i] + bK490*k490_merge[i] + bNPP*npp_merge[i] + bSST*sst_merge[i]);

    // calculate the corresponding a parameter
    a[i] = mu[i] * b[ID[i]];

    // likelihood
    dt[i] ~ gamma(a[i], b[ID[i]]);
  }
}
'

mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 3000, warmup = 1000, seed = 8675309)
# took 100 hrs to run 3000 iter for 100% of data

params <- c('mu_b','b','b0_id','bBathym','bK490','bNPP','bSST','b0_bar','sd_b0')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
print(mod1, digits_summary = 3, pars = 'bathym_miss', probs = c(0.025, 0.5, 0.975))
print(mod1, digits_summary = 3, pars = 'k490_miss', probs = c(0.025, 0.5, 0.975))
print(mod1, digits_summary = 3, pars = 'npp_miss', probs = c(0.025, 0.5, 0.975))
print(mod1, digits_summary = 3, pars = 'sst_miss', probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod1, ind = TRUE, iter = 4000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod1, params = params)

bayesplot::mcmc_neff(neff_ratio(mod1, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)

## Save stanfit object
saveRDS(mod1, "Data_products/Time_model_intercept_stanfit.rds")


#posterior prediction
y_hat <- matrix(NA, nrow = 100, ncol = nrow(dat3))

bathym = dat.list$bathym
# bathym_merge[bathym_missidx] = extract(mod1, par = 'bathym_miss')$bathym_miss %>% colMeans()
chla_merge = dat.list$chla
chla_merge[chla_missidx] = extract(mod1, par = 'chla_miss')$chla_miss %>% colMeans()
sst_merge = dat.list$sst
sst_merge[sst_missidx] = extract(mod1, par = 'sst_miss')$sst_miss %>% colMeans()
b0_hat <- extract(mod1, par = 'b0_id')$b0_id
bBathym_hat <- extract(mod1, par = 'bBathym')$bBathym
bChla_hat <- extract(mod1, par = 'bChla')$bChla
bSST_hat <- extract(mod1, par = 'bSST')$bSST
b_hat <- extract(mod1, par = 'b')$b

for (i in 1:nrow(y_hat)) {
  print(i)
  mu <- dat3$dist * exp(b0_hat[i,dat3$id] + bBathym_hat[i]*bathym + bChla_hat[i]*chla_merge + bSST_hat[i]*sst_merge)
  y_hat[i,] <- mu
}

ppc_dens_overlay(dat.list$dt, y_hat) +
  xlim(0,1000)









#################################################
### Run model w/ varying intercept and slopes ###
#################################################

# Do this step if not running the prior predictive simulation
dat2$id <- factor(dat2$id) %>%  #need to convert to integer for Stan
  as.numeric()

dat3 <- dat2 #%>%
  # slice_sample(n = 0.05*nrow(dat2)) %>%
  # arrange(id, date)

# Define objects pertaining to imputation of missing covars
# bathym_missidx <- which(is.na(dat3$bathym))
# n_bathym_miss <- length(bathym_missidx)
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
  bathym = dat3$bathym.s,
  chla = dat3$chla.s,
  sst = dat3$sst.s,
  # n_bathym_miss = n_bathym_miss,
  # bathym_missidx = bathym_missidx,
  n_chla_miss = n_chla_miss,
  chla_missidx = chla_missidx,
  n_sst_miss = n_sst_miss,
  sst_missidx = sst_missidx
)

# Replace missing covar values w/ Inf
# dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
dat.list$chla <- ifelse(is.na(dat.list$chla), Inf, dat.list$chla)
dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)



stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
  int n_chla_miss;                        // number of missing Chla vals
  int chla_missidx[n_chla_miss];          // index for missing Chla vals
  //int n_bathym_miss;                      // number of missing bathymetry vals
  //int bathym_missidx[n_bathym_miss];      // index for missing bathymetry vals
  int n_sst_miss;                         // number of missing SST vals
  int sst_missidx[n_sst_miss];            // index for missing SST vals
  vector[N] dt;                           // time interval (min)
  vector[N] dist;                         // Distance traveled for given step (m)
  vector[N] bathym;                       // Bathymetric depth (m)
  vector[N] chla;                         // Chlorophyll a concentration (mg m^-3 d^-1)
  vector[N] sst;                          // Sea surface temperature (C)

}


parameters {
    real<lower=0> mu_b;
    vector<lower=0>[nID] b;
    matrix[nID,4] betas;
    vector[4] mu_betas;
    corr_matrix[4] L_Rho;
    vector<lower=0>[4] sigma;

    vector[n_chla_miss] chla_miss;
    //vector[n_bathym_miss] bathym_miss;
    vector[n_sst_miss] sst_miss;

    //real nu_chla;
    //real nu_bathym;
    //real nu_sst;
    //real<lower=0> sigma_chla;
    //real<lower=0> sigma_bathym;
    //real<lower=0> sigma_sst;
}


model {
  vector[N] mu;
  real a[N];
  vector[N] chla_merge;
  //vector[N] bathym_merge;
  vector[N] sst_merge;

  chla_merge = chla;
  //bathym_merge = bathym;
  sst_merge = sst;

  chla_merge[chla_missidx] = chla_miss;
  //bathym_merge[bathym_missidx] = bathym_miss;
  sst_merge[sst_missidx] = sst_miss;




  // priors
  mu_betas ~ normal(0,0.5);
  sigma ~ normal(0,0.5);
  L_Rho ~ lkj_corr(2);

  for (j in 1:nID) {
    betas[j,] ~ multi_normal(mu_betas, quad_form_diag(L_Rho, sigma));
  }


  mu_b ~ normal(0, 0.5);
  b ~ normal(mu_b, 0.5);
  //[sigma_chla, sigma_bathym, sigma_sst] ~ normal(0,1);
  //[nu_chla, nu_bathym, nu_sst] ~ normal(0, 0.5);
  //target += normal_lpdf(bathym_merge | 0, 1);
  target += normal_lpdf(chla_merge | 0, 1);
  target += normal_lpdf(sst_merge | 0, 1);


  for (i in 1:N){
    // mean of gamma distribution
    mu[i] = dist[i] * exp(betas[ID[i],1] + betas[ID[i],2]*bathym[i] + betas[ID[i],3]*chla_merge[i] + betas[ID[i],4]*sst_merge[i]);

    // calculate the corresponding a[i] parameter
    a[i] = mu[i] * b[ID[i]];

    // likelihood
    dt[i] ~ gamma(a[i], b[ID[i]]);
  }
}
'


mod2 <- rstan::stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
                    seed = 8675309, refresh = 100)
#took 14 days to run 2000 iter for full dataset

params <- c('mu_b','b','mu_betas','sigma','L_Rho')
print(mod2, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
print(mod2, digits_summary = 3, pars = 'betas', probs = c(0.025, 0.5, 0.975))
print(mod2, digits_summary = 3, pars = 'chla_miss', probs = c(0.025, 0.5, 0.975))
print(mod2, digits_summary = 3, pars = 'sst_miss', probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod2, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod2, params = params)
# pairs(mod2, pars = c('L_Rho','sigma'))

bayesplot::mcmc_neff(neff_ratio(mod2, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod2, pars = 'betas')) +
  bayesplot::yaxis_text(hjust = 0)


## Save stanfit object
saveRDS(mod2, "Data_products/Time_model_intercept-slopes_stanfit.rds")
