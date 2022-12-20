
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
# library(arrow)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# dat <- vroom("Processed_data/Input for time model.csv", delim = ",")
dat <- read_csv("Processed_data/Input for time model.csv")

glimpse(dat)
summary(dat)


# Change time step from secs to mins
dat$dt <- dat$dt/60

# # Retain only observed steps
# dat <- dat %>%
#   filter(obs == 1)

# Remove all observations w/ missing values
dat2 <- dat %>%
  drop_na(bathym, k490, npp, sst)
# dat2 <- dat  #all points moved off land; NAs represent missing depths or interpolation close to land


# Center and scale covariates
dat2 <- dat2 %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         k490.s = scale(k490) %>%
           as.vector(),
         npp.s = scale(npp) %>%
           as.vector(),
         sst.s = scale(sst) %>%
           as.vector())



# Check correlation among covars
cor(dat2[,c('bathym','k490','npp','sst')]) #all corrs < 0.6





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

# dat3 <- dat2 %>%
#   drop_na(bathym, Kd490, NPP, SST)
  # slice_sample(n = 0.01*nrow(dat2)) %>%
  # arrange(id, date)

# Do this step if not running the prior predictive simulation
dat2$id <- factor(dat2$id) %>%  #need to convert to integer for Stan
  as.numeric()

table(dat2$id) #47 IDs remaining


# Define objects pertaining to imputation of missing covars
# bathym_missidx <- which(is.na(dat3$bathym))
# n_bathym_miss <- length(bathym_missidx)
# k490_missidx <- which(is.na(dat3$Kd490))
# n_k490_miss <- length(k490_missidx)
# npp_missidx <- which(is.na(dat3$NPP))
# n_npp_miss <- length(npp_missidx)
# sst_missidx <- which(is.na(dat3$SST))
# n_sst_miss <- length(sst_missidx)


# data
dat.list <- list(
  N = nrow(dat2),
  ID = dat2$id,
  nID = n_distinct(dat2$id),
  dt = dat2$dt,
  dist = dat2$dist,
  bathym = dat2$bathym,
  k490 = dat2$k490,
  npp = dat2$npp,
  sst = dat2$sst#,
  # n_bathym_miss = n_bathym_miss,
  # bathym_missidx = bathym_missidx,
  # n_k490_miss = n_k490_miss,
  # k490_missidx = k490_missidx,
  # n_npp_miss = n_npp_miss,
  # npp_missidx = npp_missidx,
  # n_sst_miss = n_sst_miss,
  # sst_missidx = sst_missidx
  )

# Replace missing covar values w/ Inf
# dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
# dat.list$k490 <- ifelse(is.na(dat.list$k490), Inf, dat.list$k490)
# dat.list$npp <- ifelse(is.na(dat.list$npp), Inf, dat.list$npp)
# dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)


stan.model <- '
data {
  int N;                                  // sample size
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
  //int n_bathym_miss;                      // number of missing bathymetry vals
  //int bathym_missidx[n_bathym_miss];      // index for missing bathymetry vals
  //int n_k490_miss;                        // number of missing K490 vals
  //int k490_missidx[n_k490_miss];          // index for missing K490 vals
  //int n_npp_miss;                         // number of missing NPP vals
  //int npp_missidx[n_npp_miss];            // index for missing NPP vals
  //int n_sst_miss;                         // number of missing SST vals
  //int sst_missidx[n_sst_miss];            // index for missing SST vals
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

    //vector[n_bathym_miss] bathym_miss;
    //vector[n_k490_miss] k490_miss;
    //vector[n_npp_miss] npp_miss;
    //vector[n_sst_miss] sst_miss;

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
  //vector[N] bathym_merge;
  //vector[N] k490_merge;
  //vector[N] npp_merge;
  //vector[N] sst_merge;
  //vector[N] nu_k490;
  //vector[N] nu_npp;

  //bathym_merge = bathym;
  //k490_merge = k490;
  //npp_merge = npp;
  //sst_merge = sst;

  //bathym_merge[bathym_missidx] = bathym_miss;
  //k490_merge[k490_missidx] = k490_miss;
  //npp_merge[npp_missidx] = npp_miss;
  //sst_merge[sst_missidx] = sst_miss;


  // priors
  b0_id ~ normal(b0_bar, sd_b0);

  b0_bar ~ normal(0,1);
  sd_b0 ~ normal(0,1);

  mu_b ~ normal(0, 0.5);
  b ~ normal(mu_b, 0.5);

  [bBathym, bK490, bNPP, bSST] ~ normal(0, 1);

  //bathym_merge ~ normal(0, 1);
  //k490_merge ~ normal(nu_k490, sigma_k490);
  //npp_merge ~ normal(nu_npp, sigma_npp);
  //k490_merge ~ normal(0, 1);
  //npp_merge ~ normal(0, 1);
  //sst_merge ~ normal(0, 1);

  //nu_k490 = a_k490 + bBk490*bathym_merge + bBk490_2*square(bathym_merge);  // mean of Kd490 is quadratic function of bathym
  //nu_npp = a_npp + bKnpp*k490_merge + bBnpp*bathym_merge;  //mean of NPP is function of bathym and Kd490
  //sigma_k490 ~ normal(0,1);
  //sigma_npp ~ normal(0,1);

  //[bBk490, bBk490_2, bKnpp, bBnpp] ~ normal(0, 0.5);
  //[a_k490, a_npp] ~ normal(0, 1);



  for (i in 1:N) {
    // mean of gamma distribution
    mu[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym[i] + bK490*k490[i] + bNPP*npp[i] + bSST*sst[i]);

    // calculate the corresponding a parameter
    a[i] = mu[i] * b[ID[i]];

    // likelihood
    dt[i] ~ gamma(a[i], b[ID[i]]);
  }
}



generated quantities {
  vector[N] y_hat;

  // mean of gamma distribution
  for (i in 1:N){
    y_hat[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym[i] + bK490*k490[i] + bNPP*npp[i] + bSST*sst[i]);
  }
}
'

mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 8675309)
# took 15 min to run 2000 iter for 100% of remaining data

params <- c('mu_b','b','b0_id','bBathym','bK490','bNPP','bSST','b0_bar','sd_b0')
print(mod1, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
# print(mod1, digits_summary = 3, pars = 'bathym_miss', probs = c(0.025, 0.5, 0.975))
# print(mod1, digits_summary = 3, pars = 'k490_miss', probs = c(0.025, 0.5, 0.975))
# print(mod1, digits_summary = 3, pars = 'npp_miss', probs = c(0.025, 0.5, 0.975))
# print(mod1, digits_summary = 3, pars = 'sst_miss', probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod1, ind = TRUE, iter = 2000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod1, params = params)

bayesplot::mcmc_neff(neff_ratio(mod1, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)

## Save stanfit object
saveRDS(mod1, "Data_products/Time_model_intercept_stanfit.rds")


#posterior prediction
y_hat <- matrix(NA, nrow = 100, ncol = nrow(dat3))

bathym <- dat.list$bathym
k490 <- dat.list$k490
npp <- dat.list$npp
sst <- dat.list$sst
# bathym_merge[bathym_missidx] = extract(mod1, par = 'bathym_miss')$bathym_miss %>% colMeans()
# chla_merge = dat.list$chla
# chla_merge[chla_missidx] = extract(mod1, par = 'chla_miss')$chla_miss %>% colMeans()
# sst_merge = dat.list$sst
# sst_merge[sst_missidx] = extract(mod1, par = 'sst_miss')$sst_miss %>% colMeans()
b0_hat <- extract(mod1, par = 'b0_id')$b0_id
bBathym_hat <- extract(mod1, par = 'bBathym')$bBathym
bK490_hat <- extract(mod1, par = 'bK490')$bK490
bNPP_hat <- extract(mod1, par = 'bNPP')$bNPP
bSST_hat <- extract(mod1, par = 'bSST')$bSST
b_hat <- extract(mod1, par = 'b')$b

for (i in 1:nrow(y_hat)) {
  print(i)
  mu <- dat3$dist * exp(b0_hat[i,dat3$id] + bBathym_hat[i]*bathym + bK490_hat[i]*k490 + bNPP_hat[i]*npp + bSST_hat[i]*sst)
  y_hat[i,] <- mu
}

ppc_dens_overlay(dat.list$dt, y_hat) +
  xlim(0,1000)









#################################################
### Run model w/ varying intercept and slopes ###
#################################################

# data
dat.list <- list(
  N = nrow(dat2),
  K = 5,  #intercept and 4 environ covariates
  ID = dat2$id,
  nID = n_distinct(dat2$id),
  dt = dat2$dt,
  dist = dat2$dist,
  bathym = dat2$bathym.s,
  k490 = dat2$k490.s,
  npp = dat2$npp.s,
  sst = dat2$sst.s#,
  # n_bathym_miss = n_bathym_miss,
  # bathym_missidx = bathym_missidx,
  # n_k490_miss = n_k490_miss,
  # k490_missidx = k490_missidx,
  # n_npp_miss = n_npp_miss,
  # npp_missidx = npp_missidx,
  # n_sst_miss = n_sst_miss,
  # sst_missidx = sst_missidx
)

# Replace missing covar values w/ Inf
# dat.list$bathym <- ifelse(is.na(dat.list$bathym), Inf, dat.list$bathym)
# dat.list$k490 <- ifelse(is.na(dat.list$k490), Inf, dat.list$k490)
# dat.list$npp <- ifelse(is.na(dat.list$npp), Inf, dat.list$npp)
# dat.list$sst <- ifelse(is.na(dat.list$sst), Inf, dat.list$sst)


stan.model <- '
data {
  int N;                                  // sample size
  int K;                                  // number of terms in linear model
  int ID[N];                              // ID label for each step
  int nID;                                // number of unique IDs
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
    //vector[nID] b0_id;
    //real bBathym;
    //real bK490;
    //real bNPP;
    //real bSST;

    //real b0_bar;
    //real<lower=0> sd_b0;

  vector[K] mu_betas;                    // Population means of coefficients
  vector<lower=0>[K] tau;                // Population scales (SDs)
  cholesky_factor_corr[K] L;             // Population correlations
  vector[K] betas[nID];                  // Coeffs per individual (int and slopes)
}


model {
  vector[N] mu;
  vector[N] a;

  // priors
  //b0_id ~ normal(b0_bar, sd_b0);

  //b0_bar ~ normal(0,1);
  //sd_b0 ~ normal(0,1);

  mu_b ~ normal(0, 0.5);
  b ~ normal(mu_b, 0.5);

  //[bBathym, bK490, bNPP, bSST] ~ normal(0, 1);

  mu_betas ~ normal(0, 1);
  tau ~ normal(0, 1);
  L ~ lkj_corr_cholesky(4);

  // Estimate means per ID
  betas ~ multi_normal_cholesky(mu_betas, diag_pre_multiply(tau, L));


  for (i in 1:N) {
    // mean of gamma distribution
    mu[i] = dist[i] * exp(betas[ID[i],1] + betas[ID[i],2]*bathym[i] + betas[ID[i],3]*k490[i] + betas[ID[i],4]*npp[i] + betas[ID[i],5]*sst[i]);

    // calculate the corresponding a parameter
    a[i] = mu[i] * b[ID[i]];

    // likelihood
    dt[i] ~ gamma(a[i], b[ID[i]]);
  }
}



//generated quantities {
  //vector[N] y_hat;

  // mean of gamma distribution
  //for (i in 1:N){
    //y_hat[i] = dist[i] * exp(b0_id[ID[i]] + bBathym*bathym[i] + bK490*k490[i] + bNPP*npp[i] + bSST*sst[i]);
  //}
//}
'


mod2 <- rstan::stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000,
                    seed = 8675309, refresh = 100, control = list(max_treedepth = 15))
#took 8.9 hrs to run 2000 iter for full dataset

params <- c('mu_b','b','mu_betas','tau','L')
print(mod2, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
print(mod2, digits_summary = 3, pars = 'betas', probs = c(0.025, 0.5, 0.975))
# print(mod2, digits_summary = 3, pars = 'chla_miss', probs = c(0.025, 0.5, 0.975))
# print(mod2, digits_summary = 3, pars = 'sst_miss', probs = c(0.025, 0.5, 0.975))

MCMCtrace(mod2, ind = TRUE, iter = 1000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod2, params = params)
# pairs(mod2, pars = c('L_Rho','sigma'))

bayesplot::mcmc_neff(neff_ratio(mod2, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod2, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_neff(neff_ratio(mod2, pars = 'betas')) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod2, pars = 'betas')) +
  bayesplot::yaxis_text(hjust = 0)


## Save stanfit object
saveRDS(mod2, "Data_products/Time_model_intercept-slopes_stanfit.rds")




### NEED TO CREATE LOG-LIK PARAMS IN GENERATED QUANTITIES BLOCKS IF WANTING USE USE INFORMATION CRITERIA (I.E., WAIC, LOO) TO PERFORM MODEL SELECTION
