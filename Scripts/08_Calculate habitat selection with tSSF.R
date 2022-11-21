
### Run time-explicit step-selection function ###

library(tidyverse)
library(lubridate)
library(rstan)
library(MCMCvis)
library(bayesplot)
library(arrow)
library(sf)
library(furrr)
library(future)
library(tictoc)
library(progressr)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


### Load model fit, model input, and tracks ###

mod <- readRDS('Data_products/Time_model_intercept_stanfit.rds')
mod.input <- read_parquet("Processed_data/Input for time model.parquet")
dat <- read_csv('Processed_data/Imputed_GoM_Cm_Tracks_SSM_2hr.csv')




### Wrangle model input to required format ###

# Change time step from secs to mins
mod.input$dt <- mod.input$dt/60

# Remove all observations w/ missing bathym values (since NA values were assigned to land)
mod.input2 <- mod.input %>%
  drop_na(bathym, Kd490, NPP, SST) %>%
  mutate(id1 = as.numeric(factor(id)), .after = id)

# Add quadratic terms for covariates
mod.input2 <- mod.input2 %>%
  mutate(bathym.2 = bathym^2,
         Kd490.2 = Kd490^2,
         NPP.2 = NPP^2,
         SST.2 = SST^2)

# Center and scale covariates
mod.input2 <- mod.input2 %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         kd490.s = scale(Kd490) %>%
           as.vector(),
         npp.s = scale(NPP) %>%
           as.vector(),
         sst.s = scale(SST) %>%
           as.vector(),
         bathym.2.s = scale(bathym.2) %>%
           as.vector(),
         kd490.2.s = scale(Kd490.2) %>%
           as.vector(),
         npp.2.s = scale(NPP.2) %>%
           as.vector(),
         sst.2.s = scale(SST.2) %>%
           as.vector())


# Only keep strata that have at least 1 used and 1 available step (since only complete cases analyzed)
mod.input3 <- mod.input2 %>%
  group_by(strata) %>%
  filter(sum(obs) == 1) %>%
  filter(n() > 1)



### Calculate probabilities of used and available steps based on time model results ###
calc_time_probs <- function(mod, dat, covar.names, p) {
  # mod = the rstan model fit for the time model
  # dat = the input for the time model that includes cols for id, dist, dt, strata, covars, and obs
  # covar names = a vector of names (in proper order) for the covars on which to make predictions
  # p = a stored 'progressr' object for creating progress bar

  # message("Calculating probabilities from time model for PTT ", names(dat$id)[1], "...")

  dat <- data.frame(dat)

  #calculate probabilities
  betas <- rstan::extract(mod, pars = c('b0_id','bBathym','bK490','bNPP','bSST')) %>%
    bind_cols()
  b <- rstan::extract(mod, pars = 'b')$b
  strat <- unique(dat$strata)
  prob.all <- vector("list", length(strat))

  for (i in 1:length(strat)) {
    # print(i)

    #get parameters
    id1 <- dat[dat$strata == strat[i], "id1"][1]
    dist1 <- dat[dat$strata == strat[i], "dist"][1]
    betas1 <- cbind(betas$b0_id[,id1], betas[,-1])
    b1 <- b[,id1]

    #define design vector
    xmat <- as.matrix(cbind(1, dat[dat$strata == strat[i], covar.names]))
    #calculate mean using linear algebra
    mean1 <- dist1 * exp(xmat %*% t(betas1))
    #calculate a and b parameters of gamma distribution
    a1 <- b1 * t(mean1)

    #calculate probabilities for the realized and potential steps
    prob <- rep(NA, ncol(a1))
    for (j in 1:ncol(a1)) {
      #calculate posterior median of gamma density
      prob[j] <- median(dgamma(dat[dat$strata == strat[i], "dt"][j], a1[,j], b1))
    }

    prob.all[[i]] <- prob
  }

  p()  #plot progress bar

  prob.all <- unlist(prob.all)
  return(prob.all)
}


mod.input.list <- split(mod.input3, mod.input3$id)

plan(multisession, workers = availableCores() - 2)

tic()
progressr::with_progress({
  #set up progress bar
  p<- progressr::progressor(steps = length(mod.input.list))

res <- future_map(mod.input.list,
                  ~calc_time_probs(mod = mod, dat = ., covar.names = c('bathym.s','kd490.s','npp.s','sst.s'), p = p),
                  .options = furrr_options(seed = 2022))
})
toc()  #took 73 min to run

plan(sequential)


mod.input3$time.prob <- unlist(res)
# mod.input3$time.prob <- ifelse(mod.input3$time.prob == 0, 1e-99, mod.input3$time.prob)



########################################################
### Fit time-explicit step-selection function (tSSF) ###
########################################################

# Prepare data for model
mod.input3 <- mod.input3 %>%
  group_by(strata) %>%
  mutate(step.id = 1:n(), .after = obs) %>%  #number each of the used and available steps for subsetting
  ungroup()

# Replace values for missing steps per strata as 0s
plan(multisession, workers = availableCores() - 2)

tic()
mod.input4 <- mod.input3 %>%
  mutate(step.id = factor(step.id)) %>%
  split(.$id1) %>%
  future_map(., ~{.x %>%
      group_by(strata) %>%
      complete(step.id) %>%
      mutate(across(obs:time.prob, replace_na, 0)) %>%
      ungroup()
    }) %>%
  bind_rows()
toc()  #took 16.5 min

plan(sequential)


# Create list of data
mod.input5 <- mod.input4 #%>%
  # slice(693391:693400) #%>%
  # slice(510301:510400)
  # slice(1:1000000)
  # filter(strata %in% sample(unique(mod.input4$strata), size = 10000, replace = FALSE))

# Remove obs w/ time gap > 12 hours
ind <- mod.input5 %>%
  filter(obs == 1) %>%
  filter(dt > 12*60) %>%
  pull(strata)

mod.input5 <- mod.input5 %>%
  filter(!strata %in% ind)  #remove any steps w/ time gaps > 12 hrs

#get names of covariates
covar.names <- c('bathym.s','bathym.2.s','kd490.s','kd490.2.s','npp.s','npp.2.s','sst.s','sst.2.s')


dat.list <- list(
  N = n_distinct(mod.input5$strata),
  K = length(covar.names),
  ID = mod.input5 %>%
    distinct(strata, .keep_all = TRUE) %>%
    pull(id1),
  nID = n_distinct(mod.input5$id, na.rm = TRUE),
  xmat_used = as.matrix(mod.input5[mod.input5$step.id == 1, covar.names]),
  xmat_avail1 = as.matrix(mod.input5[mod.input5$step.id == 2, covar.names]),
  xmat_avail2 = as.matrix(mod.input5[mod.input5$step.id == 3, covar.names]),
  xmat_avail3 = as.matrix(mod.input5[mod.input5$step.id == 4, covar.names]),
  xmat_avail4 = as.matrix(mod.input5[mod.input5$step.id == 5, covar.names]),
  pmov_used = mod.input5[mod.input5$step.id == 1,]$time.prob,
  pmov_avail1 = mod.input5[mod.input5$step.id == 2,]$time.prob,
  pmov_avail2 = mod.input5[mod.input5$step.id == 3,]$time.prob,
  pmov_avail3 = mod.input5[mod.input5$step.id == 4,]$time.prob,
  pmov_avail4 = mod.input5[mod.input5$step.id == 5,]$time.prob,
  y = rep(1, n_distinct(mod.input5$strata))
)





# Run model
stan.model <- '
data {
  int N;                                 // number of unique strata
  int K;                                 // number of covariates
  int nID;                               // number of unique IDs
  int<lower=1, upper=nID> ID[N];         // ID label for each step
  matrix[N,K] xmat_used;                 // Design matrix for used steps
  matrix[N,K] xmat_avail1;               // Design matrix for first available step
  matrix[N,K] xmat_avail2;               // Design matrix for second available step
  matrix[N,K] xmat_avail3;               // Design matrix for third available step
  matrix[N,K] xmat_avail4;               // Design matrix for fourth available step

  vector<lower=0>[N] pmov_used;          // Time model probs for used steps
  vector<lower=0>[N] pmov_avail1;        // Time model probs for first available step
  vector<lower=0>[N] pmov_avail2;        // Time model probs for second available step
  vector<lower=0>[N] pmov_avail3;        // Time model probs for third available step
  vector<lower=0>[N] pmov_avail4;        // Time model probs for fourth available step

  int<lower=1, upper=1> y[N];            // Ones for "ones trick" of Bernoulli model
}


parameters {
  //vector[K] betas;
  vector[K] mu;                          // Population means of slopes
  vector<lower=0>[K] tau;                // Population scales (SDs)
  cholesky_factor_corr[K] L;             // Population correlations
  vector[K] betas[nID];                  // Slopes per individual

}


model {
    vector[N] p_used;                     // Vector to store probs for used steps
    vector[N] p_avail1;                   // Vector to store probs for first available step
    vector[N] p_avail2;                   // Vector to store probs for second available step
    vector[N] p_avail3;                   // Vector to store probs for third available step
    vector[N] p_avail4;                   // Vector to store probs for fourth available step
    vector[N] pi;                         // Vector to store tSSF probs

// priors
  //betas ~ normal(0, 1);
  mu ~ normal(0, 1);
  tau ~ normal(0, 2);
  L ~ lkj_corr_cholesky(K);

  // Estimate means per ID
  betas ~ multi_normal_cholesky(mu, diag_pre_multiply(tau, L));



// model
  for (i in 1:N) {
    // calculation of pi
    p_used[i] = pmov_used[i] * exp(dot_product(xmat_used[i,], betas[ID[i]]));
    p_avail1[i] = pmov_avail1[i] * exp(dot_product(xmat_avail1[i,], betas[ID[i]]));
    p_avail2[i] = pmov_avail2[i] * exp(dot_product(xmat_avail2[i,], betas[ID[i]]));
    p_avail3[i] = pmov_avail3[i] * exp(dot_product(xmat_avail3[i,], betas[ID[i]]));
    p_avail4[i] = pmov_avail4[i] * exp(dot_product(xmat_avail4[i,], betas[ID[i]]));

    pi[i] = p_used[i] / (p_used[i] + p_avail1[i] + p_avail2[i] + p_avail3[i] + p_avail4[i]);


    // calculate likelihood using "ones trick"
    y[i] ~ bernoulli(pi[i]);
  }

  // calculation of pi
    //p_used = pmov_used .* exp(xmat_used * betas);
    //p_avail1 = pmov_avail1 .* exp(xmat_avail1 * betas);
    //p_avail2 = pmov_avail2 .* exp(xmat_avail2 * betas);
    //p_avail3 = pmov_avail3 .* exp(xmat_avail3 * betas);
    //p_avail4 = pmov_avail4 .* exp(xmat_avail4 * betas);

    //pi = p_used ./ (p_used + p_avail1 + p_avail2 + p_avail3 + p_avail4);


    // calculate likelihood using "ones trick"
    //target += bernoulli_lpmf(y | pi);
}

'


mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 8675309)
# took 4.8 hrs to run 2000 iter for full dataset

# params <- c('mu_b','b','b0_id','bBathym','bK490','bNPP','bSST','b0_bar','sd_b0')
print(mod1, digits_summary = 3, probs = c(0.025, 0.5, 0.975))


bayesplot::mcmc_neff(neff_ratio(mod1)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod1)) +
  bayesplot::yaxis_text(hjust = 0)
MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = c('mu','tau'))
par(mfrow=c(1,1))



## Save stanfit object
saveRDS(mod1, "Data_products/tSSF_model_GLM_stanfit.rds")
