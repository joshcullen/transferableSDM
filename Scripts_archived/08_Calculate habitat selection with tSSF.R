
### Run time-explicit step-selection function ###

library(tidyverse)
library(lubridate)
library(rstan)
library(MCMCvis)
library(bayesplot)
library(terra)
library(sf)
library(furrr)
library(future)
library(tictoc)
library(progressr)

source('Scripts/helper functions.R')


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


### Load model fit, model input, and tracks ###

mod <- readRDS('Data_products/Time_model_intercept-slopes_stanfit.rds')
# time.mod.input <- read_csv("Processed_data/Input for time model.csv")
tSSF.input <- read_csv("Processed_data/Input for tSSF.csv")
dat <- read_csv('Processed_data/Processed_GoM_Cm_Tracks_SSM_12hr_aniMotum.csv')

# glimpse(time.mod.input)
glimpse(tSSF.input)
# summary(time.mod.input)
summary(tSSF.input)


# Change time step from secs to hrs
tSSF.input$dt <- tSSF.input$dt/3600



#######################################
### Import Environmental Covariates ###
#######################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
names(cov_list) <- c('bathym', 'k490', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('k490', 'npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- resample(cov_list[[var]], cov_list$npp, method = "average")
}


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')




##############################
### Extract environ covars ###
##############################

# Extract covar values at end of each used/available step
plan(multisession, workers = availableCores() - 2)
tSSF.input2 <- extract.covars(data = tSSF.input, layers = cov_list, dyn_names = c('k490', 'npp', 'sst'),
                              along = TRUE, ind = "month.year", imputed = FALSE)
#takes 1.5 min to run on desktop (18 cores)
plan(sequential)


# Center and scale covariates
tSSF.input3 <- tSSF.input2 %>%
  mutate(strata = tSSF.input$strata,
         obs = tSSF.input$obs,
         step = tSSF.input$step) %>%
  drop_na(bathym, k490, npp, sst) %>%
  mutate(bathym.s = scale(bathym) %>%
           as.vector(),
         k490.s = scale(k490) %>%
           as.vector(),
         npp.s = scale(npp) %>%
           as.vector(),
         sst.s = scale(sst) %>%
           as.vector())

# Only keep IDs where N >= 30
tSSF.input4 <- tSSF.input3 %>%
  group_by(id) %>%
  filter(sum(obs) >= 30) %>%
  ungroup()



#############################################################
### Calculate probabilities from time model for each step ###
#############################################################


# Only keep strata that have at least 1 used and 1 available step (since only complete cases analyzed)
tSSF.input4 <- tSSF.input4 %>%
  group_by(strata) %>%
  filter(sum(obs) == 1) %>%
  filter(n() > 1) %>%
  ungroup()

# Add new ID column as consecutive numeric values
tSSF.input4 <- tSSF.input4 %>%
  mutate(id1 = as.numeric(factor(id)), .after = id)



tSSF.input.list <- split(tSSF.input4, tSSF.input4$id1)

plan(multisession, workers = availableCores() - 2)

tic()
progressr::with_progress({
  #set up progress bar
  p<- progressr::progressor(steps = length(tSSF.input.list))

res <- future_map(tSSF.input.list,
                  ~calc_time_probs(mod = mod, dat = ., covar.names = c('bathym.s','k490.s','npp.s','sst.s'), p = p),
                  .options = furrr_options(seed = 2022))
})
toc()  #took 2.5 min to run

plan(sequential)


tSSF.input4$time.prob <- unlist(res)
# time.mod.input3$time.prob <- ifelse(time.mod.input3$time.prob == 0, 1e-99, time.mod.input3$time.prob)
summary(tSSF.input4$time.prob)





########################################################
### Fit time-explicit step-selection function (tSSF) ###
########################################################

# Prepare data for model
tSSF.input4 <- tSSF.input4 %>%
  group_by(strata) %>%
  mutate(step.id = 1:n(), .after = obs) %>%  #number each of the used and available steps for subsetting
  ungroup()


# Replace values for missing steps per strata as 0s
plan(multisession, workers = availableCores() - 2)

tic()
tSSF.input5 <- tSSF.input4 %>%
  mutate(step.id = factor(step.id)) %>%
  split(.$id1) %>%
  future_map(., ~{.x %>%
      group_by(strata) %>%
      complete(step.id) %>%
      mutate(across(obs:time.prob, replace_na, 0)) %>%
      ungroup()
    }) %>%
  bind_rows()
toc()  #took 23 sec

plan(sequential)



# tSSF.input5 <- tSSF.input4 #%>%
  # slice(693391:693400) #%>%
  # slice(510301:510400)
  # slice(1:1000000)
  # filter(strata %in% sample(unique(tSSF.input4$strata), size = 1000, replace = FALSE))
  # filter(id %in% c(128352, 142659, 159781, 181796, 181800))

# Remove obs w/ time gap > 12 hours
# ind <- mod.input5 %>%
#   filter(obs == 1) %>%
#   filter(dt > 12*60) %>%
#   pull(strata)
#
# mod.input5 <- mod.input5 %>%
#   filter(!strata %in% ind)  #remove any steps w/ time gaps > 12 hrs

#get names of covariates
covar.names <- c('bathym.s','k490.s','npp.s','sst.s')


# Define arrays to store matrices of covars and time probs
xmat <- tSSF.input5 %>%
  dplyr::select(step.id, bathym.s, k490.s, npp.s, sst.s) %>%
  split(.$step.id) %>%
  map(~{.x %>%
      dplyr::select(-step.id) %>%
      as.matrix()})
xmat <- array(unlist(xmat), dim = c(n_distinct(tSSF.input5$strata), length(covar.names), max(tSSF.input4$step.id)))

pmov <- tSSF.input5 %>%
  dplyr::select(step.id, time.prob) %>%
  split(.$step.id) %>%
  map(~{.x %>%
      dplyr::select(-step.id) %>%
      as.matrix()})
pmov <- array(unlist(pmov), dim = c(n_distinct(tSSF.input5$strata), 1, max(tSSF.input4$step.id)))

# Create list of data
dat.list <- list(
  N = n_distinct(tSSF.input5$strata),
  K = length(covar.names),
  ID = tSSF.input5 %>%
    distinct(strata, .keep_all = TRUE) %>%
    pull(id1) %>%
    factor() %>% as.numeric(),
  nID = n_distinct(tSSF.input5$id, na.rm = TRUE),
  nstep = max(tSSF.input4$step.id),
  xmat = xmat,
  pmov = pmov,
  # xmat_used = as.matrix(tSSF.input5[tSSF.input5$step.id == 1, covar.names]),
  # xmat_avail1 = as.matrix(tSSF.input5[tSSF.input5$step.id == 2, covar.names]),
  # xmat_avail2 = as.matrix(tSSF.input5[tSSF.input5$step.id == 3, covar.names]),
  # xmat_avail3 = as.matrix(tSSF.input5[tSSF.input5$step.id == 4, covar.names]),
  # xmat_avail4 = as.matrix(tSSF.input5[tSSF.input5$step.id == 5, covar.names]),
  # pmov_used = tSSF.input5[tSSF.input5$step.id == 1,]$time.prob,
  # pmov_avail1 = tSSF.input5[tSSF.input5$step.id == 2,]$time.prob,
  # pmov_avail2 = tSSF.input5[tSSF.input5$step.id == 3,]$time.prob,
  # pmov_avail3 = tSSF.input5[tSSF.input5$step.id == 4,]$time.prob,
  # pmov_avail4 = tSSF.input5[tSSF.input5$step.id == 5,]$time.prob,
  y = rep(1, n_distinct(tSSF.input5$strata))
)





# Run model
stan.model <- "
data {
  int N;                                 // number of unique strata
  int K;                                 // number of covariates
  int nID;                               // number of unique IDs
  int nstep;                             // number of used and available steps
  int<lower=1, upper=nID> ID[N];         // ID label for each step
  matrix[K,nstep] xmat[N];                 // Design matrix for all used and available steps
  matrix<lower=0>[1,nstep] pmov[N];      // Time model probs for all used and available steps
  int<lower=1, upper=1> y[N];            // Ones for ones trick of Bernoulli model
}


parameters {
  //vector[K] betas;
  vector[K] mu;                          // Population means of slopes
  vector<lower=0>[K] sigma;              // Population SDs
  //vector<lower=0>[K] tau;                // Population scales (SDs)
  //cholesky_factor_corr[K] L;             // Population correlations
  //vector[K] betas[nID];                  // Slopes per individual
  vector[K] eta[nID];                    // Non-centered slopes per individual
}

transformed parameters {
  vector[K] betas[nID];

  for (i in 1:nID) {
    betas[i] = mu + eta[i] .* sigma;
  }
}


model {
    vector[nstep] prob[N];                // Matrix to store joint probabilities of time model and SSF
    vector[N] pi;                         // Vector to store tSSF probs

// priors
  //betas ~ normal(0, 1);
  mu ~ normal(0, 1);
  sigma ~ normal(0, 2);
  //tau ~ normal(0, 1);
  //L ~ lkj_corr_cholesky(4);

for (i in 1:nID) {
  eta[i] ~ normal(0, 1);
}

  // Estimate means per ID
  //betas ~ multi_normal_cholesky(mu, diag_pre_multiply(tau, L));
  //for (i in 1:nID) {
  //   betas[i] ~ normal(mu, sigma);
  // }


// model
  for (i in 1:N) {
    // calculation of pi
    prob[i] = pmov[i] * exp(xmat[i]' * betas[ID[i]]);

    pi[i] = prob[i,1] / sum(prob[i,]);


    // calculate likelihood using ones trick
    y[i] ~ bernoulli(pi[i]);
  }
}

"


mod1 <- stan(model_code = stan.model, data = dat.list, chains = 4, iter = 2000, warmup = 1000, seed = 2023,
             refresh = 100, control = list(adapt_delta = 0.99))
# took 1.25 hrs to run 2000 iter for full dataset

# params <- c('mu_b','b','b0_id','bBathym','bK490','bNPP','bSST','b0_bar','sd_b0')
print(mod1, digits_summary = 3, probs = c(0.025, 0.5, 0.975))


bayesplot::mcmc_neff(neff_ratio(mod1)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod1)) +
  bayesplot::yaxis_text(hjust = 0)
MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = c('mu','sigma'))
MCMCtrace(mod1, ind = TRUE, iter = 1000, pdf = FALSE, params = c('eta'))
par(mfrow=c(1,1))

pairs(mod1, pars = c('mu[1]','sigma[1]','eta[1,1]','betas[1,1]'))


## Save stanfit object
saveRDS(mod1, "Data_products/tSSF_model_GLMM_stanfit.rds")
