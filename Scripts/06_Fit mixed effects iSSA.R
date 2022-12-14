

### Fit iSSA w/ random effects to green turtles by behavioral state ###

library(tidyverse)
library(lubridate)
library(amt)
library(terra)
library(INLA)
library(future)
library(furrr)
library(tictoc)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

### NEED TO FILTER OUT ANY EXTREME STEP LENGTHS THAT MAY BE PRESENT (I.E., > 99TH QUANTILE)

dat <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foieGras.csv")

glimpse(dat)
summary(dat)


# Remove obs that fell w/in large (i.e., 3-day) gaps
dat2 <- dat %>%
  filter(!is.na(bout))


# Wrangle data into proper format for {amt}
dat2 <- dat2 %>%
  rename(t = date) %>%
  mutate(month.year = as_date(t), .after = 't') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))

dat_all <- dat2 %>%
  nest(data = -id) %>%
  mutate(trk = map(data, ~{
    make_track(.x, x, y, t, crs = sf::st_crs('EPSG:3395'))
  }
  ))

# Summarize time steps
dat_all %>%
  mutate(time.step = map(trk, summarize_sampling_rate)) %>%
  dplyr::select(id, time.step) %>%
  unnest(cols = time.step) %>%
  data.frame()
#most obs at exactly 2 hr time step; some large gaps due to previous filtering step



#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'k490', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('k490', 'npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
  # terra::time(cov_list[[var]]) <- as_date(names(cov_list[[var]]))
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')


# Convert to Raster format for extraction via {amt} and define time component
cov_list2 <- map(cov_list, brick)

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('k490', 'npp', 'sst')) {
  cov_list2[[var]] <- raster::setZ(cov_list2[[var]], as.Date(names(cov_list[[var]])))
}


####################################################################
### Generate random (available) steps and extract environ covars ###
####################################################################

dat_ssf <- dat_all %>%
  mutate(stps = map(trk, ~{.x %>%
      track_resample(rate = hours(2), tolerance = minutes(5)) %>%
                     steps_by_burst() %>%
                     random_steps(n_control = 10)})  #depending on how the iSSA model goes, potentially increase to 30 steps
  )

dat_ssf2 <- dat_ssf %>%
  dplyr::select(id, stps) %>%
  unnest(cols = stps) %>%
  rename(date = t1_) %>%
  mutate(month.year = as_date(date),
         .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))

names(dat_ssf2) <- gsub(pattern = "_$", replacement = "", x = names(dat_ssf2))

# Remove ID 104833 from dataset
dat_ssf2 <- dat_ssf2 %>%
  filter(id != 104833)  #no raster data available for K490 and NPP in 2011


plan(multisession, workers = availableCores() - 2)
dat_ssf3 <- extract.covars(data = dat_ssf2, layers = cov_list, dyn_names = c('k490','npp','sst'),
                           along = FALSE, ind = 'month.year')
# took 1.5 min to extract
plan(sequential)

dat_ssf4 <- dat_ssf3 %>%
  # mutate(stps2 = map(stps, ~{.x %>%
  #                     extract_covariates_var_time(covariates = cov_list2$sst,
  #                                                 when = "after",
  #                                                 where = "end",
  #                                                 max_time = days(2),
  #                                                 name_covar = names(cov_list2)[4])})
  # ) #%>%
  # select(id, stps) %>%
  # unnest(cols = stps) %>%
  mutate(y = as.numeric(case),
         id = as.numeric(factor(id)),
         step_id = paste(id, step_id, sep = "-"),
         cos_ta = cos(ta),
         log_sl = log(sl))

dat_ssf4





######################
### Fit iSSA model ###
######################


m1.ssf <- dat_ssf4 %>%
  fit_clogit(y ~ bathym + k490 + npp + sst + strata(step_id))

m2.ssf <- dat_ssf4 %>%
  fit_clogit(y ~ bathym + I(bathym^2) + k490 + I(k490^2) + npp + I(npp^2) + sst + I(sst^2) + strata(step_id))

summary(m1.ssf)
summary(m2.ssf)



m1.issa <- dat_ssf4 %>%
  fit_clogit(y ~ bathym + k490 + npp + sst + sl + cos_ta + log_sl + strata(step_id))

m2.issa <- dat_ssf4 %>%
  fit_clogit(y ~ bathym + I(bathym^2) + k490 + I(k490^2) + npp + I(npp^2) + sst + I(sst^2) + sl + cos_ta + log_sl + strata(step_id))

summary(m1.issa)
summary(m2.issa)





# To fit the model with random slopes in INLA, we need to generate new (but identical) variables of individual ID (ID cannot be used multiple times in the model formula):
dat_ssf4$id1 <- dat_ssf4$id
dat_ssf4$id2 <- dat_ssf4$id
dat_ssf4$id3 <- dat_ssf4$id
dat_ssf4$id4 <- dat_ssf4$id

# Set the model formula as for the fixed-effects model, but now add four random slope terms, namely for bathymetry, k490, npp, and sst. The priors for precision of the four random slopes are PC(3,0.05), while the intercept variance is again fixed:
formula.random <- y ~ -1 + bathym + k490 + npp + sst +
  # Breaks_Dis +
  f(step_id, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id1, bathym, values = unique(dat_ssf4$id), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
    f(id2, bathym, values = unique(dat_ssf4$id), model = "iid",
      hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
    f(id3, bathym, values = unique(dat_ssf4$id), model = "iid",
      hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
    f(id4, bathym, values = unique(dat_ssf4$id), model = "iid",
      hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05))))

# Fit the model
tic()
r.inla.random <- inla(formula.random, family = "Poisson", data = dat_ssf4,
                      control.fixed = list(
                        mean = 0,
                        prec = list(default = 1e-4)
                      ), verbose = TRUE
)
toc()  # took x min to run

# The summary for the posterior distribution of the fixed effects:
r.inla.random$summary.fixed

# Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions:
r.inla.random$summary.hyperpar


# Source R functions for calculating posterior means
# and medians of the precisions.
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(r.inla.random)
inla_mmarginal(r.inla.random)
