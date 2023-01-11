

### Fit iSSA w/ random effects to green turtles by behavioral state ###

library(tidyverse)
library(lubridate)
library(bayesmove)
library(terra)
library(INLA)
library(future)
library(furrr)
library(tictoc)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foieGras.csv")

glimpse(dat)
summary(dat)


dat <- dat %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))



###########################################################################
### Calculate step lengths, turning angles, and dt; filter observations ###
###########################################################################

dat <- prep_data(dat = dat, coord.names = c('x','y'), id = "id")
dat$speed <- dat$step / dat$dt


# Check speed for any outliers
round(quantile(dat$speed, c(0.01, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE), 2)
any(dat$step == 0, na.rm = TRUE)  #no step lengths equal to 0


# Remove observations w/ NA step length (i.e., last obs of each ID)
dat2 <- dat %>%
  filter(!is.na(step))


# Remove any observations w/ speeds faster than 3 m/s and w/ large time gaps (i.e., 3 days)
dat2 <- dat2 %>%
  filter(speed <= 3) %>%
  filter(!is.na(bout))




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

## Center and scale (by 1 SD) all covar layers
# cov_list_s <- sapply(cov_list, scale_across_covar)


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')
# cov_list_s <- map(cov_list_s, terra::project, 'EPSG:3395')




####################################################################
### Generate random (available) steps and extract environ covars ###
####################################################################

plan(multisession, workers = availableCores() - 2)
dat_ssf <- add_avail_steps(data = dat2, nsteps = 30)
plan(sequential)  #took 1 min to run


# Remove ID 104833 from dataset and remove rows where x2 and y2 are NA
dat_ssf2 <- dat_ssf %>%
  filter(id != 104833) %>%  #no raster data available for K490 and NPP in 2011
  filter(!is.na(x2) & !is.na(y2))


plan(multisession, workers = availableCores() - 2)
dat_ssf3 <- extract.covars(data = dat_ssf2, layers = cov_list, dyn_names = c('k490','npp','sst'),
                           along = FALSE, ind = 'month.year')
# took 1.5 min to extract
plan(sequential)


# Scale environ vars (before removing obs w/ missing values)
dat_ssf4 <- dat_ssf3 %>%
  mutate(bathym.s = as.numeric(scale(bathym)),
         k490.s = as.numeric(scale(k490)),
         npp.s = as.numeric(scale(npp)),
         sst.s = as.numeric(scale(sst)))


dat_ssf4 <- dat_ssf4 %>%
  mutate(step = step / 1000) %>%  #convert from m to km so on similar scale to standardized covars
  mutate(cos_ta = cos(angle),
         log_sl = log(step)) %>%
  drop_na(bathym, k490, npp, sst)  #remove all rows w/ any missing covar data

dat_ssf4


# Only keep strata that have at least 1 used and 1 available step w/ no missing covars
dat_ssf5 <- dat_ssf4 %>%
  group_by(strata) %>%
  filter(sum(obs) == 1) %>%
  filter(n() > 1) %>%
  ungroup()


# Check corr among covariates
cor(dat_ssf5[,c('bathym','k490','npp','sst')]) # all < 0.6


######################
### Fit iSSA model ###
######################

# To fit the model with random slopes in INLA, we need to generate new (but identical) variables of individual ID (ID cannot be used multiple times in the model formula):
dat_ssf5$id1 <- as.numeric(factor(dat_ssf5$id))
dat_ssf5$id2 <- dat_ssf5$id1
dat_ssf5$id3 <- dat_ssf5$id1
dat_ssf5$id4 <- dat_ssf5$id1
# dat_ssf5$id5 <- dat_ssf5$id1
# dat_ssf5$id6 <- dat_ssf5$id1
# dat_ssf5$id7 <- dat_ssf5$id1
# dat_ssf5$id8 <- dat_ssf5$id1

# Set the model formula as for the fixed-effects model, but now add four random slope terms, namely for bathymetry, k490, npp, and sst. The priors for precision of the four random slopes are PC(3,0.05), while the intercept variance is again fixed:
iSSA.formula <- obs ~ -1 + bathym.s + k490.s + npp.s + sst.s +
  step + log_sl + cos_ta +
  f(strata, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id1, bathym.s, values = unique(dat_ssf5$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
  # f(id2, I(bathym^2), values = unique(dat_ssf5$id1), model = "iid",
  #   hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id2, k490.s, values = unique(dat_ssf5$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
  # f(id4, I(k490^2), values = unique(dat_ssf5$id1), model = "iid",
  #   hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id3, npp.s, values = unique(dat_ssf5$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) +
  # f(id6, I(npp^2), values = unique(dat_ssf5$id1), model = "iid",
  #   hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id4, sst.s, values = unique(dat_ssf5$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(3,0.05)))) #+
  # f(id8, I(sst^2), values = unique(dat_ssf5$id1), model = "iid",
  #   hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))

# Fit the model
tic()
fit.iSSA <- inla(iSSA.formula, family = "Poisson", data = dat_ssf5,
                  control.fixed = list(
                    mean = 0,
                    prec = list(default = 1e-4)),
                  control.compute = list(waic = TRUE,
                                         dic = TRUE), verbose = FALSE
                  )
toc()  # took 2 min to run

summary(fit.iSSA)

# The summary for the posterior distribution of the fixed effects:
fixed.coeffs <- fit.iSSA$summary.fixed %>%
  mutate(param = factor(rownames(.), levels = rownames(.)))

ggplot(fixed.coeffs, aes(y = param)) +
  geom_vline(xintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
  geom_point(aes(x = `0.5quant`)) +
  theme_bw()


random.coeffs <- fit.iSSA$summary.random[-1] %>%
  bind_rows(.id = "param") %>%
  mutate(across(param, factor))
levels(random.coeffs$param) <- fixed.coeffs$param[1:4]

# Add population means to random effects
random.coeffs <- random.coeffs %>%
  split(.$param) %>%
  map2(.x = ., .y = fixed.coeffs$mean[1:4],
       ~mutate(.x,
               mean = mean + .y,
               `0.025quant` = `0.025quant` + .y,
               `0.975quant` = `0.975quant` + .y)) %>%
  bind_rows()
random.coeffs <- rbind(random.coeffs,
                       fixed.coeffs %>%
                         mutate(ID = "Pop", .before = mean) %>%
                         relocate(param, .before = ID)) %>%
  arrange(param)

ggplot(random.coeffs) +
  geom_hline(yintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, color = param)) +
  geom_point(aes(y = mean, x = ID, color = param)) +
  theme_bw() +
  facet_wrap(~ param, scales = "free")


# Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions
fit.iSSA$summary.hyperpar


# Source R functions for calculating posterior means
# and medians of the precisions.
inla_emarginal(fit.iSSA)
inla_mmarginal(fit.iSSA)




####################################################
### Viz marginal effects plots per environ covar ###
####################################################

fixed.coeffs2 <- fit.iSSA$summary.fixed[1:4,c("mean","0.025quant","0.975quant")] %>%
  as.matrix()

random.coeffs2 <- random.coeffs %>%
  filter(!param %in% c('step', 'log_sl', 'cos_ta')) %>%
  dplyr::select(param, ID, mean, `0.025quant`, `0.975quant`)


### Bathymetry ###

bathym.newdata <- data.frame(bathym = seq(min(dat_ssf5$bathym.s), max(dat_ssf5$bathym.s), length.out = 100),
                             k490 = 0,
                             npp = 0,
                             sst = 0) %>%
  as.matrix()


## Come back and calculate log-RSS for these results (or the avg effect of each covar) as discussed in Avgar et al 2017



random.coeffs2 <- random.coeffs %>%
  filter(!param %in% c('step', 'log_sl', 'cos_ta')) %>%
  dplyr::select(param, ID, mean, `0.025quant`, `0.975quant`)

pred.bathym <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.bathym) <- n_distinct(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- bathym.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(bathym = (bathym.newdata[,1] * sd(dat_ssf3$bathym, na.rm = T)) + mean(dat_ssf3$bathym, na.rm = T))

  pred.bathym[[i]] <- tmp
}

pred.bathym <- pred.bathym %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.bathym %>%
              filter(id != 48), aes(x = bathym, y = exp(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.bathym %>%
                filter(id == 48), aes(x = bathym, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.bathym %>%
              filter(id == 48), aes(x = bathym, y = exp(mean)), linewidth = 1.5) +
  theme_bw() +
  ylim(0,2)









### K490 ###

k490.newdata <- data.frame(bathym = 0,
                           k490 = seq(min(dat_ssf5$k490.s), max(dat_ssf5$k490.s), length.out = 100),
                           npp = 0,
                           sst = 0) %>%
  as.matrix()


## Come back and calculate log-RSS for these results (or the avg effect of each covar) as discussed in Avgar et al 2017



pred.k490 <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.k490) <- unique(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- k490.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(k490 = (k490.newdata[,2] * sd(dat_ssf3$k490, na.rm = T)) + mean(dat_ssf3$k490, na.rm = T))

  pred.k490[[i]] <- tmp
}

pred.k490 <- pred.k490 %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.k490 %>%
              filter(id != "Pop"), aes(x = k490, y = exp(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.k490 %>%
                filter(id == "Pop"), aes(x = k490, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.k490 %>%
              filter(id == "Pop"), aes(x = k490, y = exp(mean)), linewidth = 1.5) +
  theme_bw()









### NPP ###

npp.newdata <- data.frame(bathym = 0,
                           k490 = 0,
                           npp = seq(min(dat_ssf5$npp.s), max(dat_ssf5$npp.s), length.out = 100),
                           sst = 0) %>%
  as.matrix()


## Come back and calculate log-RSS for these results (or the avg effect of each covar) as discussed in Avgar et al 2017



pred.npp <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.npp) <- unique(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- npp.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(npp = (npp.newdata[,3] * sd(dat_ssf3$npp, na.rm = T)) + mean(dat_ssf3$npp, na.rm = T))

  pred.npp[[i]] <- tmp
}

pred.npp <- pred.npp %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.npp %>%
              filter(id != "Pop"), aes(x = npp, y = exp(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.npp %>%
                filter(id == "Pop"), aes(x = npp, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.npp %>%
              filter(id == "Pop"), aes(x = npp, y = exp(mean)), linewidth = 1.5) +
  theme_bw()









### SST ###

sst.newdata <- data.frame(bathym = 0,
                          k490 = 0,
                          npp = 0,
                          sst = seq(min(dat_ssf5$sst.s), max(dat_ssf5$sst.s), length.out = 100)) %>%
  as.matrix()


## Come back and calculate log-RSS for these results (or the avg effect of each covar) as discussed in Avgar et al 2017



pred.sst <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.sst) <- unique(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- sst.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(sst = (sst.newdata[,4] * sd(dat_ssf3$sst, na.rm = T)) + mean(dat_ssf3$sst, na.rm = T))

  pred.sst[[i]] <- tmp
}

pred.sst <- pred.sst %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.sst %>%
              filter(id != "Pop"), aes(x = sst, y = exp(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.sst %>%
                filter(id == "Pop"), aes(x = sst, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.sst %>%
              filter(id == "Pop"), aes(x = sst, y = exp(mean)), linewidth = 1.5) +
  theme_bw() +
  ylim(0,20)
