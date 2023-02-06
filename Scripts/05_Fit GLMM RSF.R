
### Fit RSF as GLMM ###

library(tidyverse)
library(lubridate)
library(bayesmove)
library(terra)
# library(INLA)
library(future)
library(furrr)
library(sf)
library(sfarrow)
library(tictoc)
library(amt)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

dat <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_4hr_aniMotum.csv")
gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

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


# Remove any observations w/ speeds faster than 3 m/s and w/ large time gaps (i.e., 3 days)
dat2 <- dat %>%
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




#####################################################################
### Generate random (available) points and extract environ covars ###
#####################################################################

# Remove points that intersect land
tic()
dat3 <- dat2 %>%
  split(.$id) %>%
  map(~{.x %>%
      st_as_sf(coords = c('x','y'), crs = 3395) %>%
      st_mask(gom.sf) %>%
      mutate(x = st_coordinates(.)[,1],
             y = st_coordinates(.)[,2]) %>%
      st_drop_geometry()}) %>%
  bind_rows()
toc()  #took 17 sec

# How many points per ID?
table(dat3$id)  #min of 37; max of 1953

# Convert to class for use of {amt} functions
dat4 <- make_track(dat3, .x = x, .y = y, .t = date, id = id, month.year = month.year, crs = 3395) %>%
  nest(data = -id)

# Create KDE_href for each ID and add as nested column
dat4 <- dat4 %>%
  mutate(ud = map(data, ~hr_akde(.x, levels = 0.99)))

# Add buffer to extend available area
mean(dat3$step) * 6  #check avg distance that could be covered in a day (based on 4 hr time step)

dat4 <- dat4 %>%
  mutate(ud_buff = map(ud, ~{.x %>%
      hr_isopleths() %>%
      slice(2) %>%
      sf::st_buffer(dist = 1e4)}  #add 10 km buffer (avg daily distance)
      ))

# Clip UDs by land mask
dat4 <- dat4 %>%
  mutate(ud_buff = map(ud_buff, ~st_mask(.x, gom.sf)))


# Generate available points and randomly assign to month.year per ID
dat4 <- dat4 %>%
  mutate(avail_pts = map2(.x = ud_buff,
                          .y = data,
                          ~{data.frame(geometry = st_sample(.x, size = 10 * nrow(.y))) %>%
                              mutate(month.year = sample(unique(.y$month.year), size = n(), replace = T)) %>%
                              st_sf()}
                          ))


# Viz example of available point spread by month.year
ggplot() +
  geom_sf(data = dat4$avail_pts[[45]], aes(color = month.year)) +
  theme_bw() +
  facet_wrap(~ month.year)


# Join all used and available points together
used <- dat4 %>%
  dplyr::select(id, data) %>%
  unnest(cols = data) %>%
  mutate(obs = 1) %>%
  rename(x = x_, y = y_) %>%
  dplyr::select(id, month.year, x, y, obs) %>%
  data.frame()

avail <- dat4 %>%
  dplyr::select(id, avail_pts) %>%
  mutate(avail_pts = map(avail_pts, ~{.x %>%
      mutate(x = st_coordinates(.)[,1],
             y = st_coordinates(.)[,2]) %>%
      st_drop_geometry()}
      )) %>%
  unnest(cols = avail_pts) %>%
  mutate(obs = 0)


rsf.pts <- rbind(used, avail)


# Remove ID 104833 since some covar data not available in 2011
rsf.pts <- rsf.pts %>%
  filter(id != 104833)


# Extract environ covars by month.year
plan(multisession, workers = availableCores() - 2)
rsf.pts <- extract.covars(data = rsf.pts, layers = cov_list, dyn_names = c('k490', 'npp', 'sst'),
                           along = FALSE, ind = "month.year", imputed = FALSE)
#takes 1.5 min to run on desktop (18 cores)
plan(sequential)



# Viz example of available point covar values by month.year
ggplot() +
  geom_point(data = rsf.pts %>%
            filter(id == 128352, obs == 0), aes(x, y, color = bathym)) +
  geom_point(data = rsf.pts %>%
               filter(id == 128352, obs == 1), aes(x, y), color = 'red') +
  theme_bw() +
  facet_wrap(~ month.year)

ggplot() +
  geom_point(data = rsf.pts %>%
               filter(id == 128352, obs == 0), aes(x, y, color = sst)) +
  geom_point(data = rsf.pts %>%
               filter(id == 128352, obs == 1), aes(x, y), color = 'red') +
  scale_color_viridis_c(option = 'inferno') +
  theme_bw() +
  facet_wrap(~ month.year)






####################
### Fit GLMM RSF ###
####################

# Center and scale covars; remove rows w/ incomplete observations
rsf.pts <- rsf.pts %>%
  mutate(bathym.s = as.numeric(scale(bathym)),
         k490.s = as.numeric(scale(k490)),
         npp.s = as.numeric(scale(npp)),
         sst.s = as.numeric(scale(sst))) %>%
  drop_na(bathym, k490, npp, sst)



# Infinitely-weighted logistic regression
rsf.pts$wts <- ifelse(rsf.pts$obs == 0, 5000, 1)

logis.mod <- glm(obs ~ bathym + k490 + npp + sst,
                data = rsf.pts, family = binomial(), weights = wts)

summary(logis.mod)


# Down-weighted Poisson regression
A <- (res(cov_list$bathym)[1] / 1000) ^ 2  #area of grid cells
rsf.pts$wts2 <- ifelse(rsf.pts$obs == 0, A/sum(rsf.pts$obs == 0), 1e-6)

pois.mod <- glm(obs/wts2 ~ bathym + k490 + npp + sst,
                data = rsf.pts, family = poisson(), weights = wts2)

summary(pois.mod)





## Mixed RSF via glmmTMB
library(glmmTMB)
logis.glmmTMB <- glmmTMB(obs ~ bathym.s + k490.s + npp.s + sst.s +
                        (1|id) + (0 + bathym|id)  + (0 + k490|id) + (0 + npp|id) + (0 + sst|id),
                      family = binomial(), data = rsf.pts,
                      doFit = FALSE, weights = wts)

# Fix standard deviation of random intercept to 1000 (variance is 1e6)
logis.glmmTMB$parameters$theta[1] <- log(1e3)

# Tell glmmTMB to leave sd of random intercept fixed, but estimate for random slopes
logis.glmmTMB$mapArg <- list(theta = factor(c(NA, 1, 1, 1, 1)))

# Fit model
logis.glmmTMB.fit <- fitTMB(logis.glmmTMB)
summary(logis.glmmTMB.fit)




## Mixed RSF via INLA
library(INLA)

rsf.pts$id1 <- as.numeric(factor(rsf.pts$id))
rsf.pts$id2 <- rsf.pts$id1
rsf.pts$id3 <- rsf.pts$id1
rsf.pts$id4 <- rsf.pts$id1
rsf.pts$id5 <- rsf.pts$id1
rsf.pts$id6 <- rsf.pts$id1
rsf.pts$id7 <- rsf.pts$id1
rsf.pts$id8 <- rsf.pts$id1
rsf.pts$id9 <- rsf.pts$id1

rsf.pts <- arrange(rsf.pts, id1)

RSF.formula <- obs ~ bathym.s + I(bathym.s ^ 2) + k490.s + I(k490.s ^ 2) + npp.s + I(npp.s ^ 2) + sst.s + I(sst.s ^ 2) +
  f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id2, bathym.s, values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id3, I(bathym.s ^ 2), values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id4, k490.s, values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id5, I(k490.s ^ 2), values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id6, npp.s, values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id7, I(npp.s ^ 2), values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id8, sst.s, values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id9, I(sst.s ^ 2), values = unique(rsf.pts$id1), model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))

# Fit the model
tic()
fit.RSF <- inla(RSF.formula, family = "binomial", data = rsf.pts, #weights = rsf.pts$wts,
                 control.fixed = list(
                   mean = 0,
                   prec = list(default = 1e-4)),
                 control.compute = list(waic = TRUE,
                                        dic = TRUE), verbose = FALSE
)
toc()  # took 11.5 min to run

summary(fit.RSF)




# The summary for the posterior distribution of the fixed effects:
fixed.coeffs <- fit.RSF$summary.fixed %>%
  dplyr::slice(2:n()) %>%  #remove intercept
  mutate(param = factor(rownames(.), levels = rownames(.)))

ggplot(fixed.coeffs, aes(y = param)) +
  geom_vline(xintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(xmin = `0.025quant`, xmax = `0.975quant`)) +
  geom_point(aes(x = mean)) +
  theme_bw()


random.coeffs <- fit.RSF$summary.random[-1] %>%
  bind_rows(.id = "param") %>%
  mutate(across(param, factor))
levels(random.coeffs$param) <- fixed.coeffs$param

# Add population means to random effects
random.coeffs <- random.coeffs %>%
  split(.$param) %>%
  map2(.x = ., .y = fixed.coeffs$mean,
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




####################################################
### Viz marginal effects plots per environ covar ###
####################################################

fixed.coeffs2 <- fixed.coeffs[,c("mean","0.025quant","0.975quant")] %>%
  as.matrix()

random.coeffs2 <- random.coeffs %>%
  dplyr::select(param, ID, mean, `0.025quant`, `0.975quant`)


### Bathymetry ###

bathym.newdata <- data.frame(bathym = seq(min(rsf.pts$bathym.s), max(rsf.pts$bathym.s), length.out = 100),
                             bathym.2 = seq(min(rsf.pts$bathym.s), max(rsf.pts$bathym.s), length.out = 100) ^ 2,
                             k490 = 0,
                             k490.s = 0,
                             npp = 0,
                             npp.2 = 0,
                             sst = 0,
                             sst.2 = 0) %>%
  as.matrix()


## Come back and calculate log-RSS for these results (or the avg effect of each covar) as discussed in Avgar et al 2017



pred.bathym <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.bathym) <- unique(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- bathym.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(bathym = (bathym.newdata[,1] * sd(rsf.pts$bathym, na.rm = T)) + mean(rsf.pts$bathym, na.rm = T))

  pred.bathym[[i]] <- tmp
}

pred.bathym <- pred.bathym %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.bathym %>%
              filter(id != "Pop"), aes(x = bathym, y = plogis(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.bathym %>%
                filter(id == "Pop"), aes(x = bathym, ymin = plogis(X0.025quant), ymax = plogis(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.bathym %>%
              filter(id == "Pop"), aes(x = bathym, y = plogis(mean)), linewidth = 1.5) +
  theme_bw()







### SST ###

sst.newdata <- data.frame(bathym = 0,
                             bathym.2 = 0,
                             k490 = 0,
                             k490.s = 0,
                             npp = 0,
                             npp.2 = 0,
                             sst = seq(min(rsf.pts$sst.s), max(rsf.pts$sst.s), length.out = 100),
                             sst.2 = seq(min(rsf.pts$sst.s), max(rsf.pts$sst.s), length.out = 100) ^ 2) %>%
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
    mutate(sst = (sst.newdata[,7] * sd(rsf.pts$sst, na.rm = T)) + mean(rsf.pts$sst, na.rm = T))

  pred.sst[[i]] <- tmp
}

pred.sst <- pred.sst %>%
  bind_rows(.id = "id")

# Pop mean in black; ID by color
ggplot() +
  geom_line(data = pred.sst %>%
              filter(id != "Pop"), aes(x = sst, y = plogis(mean), group = id, color = id), linewidth = 0.75, show.legend = FALSE) +
  geom_ribbon(data = pred.sst %>%
                filter(id == "Pop"), aes(x = sst, ymin = plogis(X0.025quant), ymax = plogis(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.sst %>%
              filter(id == "Pop"), aes(x = sst, y = plogis(mean)), linewidth = 1.5) +
  theme_bw()
