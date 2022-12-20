
### Explore results of time model

library(tidyverse)
library(lubridate)
library(rstan)
library(MCMCvis)
library(bayesplot)
# library(arrow)
library(terra)
library(sfarrow)
library(sf)

source('Scripts/helper functions.R')


### Load model fit, model input, and tracks ###

mod <- readRDS('Data_products/Time_model_intercept-slopes_stanfit.rds')
mod.input <- read_csv("Processed_data/Input for time model.csv")
dat <- read_csv('Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foieGras.csv')

gom <- sfarrow::st_read_parquet("Environ_data/GoM_land.parquet") %>%
  st_transform(3395)


### Load environmental raster layers to make predictions of time to cross pixel ###

files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[!grepl(pattern = "parquet", files)]  #remove any parquet files

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

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



##--------------------------------------------------------




### Wrangle tracking data ###
dat <- dat %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))


### Wrangle model input to required format ###

# Change time step from secs to mins
mod.input$dt <- mod.input$dt/60

# Retain only observed steps
# mod.input <- mod.input %>%
#   filter(obs == 1)

# Remove all observations w/ missing values (since not imputed)
mod.input2 <- mod.input %>%
  drop_na(bathym, k490, npp, sst)


# Center and scale covariates
# mod.input2 <- mod.input2 %>%
#   mutate(bathym.s = scale(bathym) %>%
#            as.vector(),
#          kd490.s = scale(Kd490) %>%
#            as.vector(),
#          npp.s = scale(NPP) %>%
#            as.vector(),
#          sst.s = scale(SST) %>%
#            as.vector())






### Explore model diagnostics and convergence ###
params <- c('mu_b','b','mu_betas','tau','L')
print(mod, digits_summary = 3, pars = params, probs = c(0.025, 0.5, 0.975))
print(mod, digits_summary = 3, pars = 'betas', probs = c(0.025, 0.5, 0.975))


MCMCtrace(mod, ind = TRUE, iter = 2000, pdf = FALSE, params = params)
par(mfrow=c(1,1))
MCMCplot(mod, params = params)

bayesplot::mcmc_neff(neff_ratio(mod, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)
bayesplot::mcmc_rhat(rhat(mod, pars = params)) +
  bayesplot::yaxis_text(hjust = 0)






### Viz distributions of parameters ###

# Posterior distribs of main coeffs
mcmc_areas(mod, regex_pars = 'betas', prob = 0.5, prob_outer = 1) +
  yaxis_text(size = 12) +
  xaxis_text(size = 12) +
  ggtitle("Posterior distributions of coefficients for environmental variables") +
  theme(plot.title = element_text(size = 20))

# Posterior distribs of intercept (including pop mean and sd)
# mcmc_intervals(mod, regex_pars = "b0", prob = 0.5, prob_outer = 0.95) +
#   yaxis_text(size = 12) +
#   xaxis_text(size = 12) +
#   ggtitle("Posterior distributions of intercept for model") +
#   theme(plot.title = element_text(size = 20))

# Posterior distribs of rate ('b') params for gamma distribution (including pop mean)
mcmc_intervals(mod, pars = "mu_b", regex_pars = "^b\\[", prob = 0.5, prob_outer = 0.95) +
  yaxis_text(size = 12) +
  xaxis_text(size = 12) +
  ggtitle("Posterior distributions of rate parameter for gamma distribution") +
  theme(plot.title = element_text(size = 20))



### Viz marginal effects plots to show relationships by ID ###

n <- 100
id1 <- unique(mod.input2$id)
dist1 <- median(mod.input2$dist)  #median step length is 474 m
b0 <- rstan::extract(mod, pars = 'betas')$betas[,,1]
bBathym <- rstan::extract(mod, pars = 'betas')$betas[,,2]
bK490 <- rstan::extract(mod, pars = 'betas')$betas[,,3]
bNPP <- rstan::extract(mod, pars = 'betas')$betas[,,4]
bSST <- rstan::extract(mod, pars = 'betas')$betas[,,5]



## Bathymetry
bathym.res<- list()
mu.bathym <- mean(terra::values(cov_list$bathym), na.rm = TRUE)
sd.bathym <- sd(terra::values(cov_list$bathym), na.rm = TRUE)

for (i in 1:length(id1)) {

  tmp<- mod.input2 %>%
    filter(id == id1[i])

  #Generate sequence along bathymetry
  rango1<- tmp %>%
    dplyr::select(bathym) %>%
    range()
  seq.bathym<- seq(rango1[1], rango1[2], length.out = n)


  # Create object to store results
  bathy.pred <- matrix(NA, n, 3)

  # Make predictions (while holding K490, NPP, and SST at 0)
  for (j in 1:n) {
    linear.func <- b0[,i] + bBathym[,i] * seq.bathym[j]
    mu <- dist1 * exp(linear.func)
    bathy.pred[j,] <- quantile(mu, c(0.025, 0.5, 0.975))
  }

  # Convert to data.frame and add additional vars
  bathy.pred <- data.frame(bathy.pred) %>%
    rename(lower = X1, med = X2, upper = X3) %>%
    mutate(bathy = (seq.bathym * sd.bathym) + mu.bathym,
           id = id1[i])

  bathym.res[[i]] <- bathy.pred
}

bathym.res.df<- bind_rows(bathym.res)

# Plot relationship
ggplot(data = bathym.res.df, aes(x = bathy)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(id)), alpha =  0.3) +
  geom_line(aes(y = med, color = factor(id)), linewidth = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Bathymetric depth (m)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.position = "none")

ggplot(data = bathym.res.df, aes(x = bathy)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = id), alpha =  0.3) +
  geom_line(aes(y = med, color = id), size = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Bathymetric depth (m)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  xlim(-500,0) +
  ylim(0,240)




## K490

k490.res<- list()

for (i in 1:length(id1)) {

  tmp<- mod.input2 %>%
    filter(id == id1[i])

  #Generate sequence along bathymetry
  rango1<- tmp %>%
    dplyr::select(k490) %>%
    range(na.rm = TRUE)
  seq.k490<- seq(rango1[1], rango1[2], length.out = n)


  # Create object to store results
  k490.pred <- matrix(NA, n, 3)

  # Make predictions (while holding bathym, NPP, and SST at 0)
  for (j in 1:n) {
    linear.func <- b0[,i] + bK490[,i] * seq.k490[j]
    mu <- dist1 * exp(linear.func)
    k490.pred[j,] <- quantile(mu, c(0.025, 0.5, 0.975))
  }

  # Convert to data.frame and add additional vars
  k490.pred <- data.frame(k490.pred) %>%
    rename(lower = X1, med = X2, upper = X3) %>%
    mutate(k490 = seq.k490 * sd(mod.input2$k490, na.rm = TRUE) + mean(mod.input2$k490, na.rm = TRUE),
           id = id1[i])

  k490.res[[i]] <- k490.pred
}

k490.res.df<- bind_rows(k490.res)

# Plot relationship
ggplot(data = k490.res.df, aes(x = k490)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(id)), alpha =  0.3) +
  geom_line(aes(y = med, color = factor(id)), linewidth = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'K490 (1/m)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.position = "none")

ggplot(data = k490.res.df, aes(x = k490)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = id), alpha =  0.3) +
  geom_line(aes(y = med, color = id), size = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'K490 (1/m)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  xlim(0,1.5) +
  ylim(0,240)




## NPP

npp.res<- list()

for (i in 1:length(id1)) {

  tmp<- mod.input2 %>%
    filter(id == id1[i])

  #Generate sequence along bathymetry
  rango1<- tmp %>%
    dplyr::select(npp) %>%
    range(na.rm = TRUE)
  seq.npp<- seq(rango1[1], rango1[2], length.out = n)


  # Create object to store results
  npp.pred <- matrix(NA, n, 3)

  # Make predictions (while holding bathym, NPP, and SST at 0)
  for (j in 1:n) {
    linear.func <- b0[,i] + bNPP[,i] * seq.npp[j]
    mu <- dist1 * exp(linear.func)
    npp.pred[j,] <- quantile(mu, c(0.025, 0.5, 0.975))
  }

  # Convert to data.frame and add additional vars
  npp.pred <- data.frame(npp.pred) %>%
    rename(lower = X1, med = X2, upper = X3) %>%
    mutate(npp = seq.npp * sd(mod.input2$npp, na.rm = TRUE) + mean(mod.input2$npp, na.rm = TRUE),
           id = id1[i])

  npp.res[[i]] <- npp.pred
}

npp.res.df<- bind_rows(npp.res)

# Plot relationship
ggplot(data = npp.res.df, aes(x = npp)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(id)), alpha =  0.3) +
  geom_line(aes(y = med, color = factor(id)), linewidth = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Net Primary Productivity (mg C m^-2 d^-1)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.position = "none")

ggplot(data = npp.res.df, aes(x = npp)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = id), alpha =  0.3) +
  geom_line(aes(y = med, color = id), size = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Net Primary Productivity (mg C m^-2 d^-1)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  xlim(0,10000) +
  ylim(0,1000)





## Sea surface temperature

sst.res<- list()

for (i in 1:length(id1)) {

  tmp<- mod.input2 %>%
    filter(id == id1[i])

  #Generate sequence along bathymetry
  rango1<- tmp %>%
    dplyr::select(sst) %>%
    range(na.rm = TRUE)
  seq.sst<- seq(rango1[1], rango1[2], length.out = n)


  # Create object to store results
  sst.pred <- matrix(NA, n, 3)

  # Make predictions (while holding bathym and SST at 0)
  for (j in 1:n) {
    linear.func <- b0[,i] + bSST[,i] * seq.sst[j]
    mu <- dist1 * exp(linear.func)
    sst.pred[j,] <- quantile(mu, c(0.025, 0.5, 0.975))
  }

  # Convert to data.frame and add additional vars
  sst.pred <- data.frame(sst.pred) %>%
    rename(lower = X1, med = X2, upper = X3) %>%
    mutate(sst = seq.sst * sd(mod.input2$sst, na.rm = TRUE) + mean(mod.input2$sst, na.rm = TRUE),
           id = id1[i])

  sst.res[[i]] <- sst.pred
}

sst.res.df<- bind_rows(sst.res)

# Plot relationship
ggplot(data = sst.res.df, aes(x = sst)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(id)), alpha =  0.3) +
  geom_line(aes(y = med, color = factor(id)), linewidth = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Sea surface temperature (°C)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18),
        legend.position = "none")

ggplot(data = sst.res.df, aes(x = sst)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = id), alpha =  0.3) +
  geom_line(aes(y = med, color = id), size = 1, alpha = 0.5) +
  scale_color_viridis_d("ID", option = "mako") +
  scale_fill_viridis_d("ID", option = "mako") +
  labs(x = 'Sea surface temperature (°C)', y = "Time to traverse 474 m (min)") +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18)) +
  ylim(0,240)






### Viz time surfaces for some IDs ###

# Inspect map to determine which IDs to highlight

ggplot(data = dat, aes(x, y, color = factor(id), group = id)) +
  geom_path(size = 0.25, alpha = 0.7) +
  theme_bw()


ggplot(data = dat %>% filter(id == 142713), aes(x, y, color = id, group = id)) +
  geom_path(size = 0.25, alpha = 0.7) +
  theme_bw() +
  coord_equal()


ggplot(data = dat %>% filter(id == 181800), aes(x, y, color = id, group = id)) +
  geom_path(size = 0.25, alpha = 0.7) +
  theme_bw() +
  coord_equal()


ggplot(data = dat %>% filter(id == 159783), aes(x, y, color = id, group = id)) +
  geom_path(size = 0.25, alpha = 0.7) +
  theme_bw() +
  coord_equal()


ggplot(data = dat %>% filter(id == 181796), aes(x, y, color = id, group = id)) +
  geom_path(size = 0.25, alpha = 0.7) +
  theme_bw() +
  coord_equal()



## Iterate predictions over IDs

focal.ids <- c(142713, 181800, 159783, 181796)
# bBathym.quant <- quantile(bBathym, c(0.025, 0.5, 0.975))
time.pred.list <- list()

for (i in 1:length(focal.ids)) {

  #define dat range for ID
  date.range <- dat %>%
    filter(id == focal.ids[i]) %>%
    mutate(month.year = as_date(month.year)) %>%
    dplyr::pull(month.year) %>%
    unique() %>%
    as.character()

  #define median step length
  med.dist <- median(mod.input2[mod.input2$id == focal.ids[i],]$dist, na.rm = TRUE)


  #subset environ rasters and scale values
  bathym <- ((cov_list$bathym - mean(mod.input2$bathym)) / sd(mod.input2$bathym))
  k490 <- ((cov_list$k490[[date.range]] - mean(mod.input2$k490, na.rm = TRUE)) / sd(mod.input2$k490, na.rm = TRUE))
  npp <- ((cov_list$npp[[date.range]] - mean(mod.input2$npp, na.rm = TRUE)) / sd(mod.input2$npp, na.rm = TRUE))
  sst <- ((cov_list$sst[[date.range]] - mean(mod.input2$sst, na.rm = TRUE)) / sd(mod.input2$sst, na.rm = TRUE))


  #calculate time
  ind <- which(id1 == focal.ids[i])
  linfunc <- median(b0[,ind]) + bathym * median(bBathym[,ind]) + k490 * median(bK490[,ind]) +
    npp * median(bNPP[,ind]) + sst * median(bSST[,ind])
  mu <- med.dist * exp(linfunc)
  names(mu) <- date.range


  # Convert to data.frame and add additional vars
  time.pred <- as.data.frame(mu, xy = TRUE) %>%
    pivot_longer(cols = -c(x, y), names_to = 'month.year', values_to = "time") %>%
    mutate(id = focal.ids[i])

  time.pred.list[[i]] <- time.pred
}

time.pred.df <- bind_rows(time.pred.list)



track.142713 <- dat %>%
  filter(id == 142713)

ggplot() +
  geom_raster(data = time.pred.df %>%
                filter(id == 142713), aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time (min)", option = 'inferno', limits = c(0,60)) +
  geom_sf(data = gom) +
  geom_path(data = track.142713, aes(x, y, group = id), color = 'dodgerblue4', size = 0.25, alpha = 0.25) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "PTT 142713") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  coord_sf(xlim = c(min(track.142713$x) - 10000, max(track.142713$x) + 10000),
         ylim = c(min(track.142713$y) - 10000, max(track.142713$y) + 10000)) +
  scale_x_continuous(breaks = c(-86.9,-86.7)) +
  facet_wrap(~ month.year)



track.181800 <- dat %>%
  filter(id == 181800)
unique(track.181800$month.year)

ggplot() +
  geom_raster(data = time.pred.df %>%
                filter(id == 181800), aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time (min)", option = 'inferno') +
  geom_sf(data = gom) +
  geom_path(data = track.181800, aes(x, y, group = id), color = 'chartreuse', size = 0.25, alpha = 0.25) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "PTT 181800") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  coord_sf(xlim = c(min(track.181800$x) - 10000, max(track.181800$x) + 10000),
           ylim = c(min(track.181800$y) - 10000, max(track.181800$y) + 10000)) +
  scale_x_continuous(breaks = c(-88, -86)) +
  facet_wrap(~ month.year, nrow = 1)



track.159783 <- dat %>%
  filter(id == 159783)
unique(track.159783$month.year)

ggplot() +
  geom_raster(data = time.pred.df %>%
                filter(id == 159783), aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time (min)", option = 'inferno') +
  geom_sf(data = gom) +
  geom_path(data = track.159783, aes(x, y, group = id), color = 'chartreuse', size = 0.25, alpha = 0.25) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "PTT 159783") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  coord_sf(xlim = c(min(track.159783$x) - 10000, max(track.159783$x) + 10000),
           ylim = c(min(track.159783$y) - 10000, max(track.159783$y) + 10000)) +
  # scale_x_continuous(breaks = c(-88, -86)) +
  facet_wrap(~ month.year)



track.181796 <- dat %>%
  filter(id == 181796)
unique(track.181796$month.year)

ggplot() +
  geom_raster(data = time.pred.df %>%
                filter(id == 181796), aes(x, y, fill = time)) +
  scale_fill_viridis_c("Time (min)", option = 'inferno') +
  geom_sf(data = gom) +
  geom_path(data = track.181796, aes(x, y, group = id), color = 'chartreuse', size = 0.25, alpha = 0.25) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", title = "PTT 181796") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 12, face = "bold")) +
  coord_sf(xlim = c(min(track.181796$x) - 10000, max(track.181796$x) + 10000),
           ylim = c(min(track.181796$y) - 10000, max(track.181796$y) + 10000)) +
  scale_x_continuous(breaks = c(-85, -84)) +
  facet_wrap(~ month.year)



