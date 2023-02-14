
### Fit RSF as HGAM ###

library(tidyverse)
library(lubridate)
library(mgcv)
library(terra)
library(INLA)
library(future)
library(furrr)
library(sf)
library(tictoc)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_30 <- read_csv("Processed_data/GoM_Cm_RSFprep_30x.csv")
rsf.pts_50 <- read_csv("Processed_data/GoM_Cm_RSFprep_50x.csv")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



####################
### Fit GLMM RSF ###
####################

# Center and scale covars; remove rows w/ incomplete observations
rsf.pts_10s <- rsf.pts_10 %>%
  drop_na(bathym, k490, npp, sst)
obs.ind_10 <- which(rsf.pts_10s$obs == 1)
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(bathym.s = (bathym - mean(bathym[obs.ind_10])) / sd(bathym),
         k490.s = (k490 - mean(k490[obs.ind_10])) / sd(k490),
         npp.s = (npp - mean(npp[obs.ind_10])) / sd(npp),
         sst.s = (sst - mean(sst[obs.ind_10])) / sd(sst))

rsf.pts_30s <- rsf.pts_30 %>%
  drop_na(bathym, k490, npp, sst)
obs.ind_30 <- which(rsf.pts_30s$obs == 1)
rsf.pts_30s <- rsf.pts_30s %>%
  mutate(bathym.s = (bathym - mean(bathym[obs.ind_30])) / sd(bathym),
         k490.s = (k490 - mean(k490[obs.ind_30])) / sd(k490),
         npp.s = (npp - mean(npp[obs.ind_30])) / sd(npp),
         sst.s = (sst - mean(sst[obs.ind_30])) / sd(sst))

rsf.pts_50s <- rsf.pts_50 %>%
  drop_na(bathym, k490, npp, sst)
obs.ind_50 <- which(rsf.pts_50s$obs == 1)
rsf.pts_50s <- rsf.pts_50s %>%
  mutate(bathym.s = (bathym - mean(bathym[obs.ind_50])) / sd(bathym),
         k490.s = (k490 - mean(k490[obs.ind_50])) / sd(k490),
         npp.s = (npp - mean(npp[obs.ind_50])) / sd(npp),
         sst.s = (sst - mean(sst[obs.ind_50])) / sd(sst))



# Infinitely-weighted logistic regression
rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, 5000, 1)
rsf.pts_30s$wts <- ifelse(rsf.pts_30s$obs == 0, 5000, 1)
rsf.pts_50s$wts <- ifelse(rsf.pts_50s$obs == 0, 5000, 1)


# Explore used vs available habitat values
rsf.pts_10s %>%
  # mutate(across(c(k490, npp), log)) %>%
  pivot_longer(cols = c(bathym, k490, npp, sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")

# Log-transform skewed covars to allow model fitting
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_30s <- rsf.pts_30s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_50s <- rsf.pts_50s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.k490 = log(k490),
         log.npp = log(npp),
         log.sst = log(sst))


# Check Pearson corrs
cor(rsf.pts_10s[,c('bathym','k490','npp','sst')])  #NPP and K490 highly corr (0.78); remove k490


# Now explore transformed distributions
rsf.pts_10s %>%
  pivot_longer(cols = c(log.bathym, log.npp, log.sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)
rsf.pts_30s$wts2 <- ifelse(rsf.pts_30s$obs == 0, A / sum(rsf.pts_30s$obs == 0), 1e-6)
rsf.pts_50s$wts2 <- ifelse(rsf.pts_50s$obs == 0, A / sum(rsf.pts_50s$obs == 0), 1e-6)


## Mixed RSF via INLA
# rsf.pts_10s$id1 <- as.numeric(factor(rsf.pts_10s$id))
# rsf.pts_10s <- arrange(rsf.pts_10s, id1)


# Generate formula in {mgcv}
rsf.pts_10s2 <- rsf.pts_10s %>%
  mutate(across(id, factor)) #%>%
  # filter(id %in% c(128352, 181800, 181796))

# rsf.pts_30s2 <- rsf.pts_30s %>%
#   mutate(across(id, factor))
#
# rsf.pts_50s2 <- rsf.pts_50s %>%
#   mutate(across(id, factor))

set.seed(2023)
tic()
fit.HGAM10_PI <- bam(obs/wts2 ~ s(log.bathym, bs = "cr", k = 5, m = 2) +
                    s(log.bathym, by = id, bs = "cr", k = 5, m = 1) +
                    s(log.npp, bs = "cr", k = 5, m = 2) +
                    s(log.npp, by = id, bs = "cr", k = 5, m = 1) +
                    s(log.sst, bs = "cr", k = 5, m = 2) +
                    s(log.sst, by = id, bs = "cr", k = 5, m = 1) +
                    s(id, bs = "re"), data = rsf.pts_10s2, method = "fREML",
                  family = poisson(), weights = wts2, discrete = TRUE)
toc()  #took 16 min to run

summary(fit.HGAM10_PI)
plot(fit.HGAM10_PI, select = 1, scale = 0, shade = TRUE, shade.col = "lightblue")
plot(fit.HGAM10_PI, select = 50, scale = 0, shade = TRUE, shade.col = "lightblue")
plot(fit.HGAM10_PI, select = 99, scale = 0, shade = TRUE, shade.col = "lightblue")
# plot(fit.HGAM10, scale = 0, shade = TRUE, shade.col = "lightblue")
gam.check(fit.HGAM10)
concurvity(fit.HGAM10)


