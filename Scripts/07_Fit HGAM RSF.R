
### Fit RSF as HGAM ###

library(tidyverse)
library(lubridate)
library(mgcv)
library(gratia)
library(terra)
library(INLA)
library(future)
library(furrr)
library(sf)
library(tictoc)
library(patchwork)

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



# Check concurvity from base model
set.seed(2023)
tic()
fit.GAM10 <- bam(obs/wts2 ~ s(log.bathym, bs = "cr", k = 5, m = 2) +
                       s(log.npp, bs = "cr", k = 5, m = 2) +
                       s(log.sst, bs = "cr", k = 5, m = 2) +
                       s(id, bs = "re"), data = rsf.pts_10s2, method = "fREML",
                     family = poisson(), weights = wts2, discrete = TRUE)
toc()  #took 1 sec to run

summary(fit.GAM10)
plot(fit.GAM10, scale = 0, shade = TRUE, shade.col = "lightblue")
gam.check(fit.GAM10)
concurvity(fit.GAM10, full = TRUE)
concurvity(fit.GAM10, full = FALSE)  #high concurvity (> 0.8) related to varying intercept; fine moving forward


# Run full model
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
gam.check(fit.HGAM10)



# Create custom partial effects plot
# evaluate the smooths
sm <- smooth_estimates(fit.HGAM10_PI) %>%
  add_confint()
sm

# add partial residuals to data
rsf.pts_10s3 <- rsf.pts_10s2 %>%
  add_partial_residuals(fit.HGAM10_PI)


p.bathym.pop <- sm %>%
  filter(smooth == "s(log.bathym)") %>%
  ggplot() +
  geom_rug(aes(x = exp(log.bathym)),
           data = rsf.pts_10s3,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = exp(log.bathym)),
              alpha = 0.5, fill = "steelblue3") +
  # geom_point(aes(x = exp(log.bathym), y = `s(log.bathym)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "steelblue3", alpha = 0.2) +
  geom_line(aes(x = exp(log.bathym), y = est), linewidth = 1.2) +

  labs(x = 'Depth (m)', y = "Partial Effect", title = "s(log.bathym)") +
  theme_bw() +
  theme(legend.position = "none")


p.bathym.id <- sm %>%
  drop_na(log.bathym, id) %>%
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.bathym),
                y = est + sm[sm$smooth == 's(log.bathym)',]$est, color = id),
            linewidth = 0.5) +
  labs(x = 'Depth (m)', y = "Partial Effect", title = "s(log.bathym)") +
  lims(x = c(0,250), y = c(-600, 600)) +
  theme_bw() +
  theme(legend.position = "none")




p.npp.pop <- sm %>%
  filter(smooth == "s(log.npp)") %>%
  ggplot() +
  geom_rug(aes(x = exp(log.npp) / 1000),
           data = rsf.pts_10s3,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = exp(log.npp) / 1000),
              alpha = 0.5, fill = 'darkgreen') +
  # geom_point(aes(x = exp(log.npp) / 1000, y = `s(log.npp)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "darkgreen", alpha = 0.2) +
  geom_line(aes(x = exp(log.npp) / 1000, y = est), linewidth = 1.2) +
  labs(x = 'Net Primary Productivity (g C m-2 day-1)', y = "Partial Effect", title = "s(log.npp)") +
  lims(x = c(0,30)) +
  theme_bw() +
  theme(legend.position = "none")


p.npp.id <- sm %>%
  drop_na(log.npp, id) %>%
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.npp) / 1000,
                y = est + sm[sm$smooth == 's(log.npp)',]$est, color = id),
            linewidth = 0.5) +
  labs(x = 'Net Primary Productivity (g C m-2 day-1)', y = "Partial Effect", title = "s(log.npp)") +
  lims(x = c(0,30), y = c(-6000, 2500)) +
  theme_bw() +
  theme(legend.position = "none")




p.sst.pop <- sm %>%
  filter(smooth == "s(log.sst)") %>%
  ggplot() +
  geom_rug(aes(x = exp(log.sst)),
           data = rsf.pts_10s3,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = exp(log.sst)),
              alpha = 0.5, fill = 'firebrick') +
  # geom_point(aes(x = exp(log.sst), y = `s(log.sst)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "darkgreen", alpha = 0.2) +
  geom_line(aes(x = exp(log.sst), y = est), linewidth = 1.2) +
  labs(x = 'Sea Surface Temperature (°C)', y = "Partial Effect", title = "s(log.sst)") +
  theme_bw() +
  theme(legend.position = "none")


p.sst.id <- sm %>%
  drop_na(log.sst, id) %>%
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.sst),
                y = est + sm[sm$smooth == 's(log.sst)',]$est, color = id),
            linewidth = 0.5) +
  labs(x = 'Sea Surface Temperature (°C)', y = "Partial Effect", title = "s(log.sst)") +
  theme_bw() +
  theme(legend.position = "none")



p.bathym.pop + p.bathym.id
p.npp.pop + p.npp.id
p.sst.pop + p.sst.id
