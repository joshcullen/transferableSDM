
### Fit RSF as GLMM w/ Gaussian Process Prior on SST ###

library(tidyverse)
library(brms)
library(bayesplot)

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)


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
# obs.ind_10 <- which(rsf.pts_10s$obs == 1)
# rsf.pts_10s <- rsf.pts_10s %>%
#   mutate(bathym.s = (bathym - mean(bathym[obs.ind_10])) / sd(bathym),
#          k490.s = (k490 - mean(k490[obs.ind_10])) / sd(k490),
#          npp.s = (npp - mean(npp[obs.ind_10])) / sd(npp),
#          sst.s = (sst - mean(sst[obs.ind_10])) / sd(sst))
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(bathym.s = as.numeric(scale(bathym)),
         k490.s = as.numeric(scale(k490)),
         npp.s = as.numeric(scale(npp)),
         sst.s = as.numeric(scale(sst))
         )

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




## Convert ID to factor for {brms}
rsf.pts_10s$id1 <- factor(rsf.pts_10s$id)
# rsf.pts_10s$id2 <- rsf.pts_10s$id1

# rsf.pts_10s <- arrange(rsf.pts_10s, id1)


# rsf.pts_30s$id1 <- as.numeric(factor(rsf.pts_30s$id))
# rsf.pts_30s <- arrange(rsf.pts_30s, id1)
#
#
# rsf.pts_50s$id1 <- as.numeric(factor(rsf.pts_50s$id))
# rsf.pts_50s <- arrange(rsf.pts_50s, id1)


tmp <- rsf.pts_10s %>%
  slice_sample(n = 5000) %>%
  arrange(id1)


# Fit simple GLMM; coeffs estimated from MVN
fit1 <- brm(obs/wts2 | weights(wts2) ~ log.bathym + log.npp + log.sst + (log.bathym + log.npp + log.sst|id1),
            data = tmp,
            iter = 500, warmup = 250, refresh = 100,
            chains = 4, cores = 4, family = "poisson", seed = 2023, backend = "cmdstanr")
# took 4 min

summary(fit1)
plot(fit1)

me1 <- conditional_effects(fit1, ndraws = 200, spaghetti = TRUE)
plot(me1, ask = TRUE, points = FALSE)





tmp2 <- rsf.pts_10s %>%
  filter(id %in% c(181800, 181796, 159776))

# Fit GP RSF w/ fixed effects
## NEED TO SET SENSIBLE PRIORS!!!
get_prior(obs | weights(wts) ~ bathym.s + npp.s + sst.s +
            gp(bathym.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
            gp(npp.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
            gp(sst.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
            (1 | id1), #+
            # (1 + bathym.s || id1),
          data = tmp2,
          family = "bernoulli")


fit.gp <- brm(obs | weights(wts) ~ #bathym.s + npp.s + sst.s +
                gp(bathym.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
                gp(npp.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
                gp(sst.s, by = id1, gr = T, scale = FALSE, k = 5, c = 5/4) +
                (1 | id1),
                # (1 + bathym.s || id1),
              prior = c(prior(normal(0,10), class = Intercept),
                        # prior(normal(0,10), class = b),
                        prior(constant(1000), class = sd, coef = Intercept, group = id1),
                        # prior(inv_gamma(2.874624, 2.941204), class = lscale),
                        prior(exponential(1), class = sdgp)#,
                        # prior(exponential(1), class = sd, coef = bathym.s, group = id1)
                        ),
              data = tmp2,
              iter = 2000, warmup = 1000, refresh = 100,
              chains = 3, cores = 3, family = "bernoulli", seed = 2023,
              backend = "cmdstanr", threads = threading(3),
              control = list(adapt_delta = 0.9, max_treedepth = 15))
# took 1 hr

summary(fit.gp)
plot(fit.gp)
rhat_vals <- rhat(fit.gp)
mcmc_rhat_data(rhat_vals)
mcmc_rhat(rhat_vals) + theme_bw()

me1 <- conditional_effects(fit.gp, ndraws = 200, spaghetti = TRUE)
plot(me1, ask = TRUE, points = FALSE)

# post_ranef <- as_draws(fit.gp_thread, add_chain = T)
# mcmc_trace(post_ranef)
