
### Fit RSF as GLMM ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(sf)
library(tictoc)
library(lubridate)
library(BRRR)

source('Scripts/helper functions.R')



#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv")
rsf.pts_30 <- read_csv("Processed_data/GoM_Cm_RSFprep_30x.csv")
rsf.pts_50 <- read_csv("Processed_data/GoM_Cm_RSFprep_50x.csv")

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



####################
### Fit GLMM RSF ###
####################

# Center and scale covars; remove rows w/ incomplete observations
rsf.pts_10s <- rsf.pts_10 %>%
  drop_na(bathym, npp, sst)
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
  drop_na(bathym, npp, sst)
# obs.ind_30 <- which(rsf.pts_30s$obs == 1)
# rsf.pts_30s <- rsf.pts_30s %>%
#   mutate(bathym.s = (bathym - mean(bathym[obs.ind_30])) / sd(bathym),
#          k490.s = (k490 - mean(k490[obs.ind_30])) / sd(k490),
#          npp.s = (npp - mean(npp[obs.ind_30])) / sd(npp),
#          sst.s = (sst - mean(sst[obs.ind_30])) / sd(sst))

rsf.pts_50s <- rsf.pts_50 %>%
  drop_na(bathym, npp, sst)
# obs.ind_50 <- which(rsf.pts_50s$obs == 1)
# rsf.pts_50s <- rsf.pts_50s %>%
#   mutate(bathym.s = (bathym - mean(bathym[obs.ind_50])) / sd(bathym),
#          k490.s = (k490 - mean(k490[obs.ind_50])) / sd(k490),
#          npp.s = (npp - mean(npp[obs.ind_50])) / sd(npp),
#          sst.s = (sst - mean(sst[obs.ind_50])) / sd(sst))



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
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_30s <- rsf.pts_30s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))

rsf.pts_50s <- rsf.pts_50s %>%
  mutate(log.bathym = log(abs(bathym)),
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
# A <- as.numeric(st_area(bbox_mask) / 1e6)  #in km^2
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)
rsf.pts_30s$wts2 <- ifelse(rsf.pts_30s$obs == 0, A / sum(rsf.pts_30s$obs == 0), 1e-6)
rsf.pts_50s$wts2 <- ifelse(rsf.pts_50s$obs == 0, A / sum(rsf.pts_50s$obs == 0), 1e-6)
#
# rsf.pts_10s2 <- rsf.pts_10s %>%
#   nest(data = -id) %>%
#   mutate(glm = map(data, ~glm(obs ~ sqrt.bathym + I(sqrt.bathym^2) + sqrt.npp + I(sqrt.npp^2) + sqrt.sst + I(sqrt.sst^2),
#                               data = .x, family = binomial(), weights = wts))
#   )
# rsf.pts_50s3 <- rsf.pts_50s2 %>%
#   mutate(coeffs = map(glm, ~coef(.x)[-1])) %>%
#   dplyr::select(coeffs) %>%
#   unnest(cols = coeffs) %>%
#   bind_rows()

logis.mod <- glm(obs ~ log.bathym + log.npp + log.sst,
                data = rsf.pts_10s, family = binomial(), weights = wts)

summary(logis.mod)
car::vif(logis.mod)  #VIF looks good; low
#
#
# # Down-weighted Poisson regression
# A <- (res(cov_list$bathym)[1] / 1000) ^ 2  #area of grid cells
# rsf.pts$wts2 <- ifelse(rsf.pts$obs == 0, A/sum(rsf.pts$obs == 0), 1e-6)
#
# pois.mod <- glm(obs/wts2 ~ bathym + k490 + npp + sst,
#                 data = rsf.pts, family = poisson(), weights = wts2)
#
# summary(pois.mod)





## Mixed RSF via INLA
rsf.pts_10s$id1 <- as.numeric(factor(rsf.pts_10s$id))
rsf.pts_10s$id2 <- rsf.pts_10s$id1
rsf.pts_10s$id3 <- rsf.pts_10s$id1
rsf.pts_10s$id4 <- rsf.pts_10s$id1
rsf.pts_10s$id5 <- rsf.pts_10s$id1
rsf.pts_10s$id6 <- rsf.pts_10s$id1
rsf.pts_10s$id7 <- rsf.pts_10s$id1

rsf.pts_10s <- arrange(rsf.pts_10s, id1)


rsf.pts_30s$id1 <- as.numeric(factor(rsf.pts_30s$id))
rsf.pts_30s$id2 <- rsf.pts_30s$id1
rsf.pts_30s$id3 <- rsf.pts_30s$id1
rsf.pts_30s$id4 <- rsf.pts_30s$id1
rsf.pts_30s$id5 <- rsf.pts_30s$id1
rsf.pts_30s$id6 <- rsf.pts_30s$id1
rsf.pts_30s$id7 <- rsf.pts_30s$id1

rsf.pts_30s <- arrange(rsf.pts_30s, id1)


rsf.pts_50s$id1 <- as.numeric(factor(rsf.pts_50s$id))
rsf.pts_50s$id2 <- rsf.pts_50s$id1
rsf.pts_50s$id3 <- rsf.pts_50s$id1
rsf.pts_50s$id4 <- rsf.pts_50s$id1
rsf.pts_50s$id5 <- rsf.pts_50s$id1
rsf.pts_50s$id6 <- rsf.pts_50s$id1
rsf.pts_50s$id7 <- rsf.pts_50s$id1

rsf.pts_50s <- arrange(rsf.pts_50s, id1)

# create vector of ID values
id.vals <- unique(rsf.pts_10s$id1)

RSF.formula <- obs/wts2 ~ log.bathym + I(log.bathym ^ 2) + log.npp + I(log.npp ^ 2) + log.sst + I(log.sst ^ 2) +
  f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id2, log.bathym, values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id3, I(log.bathym ^ 2), values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id4, log.npp, values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id5, I(log.npp ^ 2), values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id6, log.sst, values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  f(id7, I(log.sst ^ 2), values = id.vals, model = "iid",
    hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))
# RSF.formula <- obs/wts2 ~ bathym.s + I(bathym.s ^ 2) + npp.s + I(npp.s ^ 2) + sst.s + I(sst.s ^ 2) +
#   f(id1, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
#   f(id2, bathym.s, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id3, I(bathym.s ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id4, npp.s, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id5, I(npp.s ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id6, sst.s, values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
#   f(id7, I(sst.s ^ 2), values = id.vals, model = "iid",
#     hyper = list(theta = list(fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))


# Fit the model
set.seed(2023)
tic()
fit.RSF_10 <- inla(RSF.formula, family = "Poisson", data = rsf.pts_10s, weights = rsf.pts_10s$wts2,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = 1e-3)),
                   control.compute = list(waic = TRUE,
                                          dic = TRUE), verbose = FALSE
)
toc()  # took 45 sec to run for x10; 76 sec on laptop

summary(fit.RSF_10)


set.seed(2023)
tic()
fit.RSF_30 <- inla(RSF.formula, family = "Poisson", data = rsf.pts_30s, weights = rsf.pts_30s$wts2,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = 1e-3)),
                   control.compute = list(waic = TRUE,
                                          dic = TRUE), verbose = FALSE
)
toc()  # took 3 min to run for 30x; 6.5 min on laptop

summary(fit.RSF_30)


set.seed(2023)
tic()
fit.RSF_50 <- inla(RSF.formula, family = "Poisson", data = rsf.pts_50s, weights = rsf.pts_50s$wts2,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = 1e-3)),
                   control.compute = list(waic = TRUE,
                                          dic = TRUE), verbose = FALSE
)
toc()  # took 4.5 min to run for 50x; 23 min on laptop

summary(fit.RSF_50)



# The summary for the posterior distribution of the fixed effects:
# fixed.coeffs <- cbind(fit.RSF_10$summary.fixed[-1,],
#                       sigma = inla_emarginal(fit.RSF_10) %>%
#                         sqrt()) %>%
#   mutate(lo = mean - sigma, hi = mean + sigma) %>%
#   mutate(param = factor(rownames(.), levels = rownames(.)))
# fixed.coeffs <- fit.RSF_10$summary.fixed[-1,] %>%
#   mutate(param = factor(rownames(.), levels = rownames(.)))


fixed.coeffs.list <- list(x10 = cbind(fit.RSF_10$summary.fixed[-1,],
                                      sigma = inla_emarginal(fit.RSF_10) %>%
                                        sqrt()),
                          x30 = cbind(fit.RSF_30$summary.fixed[-1,],
                                      sigma = inla_emarginal(fit.RSF_30) %>%
                                        sqrt()),
                          x50 = cbind(fit.RSF_50$summary.fixed[-1,],
                                      sigma = inla_emarginal(fit.RSF_50) %>%
                                        sqrt()))
fixed.coeffs <- fixed.coeffs.list %>%
  map(~{.x %>%
      mutate(lo = mean - sigma, hi = mean + sigma) %>%
      mutate(param = factor(rownames(.x), levels = rownames(.x)))}) %>%
  bind_rows(.id = "dataset")
# dplyr::slice(2:n()) %>%  #remove intercept


ggplot(fixed.coeffs, aes(x = param)) +
  geom_hline(yintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(ymin = lo, ymax = hi, color = dataset), position = position_dodge(width = 0.1)) +
  geom_point(aes(y = mean, color = dataset), position = position_dodge(width = 0.1)) +
  theme_bw()



random.coeffs.list <- list(x10 = fit.RSF_10$summary.random[-1] %>%
                             bind_rows(.id = "param") %>%
                             mutate(across(param, factor)),
                           x30 = fit.RSF_30$summary.random[-1] %>%
                             bind_rows(.id = "param") %>%
                             mutate(across(param, factor)),
                           x50 = fit.RSF_50$summary.random[-1] %>%
                             bind_rows(.id = "param") %>%
                             mutate(across(param, factor)))

random.coeffs <- random.coeffs.list %>%
  bind_rows(.id = "dataset")
levels(random.coeffs$param) <- unique(fixed.coeffs$param)

# random.coeffs <- fit.RSF_10$summary.random[-1] %>%
#   bind_rows(.id = "param") %>%
#   mutate(across(param, factor))
# levels(random.coeffs$param) <- fixed.coeffs$param

# Add population means to random effects
random.coeffs <- random.coeffs %>%
  group_by(dataset, param) %>%
  group_split() %>%
  map2(.x = .,
       .y = fixed.coeffs %>%
         group_by(dataset, param) %>%
         group_split(),
       ~mutate(.x,
               mean = mean + .y$mean,
               `0.025quant` = `0.025quant` + .y$mean,
               `0.975quant` = `0.975quant` + .y$mean)) %>%
  bind_rows()
# random.coeffs <- rbind(random.coeffs,
#                        fixed.coeffs %>%
#                          mutate(ID = "Pop", .before = mean) %>%
#                          relocate(param, .before = ID)) %>%
#   arrange(param)

ggplot(random.coeffs) +
  geom_hline(yintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, color = param)) +
  geom_point(aes(y = mean, x = ID, color = param)) +
  geom_linerange(data = fixed.coeffs, aes(ymin = lo, ymax = hi), x = 50) +
  geom_point(data = fixed.coeffs, aes(y = mean), x = 50) +
  theme_bw() +
  facet_grid(dataset ~ param, scales = "free")




####################################################
### Viz marginal effects plots per environ covar ###
####################################################

# fixed.coeffs2 <- fixed.coeffs[,c("mean","0.025quant","0.975quant")] %>%
#   as.matrix()

random.coeffs2 <- random.coeffs %>%
  filter(dataset == 'x10') %>%
  dplyr::select(param, ID, mean, `0.025quant`, `0.975quant`)


### Bathymetry ###

bathym.newdata <- data.frame(bathym = seq(min(rsf.pts_10s$log.bathym), max(rsf.pts_10s$log.bathym), length.out = 500),
                             bathym.2 = seq(min(rsf.pts_10s$log.bathym), max(rsf.pts_10s$log.bathym), length.out = 500) ^ 2,
                             npp = 0,
                             npp.2 = 0,
                             sst = 0,
                             sst.2 = 0) %>%
  as.matrix()


# Generate samples for marginal effects from posterior
# set.seed(2023)
# foo <- fixed.coeffs %>%
#   split(.$dataset) %>%
#   pluck(1) %>%
#   filter(str_detect(param, "bathym")) %>%
#   dplyr::select(mean, sigma) %>%
#   pmap(., function(mean, sigma) rnorm(1000, mean = mean, sd = sigma)) %>%
#   map(data.frame) %>%
#   bind_cols() %>%
#   mutate(npp = 0, npp.2 = 0, sst = 0, sst.2 = 0) %>%
#   t()
#
# tmp <- bathym.newdata %*% foo %>%
#   data.frame() %>%
#   mutate(bathym = exp(bathym.newdata[,1]))




pred.bathym <- vector("list", length = n_distinct(random.coeffs2$ID))
names(pred.bathym) <- unique(random.coeffs2$ID)

for (i in 1:n_distinct(random.coeffs2$ID)) {

  coeff1 <- random.coeffs2 %>%
    filter(ID == unique(random.coeffs2$ID)[i]) %>%
    dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
    as.matrix()

  tmp <- bathym.newdata %*% coeff1 %>%
    data.frame() %>%
    mutate(bathym = exp(bathym.newdata[,1]))

  pred.bathym[[i]] <- tmp
}
pred.bathym <- pred.bathym %>%
  bind_rows(.id = "id")

tmp <- fixed.coeffs %>%
  filter(dataset == "x10") %>%
  dplyr::select(mean,`0.025quant`,`0.975quant`) %>%
  as.matrix()

pred.bathym.pop <- bathym.newdata %*% tmp %>%
  data.frame() %>%
  mutate(bathym = exp(bathym.newdata[,1]))



# Pop mean in black; ID by color
ggplot() +
  # geom_ribbon(data = pred.bathym.pop, aes(x = bathym, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.bathym.pop, aes(x = bathym, y = exp(mean)), linewidth = 1) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

ggplot() +
  geom_line(data = pred.bathym, aes(x = bathym, y = exp(mean), group = id, color = id),
            linewidth = 1, show.legend = FALSE) +
  theme_bw() +
  lims(x = c(0,300), y = c(0,1500)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())


# ggplot() +
#   geom_density(data = rsf.pts, aes(bathym, fill = factor(obs)), alpha = 0.5) +
#   xlim(-1000,0) +
#   theme_bw()
# rsf.pts %>%
#   filter(obs == 0) %>%
#   pull(bathym) %>%
#   summary()





### SST ###

sst.newdata <- data.frame(bathym = 0,
                          bathym.2 = 0,
                          npp = 0,
                          npp.2 = 0,
                          sst = seq(min(rsf.pts_10s$log.sst), max(rsf.pts_10s$log.sst), length.out = 500),
                          sst.2 = seq(min(rsf.pts_10s$log.sst), max(rsf.pts_10s$log.sst), length.out = 500) ^ 2) %>%
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
    # mutate(sst = (sst.newdata[,7] * sd(rsf.pts_10s$sst, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$sst, na.rm = T))
    mutate(sst = exp(sst.newdata[,5]))

  pred.sst[[i]] <- tmp
}

pred.sst <- pred.sst %>%
  bind_rows(.id = "id")


tmp <- fixed.coeffs %>%
  filter(dataset == "x10") %>%
  dplyr::select(mean) %>%
  as.matrix()

pred.sst.pop <- sst.newdata %*% tmp %>%
  data.frame() %>%
  # mutate(sst = (sst.newdata[,1] * sd(rsf.pts_10s$sst, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$sst, na.rm = T))
  mutate(sst = exp(sst.newdata[,5]))



# Pop mean in black; ID by color
ggplot() +
  # geom_ribbon(data = pred.sst.pop), aes(x = sst, ymin = plogis(X0.025quant), ymax = plogis(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.sst.pop, aes(x = sst, y = exp(mean)), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

ggplot() +
  geom_line(data = pred.sst, aes(x = sst, y = exp(mean), group = id, color = id),
            linewidth = 1, show.legend = FALSE) +
  theme_bw() +
  lims(y = c(0,2e140)) +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())





### NPP ###

npp.newdata <- data.frame(bathym = 0,
                          bathym.2 = 0,
                          npp = seq(min(rsf.pts_10s$log.npp), max(rsf.pts_10s$log.npp), length.out = 500),
                          npp.2 = seq(min(rsf.pts_10s$log.npp), max(rsf.pts_10s$log.npp), length.out = 500) ^ 2,
                          sst = 0,
                          sst.2 = 0) %>%
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
    # mutate(npp = (npp.newdata[,7] * sd(rsf.pts_10s$npp, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$npp, na.rm = T))
    mutate(npp = exp(npp.newdata[,3]))

  pred.npp[[i]] <- tmp
}

pred.npp <- pred.npp %>%
  bind_rows(.id = "id")


tmp <- fixed.coeffs %>%
  filter(dataset == "x10") %>%
  dplyr::select(mean) %>%
  as.matrix()

pred.npp.pop <- npp.newdata %*% tmp %>%
  data.frame() %>%
  # mutate(npp = (npp.newdata[,1] * sd(rsf.pts_10s$npp, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$npp, na.rm = T))
  mutate(npp = exp(npp.newdata[,3]))

# Pop mean in black; ID by color
ggplot() +
  # geom_ribbon(data = pred.npp.pop, aes(x = npp, ymin = exp(X0.025quant), ymax = exp(X0.975quant)), alpha = 0.4) +
  geom_line(data = pred.npp.pop, aes(x = npp / 1000, y = exp(mean)), linewidth = 1.5) +
  theme_bw() +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))

ggplot() +
  geom_line(data = pred.npp, aes(x = npp / 1000, y = exp(mean), group = id, color = id),
            linewidth = 1, show.legend = FALSE) +
  theme_bw() +
  lims(y = c(0,1000)) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())






#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
files <- files[!grepl(pattern = "Kd490", files)]  #remove Kd490 datasets
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list <- sapply(files, rast)
cov_list

names(cov_list) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list[[var]]) <- gsub(names(cov_list[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list[["bathym"]][cov_list[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')




####################################
### Generate predictive surfaces ###
####################################

# For pop mean in 2020-09
newdat <- data.frame(log.bathym = terra::values(cov_list$bathym) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = terra::values(cov_list$bathym) %>%
                       abs() %>%
                       log() %>%
                       apply(., 2, function(x) x^2),
                     log.npp = terra::values(cov_list$npp$`2020-09-01`) %>%
                       log(),
                     log.npp2 = terra::values(cov_list$npp$`2020-09-01`) %>%
                       log() %>%
                       apply(., 2, function(x) x^2),
                     log.sst = terra::values(cov_list$sst$`2020-09-01`) %>%
                       log(),
                     log.sst2 = terra::values(cov_list$sst$`2020-09-01`) %>%
                       log() %>%
                       apply(., 2, function(x) x^2))
names(newdat) <- c('log.bathym','log.bathym2','log.npp','log.npp2','log.sst','log.sst2')
summary(newdat)

coeff1 <- fit.RSF_10$summary.fixed$mean[-1]
mean.pred <- as.matrix(newdat) %*% coeff1  #make predictions


rast.pred <- cov_list$bathym
terra::values(rast.pred) <- exp(mean.pred)
plot(rast.pred)


rast.pred.df <- as.data.frame(rast.pred, xy = TRUE)
names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)

ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = log(pred))) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  labs(x="",y="", title = "Population Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20))






# For ID 181796 in 2020-08
rsf.pts_10s %>%
  filter(id == 181796, obs == 1) %>%
  group_by(month.year) %>%
  count()  # most obs are in August

newdat <- data.frame(log.bathym = terra::values(cov_list$bathym) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = terra::values(cov_list$bathym) %>%
                       abs() %>%
                       log() %>%
                       apply(., 2, function(x) x^2),
                     log.npp = terra::values(cov_list$npp$`2020-08-01`) %>%
                       log(),
                     log.npp2 = terra::values(cov_list$npp$`2020-08-01`) %>%
                       log() %>%
                       apply(., 2, function(x) x^2),
                     log.sst = terra::values(cov_list$sst$`2020-08-01`) %>%
                       log(),
                     log.sst2 = terra::values(cov_list$sst$`2020-08-01`) %>%
                       log() %>%
                       apply(., 2, function(x) x^2))
names(newdat) <- c('log.bathym','log.bathym2','log.npp','log.npp2','log.sst','log.sst2')
summary(newdat)

# Specify terms for ID 181796
ind <- which(unique(rsf.pts_10s$id) == 181796)
coeff.181796 <- random.coeffs2 %>%
  filter(ID == ind) %>%
  pull(mean)
x181796.pred <- as.matrix(newdat) %*% coeff.181796  #make predictions


rast.pred.181796 <- cov_list$bathym
terra::values(rast.pred.181796) <- exp(x181796.pred)


rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
names(rast.pred.181796.df)[3] <- "pred"

ggplot() +
  geom_raster(data = rast.pred.181796.df, aes(x, y, fill = log(pred))) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_path(data = rsf.pts_10s %>%
              filter(id == 181796, month.year == '2020-08-01', obs == 1), aes(x, y),
            color = "chartreuse", linewidth = 0.5) +
  labs(x="",y="", title = "ID 181796: August 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20))









#####################
### Validate HGAM ###
#####################

dat.br <- read_csv("Processed_data/Brazil_Cm_Tracks_behav.csv")
dat.qa <- read_csv("Processed_data/Qatar_Cm_Tracks_behav.csv")

# Create indexing column "month.year"
dat.br <- dat.br %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))

dat.qa <- dat.qa %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01"))






### Brazil ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Brazil", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_br <- sapply(files, rast)
# cov_list_br

names(cov_list_br) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_br[[var]]) <- gsub(names(cov_list_br[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list_br[[var]] <- terra::resample(cov_list_br[[var]], cov_list_br$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list_br <- map(cov_list_br, terra::project, 'EPSG:3395')






### Validate model using CBI ###

my.ind.br <- names(cov_list_br$npp)
br.rast.pred <- rep(cov_list_br$bathym, nlyr(cov_list_br$npp))
names(br.rast.pred) <- my.ind.br

# Define coeff values from GLMM
coeff1 <- fit.RSF_10$summary.fixed$mean[-1]

tic()
for (i in 1:nlyr(cov_list_br$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = as.vector(terra::values(cov_list_br$bathym)) %>%
                       abs() %>%
                       log() %>%
                       sapply(., function(x) x^2),
                     log.npp = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log(),
                     log.npp2 = as.vector(terra::values(cov_list_br$npp[[my.ind.br[i]]])) %>%
                       log() %>%
                       sapply(.,function(x) x^2),
                     log.sst = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log(),
                     log.sst2 = as.vector(terra::values(cov_list_br$sst[[my.ind.br[i]]])) %>%
                       log() %>%
                       sapply(., function(x) x^2))


  # Make predictions on intensity of use from model
  br.pred <- as.matrix(vars) %*% coeff1  #make predictions

  # Store results in raster stack
  terra::values(br.rast.pred[[i]]) <- br.pred  #keep on log-scale since response scale results in crazy large values

}
skrrrahh('khaled2')
toc()  #took 1 min



cbi.br <- vector("list", nlyr(br.rast.pred))
tic()
for (i in 1:nlyr(br.rast.pred)) {

  # Subset tracks by month.year
  obs <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  cbi.br[[i]] <- cbi(fit = as(br.rast.pred[[i]], "Raster"),
                     obs = obs,
                     nclass = 10,
                     PEplot = FALSE)
}
skrrrahh("khaled3")
toc()  #took 34 sec


cbi.br <- cbi.br %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Brazil")








### Qatar ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Qatar", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_qa <- sapply(files, rast)
cov_list_qa

names(cov_list_qa) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_qa[[var]]) <- gsub(names(cov_list_qa[[var]]), pattern = "-..$", replacement = "-01")
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list_qa[[var]] <- terra::resample(cov_list_qa[[var]], cov_list_qa$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list_qa <- map(cov_list_qa, terra::project, 'EPSG:3395')






### Validate model using CBI ###

my.ind.qa <- names(cov_list_qa$npp)
qa.rast.pred <- rep(cov_list_qa$bathym, nlyr(cov_list_qa$npp))
names(qa.rast.pred) <- my.ind.qa

# Define coeff values from GLMM
coeff1 <- fit.RSF_10$summary.fixed$mean[-1]

tic()
for (i in 1:nlyr(cov_list_qa$npp)) {

  # Subset covars by month.year
  vars <- data.frame(log.bathym = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log(),
                     log.bathym2 = as.vector(terra::values(cov_list_qa$bathym)) %>%
                       abs() %>%
                       log() %>%
                       sapply(., function(x) x^2),
                     log.npp = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log(),
                     log.npp2 = as.vector(terra::values(cov_list_qa$npp[[my.ind.qa[i]]])) %>%
                       log() %>%
                       sapply(.,function(x) x^2),
                     log.sst = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log(),
                     log.sst2 = as.vector(terra::values(cov_list_qa$sst[[my.ind.qa[i]]])) %>%
                       log() %>%
                       sapply(., function(x) x^2))


  # Make predictions on intensity of use from model
  qa.pred <- as.matrix(vars) %*% coeff1  #make predictions

  # Store results in raster stack
  terra::values(qa.rast.pred[[i]]) <- qa.pred

}
skrrrahh('khaled2')
toc()  #took 1 sec



cbi.qa <- vector("list", nlyr(qa.rast.pred))
tic()
for (i in 1:nlyr(qa.rast.pred)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  cbi.qa[[i]] <- cbi(fit = as(qa.rast.pred[[i]], "Raster"),
                     obs = obs,
                     nclass = 10,
                     PEplot = FALSE)
}
skrrrahh("khaled3")
toc()  #took 6 sec


cbi.qa <- cbi.qa %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(cor = .,
             Region = "Qatar")






### Summarize Validation Results ###

cbi.fit <- rbind(cbi.br, cbi.qa)

cbi.mean <- cbi.fit %>%
  group_by(Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE))

set.seed(2023)
ggplot(data = cbi.fit, aes(Region, cor)) +
  # geom_boxplot(aes(fill = Region), outlier.color = NA) +
  geom_jitter(aes(fill = Region), pch = 21, height = 0, width = 0.2, alpha = 0.7, size = 5) +
  scale_fill_manual(values = c("#1DB100", "#00A2FF"), guide = "none") +
  geom_point(data = cbi.mean, aes(Region, mean), fill = "black", pch = 21, size = 6) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Continuous Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24))
