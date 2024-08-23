
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

gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet")

glimpse(rsf.pts_10)
summary(rsf.pts_10)



####################
### Fit HGLM RSF ###
####################

# Remove rows w/ incomplete observations
rsf.pts_10s <- rsf.pts_10 %>%
  drop_na(bathym, npp, sst)


# Explore used vs available habitat values
rsf.pts_10s %>%
  pivot_longer(cols = c(bathym, npp, sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")


# Log-transform skewed covars to allow model fitting
rsf.pts_10s <- rsf.pts_10s %>%
  mutate(log.bathym = log(abs(bathym)),
         log.npp = log(npp),
         log.sst = log(sst))


# Check Pearson corrs
cor(rsf.pts_10s[,c('bathym','npp','sst')])  #All low (< 0.17)


# Now explore transformed distributions
rsf.pts_10s %>%
  pivot_longer(cols = c(log.bathym, log.npp, log.sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free")


# Down-weighted Poisson regression
A <- 4759.836 ^ 2  #in m^2; pixel res is 4759.836 m
rsf.pts_10s$wts <- ifelse(rsf.pts_10s$obs == 0, 5000, 1)
rsf.pts_10s$wts2 <- ifelse(rsf.pts_10s$obs == 0, A / sum(rsf.pts_10s$obs == 0), 1e-6)


logis.mod <- glm(obs ~ log.bathym + log.npp + log.sst,
                data = rsf.pts_10s, family = binomial(), weights = wts)

summary(logis.mod)

pois.mod <- glm(obs/wts2 ~ log.bathym + log.npp + log.sst,
                 data = rsf.pts_10s, family = poisson(), weights = wts2)

summary(pois.mod)
car::vif(pois.mod)  #VIF looks good; low (< 3)






## HGLM RSF via INLA
rsf.pts_10s$id1 <- as.numeric(factor(rsf.pts_10s$id))
rsf.pts_10s$id2 <- rsf.pts_10s$id1
rsf.pts_10s$id3 <- rsf.pts_10s$id1
rsf.pts_10s$id4 <- rsf.pts_10s$id1
rsf.pts_10s$id5 <- rsf.pts_10s$id1
rsf.pts_10s$id6 <- rsf.pts_10s$id1
rsf.pts_10s$id7 <- rsf.pts_10s$id1

rsf.pts_10s <- arrange(rsf.pts_10s, id1)




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


# Fit correlative HGLM
set.seed(2023)
tic()
fit.RSF_corr <- inla(RSF.formula, family = "Poisson", data = rsf.pts_10s, weights = rsf.pts_10s$wts2,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = 1e-3)),
                   control.compute = list(waic = TRUE,
                                          dic = TRUE)#, num.threads = 1:1
)
toc()  # took 25 sec; took 50 sec to run w/ threads = 1

summary(fit.RSF_corr)


# Fit hybrid HGLM
set.seed(2023)
tic()
fit.RSF_hybrid <- inla(RSF.formula, family = "Poisson", data = rsf.pts_10s, weights = rsf.pts_10s$wts2,
                     control.fixed = list(
                       mean = list(log.sst = 30*6.592, `I(log.sst^2)` = 31*-1),
                       prec = list(log.sst = 0.005, `I(log.sst^2)` = 0.005)
                       ),
                     control.compute = list(waic = TRUE,
                                            dic = TRUE)#, num.threads = 1:1
)
toc()  # took 25 sec; took 50 sec to run w/ threads = 1

summary(fit.RSF_hybrid)



# The summary for the posterior distribution of the fixed effects:
# fixed.coeffs <- cbind(fit.RSF_10$summary.fixed[-1,],
#                       sigma = inla_emarginal(fit.RSF_10) %>%
#                         sqrt()) %>%
#   mutate(lo = mean - sigma, hi = mean + sigma) %>%
#   mutate(param = factor(rownames(.), levels = rownames(.)))
# fixed.coeffs <- fit.RSF_10$summary.fixed[-1,] %>%
#   mutate(param = factor(rownames(.), levels = rownames(.)))


fixed.coeffs.list <- list(corr = fit.RSF_corr$summary.fixed[-1,],
                          hybrid = fit.RSF_hybrid$summary.fixed[-1,])
fixed.coeffs <- fixed.coeffs.list %>%
  map(~{.x %>%
      # mutate(lo = mean - sigma, hi = mean + sigma) %>%
      mutate(param = factor(rownames(.x), levels = rownames(.x)))}) %>%
  bind_rows(.id = "dataset")
# dplyr::slice(2:n()) %>%  #remove intercept


ggplot(fixed.coeffs, aes(x = param)) +
  geom_hline(yintercept = 0, linewidth = 0.75) +
  geom_linerange(aes(ymin = `0.025quant`, ymax = `0.975quant`, color = dataset), position = position_dodge(width = 0.1)) +
  geom_point(aes(y = mean, color = dataset), position = position_dodge(width = 0.1)) +
  theme_bw()



# random.coeffs.list <- list(x10 = fit.RSF_10$summary.random[-1] %>%
#                              bind_rows(.id = "param") %>%
#                              mutate(across(param, factor)),
#                            x30 = fit.RSF_30$summary.random[-1] %>%
#                              bind_rows(.id = "param") %>%
#                              mutate(across(param, factor)),
#                            x50 = fit.RSF_50$summary.random[-1] %>%
#                              bind_rows(.id = "param") %>%
#                              mutate(across(param, factor)))
#
# random.coeffs <- random.coeffs.list %>%
#   bind_rows(.id = "dataset")
# levels(random.coeffs$param) <- unique(fixed.coeffs$param)


# Add population means to random effects
# random.coeffs <- random.coeffs %>%
#   group_by(dataset, param) %>%
#   group_split() %>%
#   map2(.x = .,
#        .y = fixed.coeffs %>%
#          group_by(dataset, param) %>%
#          group_split(),
#        ~mutate(.x,
#                mean = mean + .y$mean,
#                `0.025quant` = `0.025quant` + .y$mean,
#                `0.975quant` = `0.975quant` + .y$mean)) %>%
#   bind_rows()

# ggplot(random.coeffs) +
#   geom_hline(yintercept = 0, linewidth = 0.75) +
#   geom_linerange(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`, color = param)) +
#   geom_point(aes(y = mean, x = ID, color = param)) +
#   geom_linerange(data = fixed.coeffs, aes(ymin = `0.025quant`, ymax = `0.975quant`), x = 50) +
#   geom_point(data = fixed.coeffs, aes(y = mean), x = 50) +
#   theme_bw() +
#   facet_grid(dataset ~ param, scales = "free")






####################################################
### Viz marginal effects plots per environ covar ###
####################################################


# random.coeffs2 <- random.coeffs %>%
#   filter(dataset == 'x10') %>%
#   dplyr::select(param, ID, mean, `0.025quant`, `0.975quant`)


### Bathymetry ###

bathym.newdata <- data.frame(bathym = seq(min(rsf.pts_10s$log.bathym), max(rsf.pts_10s$log.bathym), length.out = 500),
                             bathym.2 = seq(min(rsf.pts_10s$log.bathym), max(rsf.pts_10s$log.bathym), length.out = 500) ^ 2,
                             npp = 0,
                             npp.2 = 0,
                             sst = 0,
                             sst.2 = 0) %>%
  as.matrix()




# pred.bathym <- vector("list", length = n_distinct(random.coeffs2$ID))
# names(pred.bathym) <- unique(random.coeffs2$ID)
#
# for (i in 1:n_distinct(random.coeffs2$ID)) {
#
#   coeff1 <- random.coeffs2 %>%
#     filter(ID == unique(random.coeffs2$ID)[i]) %>%
#     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
#     as.matrix()
#
#   tmp <- bathym.newdata %*% coeff1 %>%
#     data.frame() %>%
#     mutate(bathym = exp(bathym.newdata[,1]))
#
#   pred.bathym[[i]] <- tmp
# }
# pred.bathym <- pred.bathym %>%
#   bind_rows(.id = "id")

coeffs <- fixed.coeffs %>%
  mutate(coeff = rep(rownames(fit.RSF_corr$summary.fixed[-1,]), 2)) %>%
  select(dataset, mean, coeff) %>%
  pivot_wider(names_from = dataset, values_from = mean, id_cols = coeff) %>%
  select(-coeff) %>%
  as.matrix()

pred.bathym.pop <- bathym.newdata %*% coeffs %>%
  data.frame() %>%
  mutate(bathym = exp(bathym.newdata[,1])) %>%
  pivot_longer(cols = -bathym, names_to = "method", values_to = "est")



# Pop mean estimate
ggplot() +
  geom_line(data = pred.bathym.pop, aes(x = bathym, y = exp(est), color = method), linewidth = 1) +
  theme_bw() +
  lims(x = c(0,300)) +
  labs(x = "Depth (m)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ method)

# ggplot() +
#   geom_line(data = pred.bathym, aes(x = bathym, y = exp(mean), group = id, color = id),
#             linewidth = 1, show.legend = FALSE) +
#   theme_bw() +
#   lims(x = c(0,300), y = c(0,1500)) +
#   labs(x = "Depth (m)", y = "Relative Intensity of Use") +
#   theme(axis.title = element_text(size = 30),
#         axis.text = element_text(size = 24),
#         panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank())







### SST ###

sst.newdata <- data.frame(bathym = 0,
                          bathym.2 = 0,
                          npp = 0,
                          npp.2 = 0,
                          sst = seq(min(rsf.pts_10s$log.sst), max(rsf.pts_10s$log.sst), length.out = 500),
                          sst.2 = seq(min(rsf.pts_10s$log.sst), max(rsf.pts_10s$log.sst), length.out = 500) ^ 2) %>%
  as.matrix()






# pred.sst <- vector("list", length = n_distinct(random.coeffs2$ID))
# names(pred.sst) <- unique(random.coeffs2$ID)
#
# for (i in 1:n_distinct(random.coeffs2$ID)) {
#
#   coeff1 <- random.coeffs2 %>%
#     filter(ID == unique(random.coeffs2$ID)[i]) %>%
#     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
#     as.matrix()
#
#   tmp <- sst.newdata %*% coeff1 %>%
#     data.frame() %>%
#     # mutate(sst = (sst.newdata[,7] * sd(rsf.pts_10s$sst, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$sst, na.rm = T))
#     mutate(sst = exp(sst.newdata[,5]))
#
#   pred.sst[[i]] <- tmp
# }
#
# pred.sst <- pred.sst %>%
#   bind_rows(.id = "id")



pred.sst.pop <- sst.newdata %*% coeffs %>%
  data.frame() %>%
  mutate(sst = exp(sst.newdata[,5])) %>%
  pivot_longer(cols = -sst, names_to = "method", values_to = "est")


# Pop mean estimate
ggplot() +
  geom_line(data = pred.sst.pop, aes(x = sst, y = exp(est), color = method), linewidth = 1.5) +
  theme_bw() +
  labs(x = "SST (°C)", y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ method, scales = "free_y")

# ggplot() +
#   geom_line(data = pred.sst, aes(x = sst, y = exp(mean), group = id, color = id),
#             linewidth = 1, show.legend = FALSE) +
#   theme_bw() +
#   lims(y = c(0,2e140)) +
#   labs(x = "SST (°C)", y = "Relative Intensity of Use") +
#   theme(axis.title = element_text(size = 30),
#         axis.text = element_text(size = 24),
#         panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank())





### NPP ###

npp.newdata <- data.frame(bathym = 0,
                          bathym.2 = 0,
                          npp = seq(min(rsf.pts_10s$log.npp), max(rsf.pts_10s$log.npp), length.out = 500),
                          npp.2 = seq(min(rsf.pts_10s$log.npp), max(rsf.pts_10s$log.npp), length.out = 500) ^ 2,
                          sst = 0,
                          sst.2 = 0) %>%
  as.matrix()





# pred.npp <- vector("list", length = n_distinct(random.coeffs2$ID))
# names(pred.npp) <- unique(random.coeffs2$ID)
#
# for (i in 1:n_distinct(random.coeffs2$ID)) {
#
#   coeff1 <- random.coeffs2 %>%
#     filter(ID == unique(random.coeffs2$ID)[i]) %>%
#     dplyr::select(mean, `0.025quant`, `0.975quant`) %>%
#     as.matrix()
#
#   tmp <- npp.newdata %*% coeff1 %>%
#     data.frame() %>%
#     # mutate(npp = (npp.newdata[,7] * sd(rsf.pts_10s$npp, na.rm = T)) + mean(rsf.pts_10s[obs.ind_10,]$npp, na.rm = T))
#     mutate(npp = exp(npp.newdata[,3]))
#
#   pred.npp[[i]] <- tmp
# }
#
# pred.npp <- pred.npp %>%
#   bind_rows(.id = "id")



pred.npp.pop <- npp.newdata %*% coeffs %>%
  data.frame() %>%
  mutate(npp = exp(npp.newdata[,3])) %>%
  pivot_longer(cols = -npp, names_to = "method", values_to = "est")

# Pop mean estimate
ggplot() +
  geom_line(data = pred.npp.pop, aes(x = npp / 1000, y = exp(est), color = method), linewidth = 1.5) +
  theme_bw() +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24)) +
  facet_wrap(~ method)

# ggplot() +
#   geom_line(data = pred.npp, aes(x = npp / 1000, y = exp(mean), group = id, color = id),
#             linewidth = 1, show.legend = FALSE) +
#   theme_bw() +
#   lims(y = c(0,1000)) +
#   labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
#   theme(axis.title = element_text(size = 30),
#         axis.text = element_text(size = 24),
#         panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank())






#####################################
### Load environmental covariates ###
#####################################

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "GoM", full.names = TRUE)
# files <- files[!grepl(pattern = "example", files)]  #remove any example datasets
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

# coeff1 <- fit.RSF_10$summary.fixed$mean[-1]
mean.pred <- as.matrix(newdat) %*% coeffs  #make predictions


rast.pred <- c(cov_list$bathym, cov_list$bathym)
names(rast.pred) <- colnames(mean.pred)
terra::values(rast.pred[[1]]) <- mean.pred[,1]
terra::values(rast.pred[[2]]) <- mean.pred[,2]
plot(rast.pred)


rast.pred.df <- as.data.frame(rast.pred, xy = TRUE) |>
  pivot_longer(cols = c(corr, hybrid), names_to = "method", values_to = "pred")
# names(rast.pred.df)[3] <- "pred"
bbox <- ext(rast.pred)

tmp.pts <- rsf.pts_10s %>%
  filter(month.year == "2020-09-01", obs == 1)

ggplot() +
  geom_raster(data = rast.pred.df, aes(x, y, fill = pred)) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_point(data = tmp.pts, aes(x, y), color = "forestgreen", alpha = 0.7, size = 2) +
  labs(x="",y="", title = "Population Mean: September 2020") +
  theme_bw() +
  coord_sf(xlim = c(bbox[1], bbox[2]),
           ylim = c(bbox[3], bbox[4]),
           expand = FALSE) +
  theme(axis.text = element_text(size = 12),
        plot.title = element_text(size = 26, face = "bold"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 20)) +
  facet_wrap(~ method)






# For ID 181796 in 2020-08
# rsf.pts_10s %>%
#   filter(id == 181796, obs == 1) %>%
#   group_by(month.year) %>%
#   count()  # most obs are in August
#
# newdat <- data.frame(log.bathym = terra::values(cov_list$bathym) %>%
#                        abs() %>%
#                        log(),
#                      log.bathym2 = terra::values(cov_list$bathym) %>%
#                        abs() %>%
#                        log() %>%
#                        apply(., 2, function(x) x^2),
#                      log.npp = terra::values(cov_list$npp$`2020-08-01`) %>%
#                        log(),
#                      log.npp2 = terra::values(cov_list$npp$`2020-08-01`) %>%
#                        log() %>%
#                        apply(., 2, function(x) x^2),
#                      log.sst = terra::values(cov_list$sst$`2020-08-01`) %>%
#                        log(),
#                      log.sst2 = terra::values(cov_list$sst$`2020-08-01`) %>%
#                        log() %>%
#                        apply(., 2, function(x) x^2))
# names(newdat) <- c('log.bathym','log.bathym2','log.npp','log.npp2','log.sst','log.sst2')
# summary(newdat)
#
# # Specify terms for ID 181796
# ind <- which(unique(rsf.pts_10s$id) == 181796)
# coeff.181796 <- random.coeffs2 %>%
#   filter(ID == ind) %>%
#   pull(mean)
# x181796.pred <- as.matrix(newdat) %*% coeff.181796  #make predictions
#
#
# rast.pred.181796 <- cov_list$bathym
# terra::values(rast.pred.181796) <- exp(x181796.pred)
#
#
# rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
# names(rast.pred.181796.df)[3] <- "pred"
#
# ggplot() +
#   geom_raster(data = rast.pred.181796.df, aes(x, y, fill = log(pred))) +
#   scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
#   geom_sf(data = gom.sf) +
#   geom_path(data = rsf.pts_10s %>%
#               filter(id == 181796, month.year == '2020-08-01', obs == 1), aes(x, y),
#             color = "chartreuse", linewidth = 0.5) +
#   labs(x="",y="", title = "ID 181796: August 2020") +
#   theme_bw() +
#   coord_sf(xlim = c(bbox[1], bbox[2]),
#            ylim = c(bbox[3], bbox[4]),
#            expand = FALSE) +
#   theme(axis.text = element_text(size = 12),
#         plot.title = element_text(size = 26, face = "bold"),
#         legend.title = element_text(size = 18),
#         legend.text = element_text(size = 16)) +
#   guides(fill = guide_colourbar(barwidth = 2, barheight = 20))






###########################
### Export model object ###
###########################

saveRDS(fit.RSF_corr, "Data_products/HGLM_corr_model_fit.rds")
saveRDS(fit.RSF_hybrid, "Data_products/HGLM_hybrid_model_fit.rds")


