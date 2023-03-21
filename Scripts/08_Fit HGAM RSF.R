
### Fit RSF as HGAM ###

library(tidyverse)
library(lubridate)
library(mgcv)
library(gratia)
library(bayesmove)
library(terra)
library(sf)
library(sfarrow)
library(tictoc)
library(patchwork)
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
  pivot_longer(cols = c(bathym, npp, sst), names_to = "covar", values_to = "value") %>%
  ggplot() +
  geom_density(aes(value, fill = factor(obs))) +
  theme_bw() +
  facet_wrap(~ covar, scales = "free", ncol = 1)

# Discretize covars to increase speed and efficiency of model fit
# bin.lims <- sapply(rsf.pts_10s[,c('bathym','npp','sst')], function(x) list(quantile(x, seq(0,1,length.out=25))))
# rsf.pts_10s <- rsf.pts_10s %>%
#   bayesmove::discrete_move_var(.,
#                                lims = bin.lims,
#                                varIn = c("bathym","npp","sst"),
#                                varOut = c("bathym.d","npp.d","sst.d")) %>%
#   mutate(across(bathym.d:sst.d, factor))
#
# bin.labs <- sapply(bin.lims, function(x) x[-1] - (diff(x) / 2)) %>%
#   data.frame()
# levels(rsf.pts_10s$bathym.d) <- bin.labs$bathym
# levels(rsf.pts_10s$npp.d) <- bin.labs$npp
# levels(rsf.pts_10s$sst.d) <- bin.labs$sst
#
# rsf.pts_10s <- rsf.pts_10s %>%
#   mutate(across(bathym.d:sst.d, function(x) as.numeric(as.character(x))))


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
cor(rsf.pts_10s[,c('bathym','npp','sst')])  #all corr low (< 0.177)


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

rsf.pts_30s2 <- rsf.pts_30s %>%
  mutate(across(id, factor))

rsf.pts_50s2 <- rsf.pts_50s %>%
  mutate(across(id, factor))



# Check concurvity from base model
# set.seed(2023)
# tic()
# fit.GAM10 <- bam(obs/wts2 ~ s(log.bathym, bs = "cr", k = 5, m = 2) +
#                        s(log.npp, bs = "cr", k = 5, m = 2) +
#                        s(log.sst, bs = "cr", k = 5, m = 2) +
#                        s(id, bs = "re"), data = rsf.pts_10s2, method = "fREML",
#                      family = poisson(), weights = wts2, discrete = TRUE)
# toc()  #took 1 sec to run

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
                  family = poisson(), weights = rsf.pts_10s2$wts2, discrete = 500)
toc()  #took 13 min to run

summary(fit.HGAM10_PI)
fit.HGAM10_PI$iter
plot(fit.HGAM10_PI, select = 1, scale = 0, shade = TRUE, shade.col = "lightblue")
plot(fit.HGAM10_PI, select = 51, scale = 0, shade = TRUE, shade.col = "lightblue")
plot(fit.HGAM10_PI, select = 101, scale = 0, shade = TRUE, shade.col = "lightblue")
gam.check(fit.HGAM10)



# Create custom partial effects plot
# evaluate the smooths
sm <- smooth_estimates(fit.HGAM10_PI, n = 500) %>%
  add_confint()
sm

# add partial residuals to data
# rsf.pts_10s3 <- rsf.pts_10s2 %>%
#   add_partial_residuals(fit.HGAM10_PI)


p.bathym.pop <- sm %>%
  filter(smooth == "s(log.bathym)") %>%
  ggplot() +
  # geom_rug(aes(x = exp(log.bathym)),
  #          data = rsf.pts_10s3,
  #          sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  # geom_ribbon(aes(ymin = exp(lower_ci), ymax = exp(upper_ci), x = exp(log.bathym)),
  #             alpha = 0.5, fill = "steelblue3") +
  # geom_point(aes(x = exp(log.bathym), y = `s(log.bathym)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "steelblue3", alpha = 0.2) +
  geom_line(aes(x = exp(log.bathym), y = exp(est)), linewidth = 1.2) +

  labs(x = 'Depth (m)', y = "Relative Intensity of Use") +
  lims(x = c(0,300)) +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))


p.bathym.id <- sm %>%
  drop_na(log.bathym, id) %>%
  ggplot() +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.bathym),
                y = exp(est + sm[sm$smooth == 's(log.bathym)',]$est), color = id),
            linewidth = 1) +
  labs(x = 'Depth (m)', y = "Relative Intensity of Use") +
  # lims(x = c(0,300), y = c(-600, 600)) +
  lims(y = c(0, 1e15)) +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        legend.position = "none",
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())




p.npp.pop <- sm %>%
  filter(smooth == "s(log.npp)") %>%
  ggplot() +
  # geom_rug(aes(x = exp(log.npp) / 1000),
  #          data = rsf.pts_10s3,
  #          sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  # geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = exp(log.npp) / 1000),
  #             alpha = 0.5, fill = 'darkgreen') +
  # geom_point(aes(x = exp(log.npp) / 1000, y = `s(log.npp)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "darkgreen", alpha = 0.2) +
  geom_line(aes(x = exp(log.npp) / 1000, y = exp(est)), linewidth = 1.2) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  # lims(x = c(0,30)) +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))


p.npp.id <- sm %>%
  drop_na(log.npp, id) %>%
  ggplot() +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.npp) / 1000,
                y = exp(est + sm[sm$smooth == 's(log.npp)',]$est), color = id),
            linewidth = 1) +
  labs(x = expression(paste("Net Primary Productivity (", g~C~m^-2~d^-1, ")")), y = "Relative Intensity of Use") +
  lims(y = c(0, 3000)) +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        legend.position = "none",
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())




p.sst.pop <- sm %>%
  filter(smooth == "s(log.sst)") %>%
  ggplot() +
  # geom_rug(aes(x = exp(log.sst)),
  #          data = rsf.pts_10s3,
  #          sides = "b", length = grid::unit(0.02, "npc")) +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  # geom_ribbon(aes(ymin = exp(lower_ci), ymax = exp(upper_ci), x = exp(log.sst)),
  #             alpha = 0.5, fill = 'firebrick') +
  # geom_point(aes(x = exp(log.sst), y = `s(log.sst)`),
  #            data = rsf.pts_10s3, cex = 1.5, colour = "darkgreen", alpha = 0.2) +
  geom_line(aes(x = exp(log.sst), y = exp(est)), linewidth = 1.2) +
  labs(x = 'SST (°C)', y = "Relative Intensity of Use") +
  theme_bw() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24))


p.sst.id <- sm %>%
  drop_na(log.sst, id) %>%
  ggplot() +
  # geom_hline(yintercept = 0, linewidth = 1, linetype = "dashed") +
  geom_line(aes(x = exp(log.sst),
                y = exp(est + sm[sm$smooth == 's(log.sst)',]$est), color = id),
            linewidth = 1) +
  labs(x = 'SST (°C)', y = "Relative Intensity of Use") +
  theme_bw() +
  lims(y = c(0,650)) +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        legend.position = "none",
        panel.background = element_rect(fill = "black"),
        panel.grid = element_blank())



p.bathym.pop + p.bathym.id
p.npp.pop + p.npp.id
p.sst.pop + p.sst.id





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
}


## Set all positive bathymetric values (i.e., elevation) as NA
cov_list[["bathym"]][cov_list[["bathym"]] > 0] <- NA


## Transform raster layers to match coarsest spatial resolution (i.e., NPP/Kd490)
for (var in c("bathym", "sst")) {
  cov_list[[var]] <- terra::resample(cov_list[[var]], cov_list$npp, method = "average")
}

## Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list[["bathym"]][cov_list[["bathym"]] > -1e-9] <- NA


## Transform CRS to match tracks
cov_list <- map(cov_list, terra::project, 'EPSG:3395')




##########################################
### Create example predictive surfaces ###
##########################################

# For pop mean in 2020-09
newdat <- data.frame(log.bathym = terra::values(cov_list$bathym) %>%
                        abs() %>%
                        log(),
                      log.npp = terra::values(cov_list$npp$`2020-09-01`) %>%
                        log(),
                      log.sst = terra::values(cov_list$sst$`2020-09-01`) %>%
                        log(),
                     id = unique(rsf.pts_10s2$id)[45])
names(newdat)[1:3] <- c('log.bathym','log.npp','log.sst')
summary(newdat)

# Specify terms to exclude
# var.eff <- grep(x = smooths(fit.HGAM10_PI), pattern = "id", value = T)


mean.pred <- predict.bam(fit.HGAM10_PI, newdata = newdat, type = "terms", terms = c("s(log.bathym)","s(log.npp)","s(log.sst)"),
                         discrete = FALSE, na.action = na.pass, se.fit = TRUE)
# mean.pred <- predict.bam(fit.HGAM10_PI, newdata = newdat, type = "link", exclude = var.eff,
#                          discrete = FALSE, na.action = na.pass)
summary(mean.pred$fit)

rast.pred <- cov_list$bathym
terra::values(rast.pred) <- exp(rowSums(mean.pred$fit))
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
rsf.pts_10s2 %>%
  filter(id == 181796, obs == 1) %>%
  group_by(month.year) %>%
  count()  # most obs are in August

newdat <- data.frame(log.bathym = terra::values(cov_list$bathym) %>%
                       abs() %>%
                       log(),
                     log.npp = terra::values(cov_list$npp$`2020-08-01`) %>%
                       log(),
                     log.sst = terra::values(cov_list$sst$`2020-08-01`) %>%
                       log(),
                     id = unique(rsf.pts_10s2$id)[45])
names(newdat)[1:3] <- c('log.bathym','log.npp','log.sst')
summary(newdat)

# Specify terms to include
var.eff <- grep(x = smooths(fit.HGAM10_PI), pattern = "181796", value = T)

x181796.pred <- predict.bam(fit.HGAM10_PI, newdata = newdat, type = "terms", terms = var.eff[1],
                         discrete = FALSE, na.action = na.pass, se.fit = TRUE)

summary(x181796.pred$fit)

rast.pred.181796 <- cov_list$bathym
terra::values(rast.pred.181796) <- exp(rowSums(x181796.pred$fit))


rast.pred.181796.df <- as.data.frame(rast.pred.181796, xy = TRUE)
names(rast.pred.181796.df)[3] <- "pred"

ggplot() +
  geom_raster(data = rast.pred.181796.df, aes(x, y, fill = log(pred))) +
  scale_fill_viridis_c("log(Intensity)", option = 'inferno') +
  geom_sf(data = gom.sf) +
  geom_path(data = rsf.pts_10s2 %>%
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








###########################
### Export model object ###
###########################

saveRDS(fit.HGAM10_PI, "Data_products/HGAM_model_fit.rds")



