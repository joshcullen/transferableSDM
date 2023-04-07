
### Compare model transferability by accounting for age class (or not) ###

library(tidyverse)
library(INLA)
library(terra)
library(sfarrow)
library(tictoc)
library(lubridate)
library(BRRR)
library(tidyterra)
library(patchwork)
library(future)
library(furrr)

source('Scripts/helper functions.R')


##########################
### Load fitted models ###
##########################

hgpr.fit <- readRDS("Data_products/HGPR_model_fit.rds")
hgpr.age <- readRDS("Data_products/HGPR_model_fit_scale_age.rds")




################################
### Load validation datasets ###
################################

dat.br <- read_csv("Processed_data/Brazil_Cm_Tracks_behav.csv")
dat.qa <- read_csv("Processed_data/Qatar_Cm_Tracks_behav.csv")
dat.age <- read_csv("Processed_data/Prefiltered_Cm_Tracks.csv") %>%
  filter(Region != "GoM") %>%
  dplyr::select(Age, Ptt) %>%
  distinct()

# Create indexing column "month.year" and only retain Resident locs
dat.br <- dat.br %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  filter(behav == 'Resident') %>%
  left_join(dat.age, by = c("id" = "Ptt"))

dat.qa <- dat.qa %>%
  mutate(month.year = as_date(date), .after = 'date') %>%
  mutate(month.year = str_replace(month.year, pattern = "..$", replacement = "01")) %>%
  mutate(across(id, as.character)) %>%
  filter(behav == 'Resident') %>%
  left_join(dat.age, by = c("id" = "Ptt"))


# Load spatial land layers
br.sf <- st_read_parquet("Environ_data/Brazil_land.parquet")
qa.sf <- st_read_parquet("Environ_data/Qatar_land.parquet")




########################################
### Load environmental raster layers ###
########################################

### Brazil ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Brazil", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_br <- sapply(files, rast)
names(cov_list_br) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_br[[var]]) <- gsub(names(cov_list_br[[var]]), pattern = "-..$", replacement = "-01")
}

# Set all positive bathymetric values (i.e., elevation) as NA
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > 0] <- NA

# Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list_br[[var]] <- terra::resample(cov_list_br[[var]], cov_list_br$npp, method = "average")
}

# Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_br[["bathym"]][cov_list_br[["bathym"]] > -1e-9] <- NA


# Transform CRS to match tracks
cov_list_br <- map(cov_list_br, terra::project, 'EPSG:3395')






### Qatar ###

## Load in environ rasters
files <- list.files(path = 'Environ_data', pattern = "Qatar", full.names = TRUE)
files <- files[grepl(pattern = "tif", files)]  #only keep GeoTIFFs

# Merge into list; each element is a different covariate
cov_list_qa <- sapply(files, rast)
names(cov_list_qa) <- c('bathym', 'npp', 'sst')

# Change names for dynamic layers to match YYYY-MM-01 format
for (var in c('npp', 'sst')) {
  names(cov_list_qa[[var]]) <- gsub(names(cov_list_qa[[var]]), pattern = "-..$", replacement = "-01")
}

# Set all positive bathymetric values (i.e., elevation) as NA
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > 0] <- NA

# Transform raster layers to match coarsest spatial resolution (i.e., NPP)
for (var in c("bathym", "sst")) {
  cov_list_qa[[var]] <- terra::resample(cov_list_qa[[var]], cov_list_qa$npp, method = "average")
}

# Deal w/ bathym depth exactly equal to 0 (since a problem on log scale)
cov_list_qa[["bathym"]][cov_list_qa[["bathym"]] > -1e-9] <- NA


# Transform CRS to match tracks
cov_list_qa <- map(cov_list_qa, terra::project, 'EPSG:3395')







#############################
### Validate HGPR w/o age ###
#############################

# Define vector of covar names
covars <- c("log.bathym","log.npp","log.sst")

# Define 1D meshes to be used for prediction across sites
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)




### Brazil ###
my.ind.br <- names(cov_list_br$npp)

tic()
br.rast.no_age <- predict.hgpr(cov_list = cov_list_br, model_fit = hgpr.fit, covars = covars,
                                   mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE)
skrrrahh('khaled2')
toc()  #took 40 sec to run


# Assess model performance via Continuous Boyce Index
boyce.br.full.no_age <- vector("list", nlyr(br.rast.no_age))
boyce.br.sub.no_age <- vector("list", nlyr(br.rast.no_age))
tic()
for (i in 1:nlyr(br.rast.no_age)) {

  # Subset tracks by month.year
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3650000) %>%
    dplyr::select(x, y)

  boyce.br.full.no_age[[i]] <- boyce(fit = br.rast.no_age[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = TRUE,
                                   method = "spearman")

  boyce.br.sub.no_age[[i]] <- boyce(fit = br.rast.no_age[[i]],
                                  obs = obs_sub,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = TRUE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 5 sec



perc.use.br.full.no_age <- boyce.br.full.no_age %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.full.no_age, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #1 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.no_age %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



perc.use.br.sub.no_age <- boyce.br.sub.no_age %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.br.sub.no_age, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #1.5 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.no_age %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  ylim(0,1) +
  theme_bw()

boyce.br.full.no_age2 <- boyce.br.full.no_age %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(Age.Class = "Pop",
             cor = .,
             Region = "Brazil_all",
             Method = "No_Age")

boyce.br.sub.no_age2 <- boyce.br.sub.no_age %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(Age.Class = "Pop",
             cor = .,
             Region = "Brazil_sub",
             Method = "No_Age")





### Qatar ###
my.ind.qa <- names(cov_list_qa$npp)

tic()
qa.rast.no_age <- predict.hgpr(cov_list = cov_list_qa, model_fit = hgpr.fit, covars = covars,
                               mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE)
skrrrahh('khaled2')
toc()  #took 2 sec to run


# Assess model performance via Continuous Boyce Index
boyce.qa.no_age <- vector("list", nlyr(qa.rast.no_age))
tic()
for (i in 1:nlyr(qa.rast.no_age)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.no_age[[i]] <- boyce(fit = qa.rast.no_age[[i]],
                                     obs = obs,
                                     nbins = 10,
                                     bin.method = "seq",
                                     PEplot = FALSE,
                                     rm.duplicate = TRUE,
                                     method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec



perc.use.qa.no_age <- boyce.qa.no_age %>%
  map(., pluck, "perc.use") %>%
  set_names(1:length(.)) %>%
  bind_rows() %>%
  janitor::remove_empty(which = "cols") %>%
  apply(., 2, function(x) cumsum(rev(x)))

# check fewest bins that contain >=90% of all obs
apply(perc.use.qa.no_age, 2, function(x) which(x >= 0.9)[1]) %>%
  mean()  #2.2 bins

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.no_age %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1)) %>%
  pivot_longer(cols = -bin, names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw()



boyce.qa.no_age2 <- boyce.qa.no_age %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(Age.Class = "Pop",
             cor = .,
             Region = "Qatar",
             Method = "No_Age")







############################
### Validate HGPR w/ age ###
############################

# Define vector of covar names
covars <- c("log.bathym","log.npp","log.sst")

# Define 1D meshes to be used for prediction across sites
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,38)) %>%
  map(log)




### Brazil ###

tic()
br.rast.age <- predict.hgpr(cov_list = cov_list_br, model_fit = hgpr.age, covars = covars,
                               mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = TRUE)
skrrrahh('khaled2')
toc()  #took 1 min to run


# Assess model performance via Continuous Boyce Index
boyce.br.full.age <- boyce.br.sub.age <- vector("list", length(br.rast.age)) %>%
  set_names(names(br.rast.age))

my.ind.br <- names(cov_list_br$npp)
age.class <- c("Juv","Adult")

tic()
for (j in 1:length(br.rast.age)) {
  for (i in 1:nlyr(br.rast.age[[j]])) {

    # Subset tracks by month.year
    obs_full <- dat.br %>%
      filter(month.year == my.ind.br[i], Age == age.class[j]) %>%
      dplyr::select(x, y)

    obs_sub <- dat.br %>%
      filter(month.year == my.ind.br[i], x < -3650000, Age == age.class[j]) %>%
      dplyr::select(x, y)

    boyce.br.full.age[[j]][[i]] <- boyce(fit = br.rast.age[[j]][[i]],
                                     obs = obs_full,
                                     nbins = 10,
                                     bin.method = "seq",
                                     PEplot = FALSE,
                                     rm.duplicate = TRUE,
                                     method = "spearman")

    boyce.br.sub.age[[j]][[i]] <- boyce(fit = br.rast.age[[j]][[i]],
                                    obs = obs_sub,
                                    nbins = 10,
                                    bin.method = "seq",
                                    PEplot = FALSE,
                                    rm.duplicate = TRUE,
                                    method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 10 sec



perc.use.br.full.age <- boyce.br.full.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
  ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      # t() %>%
      data.frame()}
  )

# check fewest bins that contain >=90% of all obs
map(perc.use.br.full.age, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
    mean()}
)  #2.5 for Juv; 3.5 for Adult

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.full.age %>%
  map(~mutate(.x, bin = factor(10:1, levels = 10:1))) %>%
  bind_rows(.id = "Age") %>%
  mutate(across(Age, factor, levels = age.class)) %>%
  pivot_longer(cols = -c(Age, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Brazil_all") +
  facet_wrap(~Age)



perc.use.br.sub.age <- boyce.br.sub.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
  ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      # t() %>%
      data.frame()}
  )

# check fewest bins that contain >=90% of all obs
map(perc.use.br.sub.age, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
    mean()}
)  #2.5 for Juv; 2.5 for Adult

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.br.sub.age %>%
  map(~mutate(.x, bin = factor(10:1, levels = 10:1))) %>%
  bind_rows(.id = "Age") %>%
  mutate(across(Age, factor, levels = age.class)) %>%
  pivot_longer(cols = -c(Age, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Brazil_sub") +
  facet_wrap(~Age)



boyce.br.full.age2 <- boyce.br.full.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
  ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Brazil_all",
                 Method = "Age")}
  ) %>%
  bind_rows(.id = "Age.Class")

boyce.br.sub.age2 <- boyce.br.sub.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
  ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Brazil_sub",
                 Method = "Age")}
  ) %>%
  bind_rows(.id = "Age.Class")





### Qatar ###

tic()
qa.rast.age <- predict.hgpr(cov_list = cov_list_qa, model_fit = hgpr.age, covars = covars,
                            mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = TRUE)
skrrrahh('khaled2')
toc()  #took 2 sec to run



# Assess model performance via Continuous Boyce Index
boyce.qa.age <- vector("list", length(qa.rast.age)) %>%
  set_names(names(qa.rast.age))

my.ind.qa <- names(cov_list_qa$npp)

tic()
for (j in 1:length(qa.rast.age)) {
  for (i in 1:nlyr(qa.rast.age[[j]])) {

    # Subset tracks by month.year
    obs <- dat.qa %>%
      filter(month.year == my.ind.qa[i], Age == age.class[j]) %>%
      dplyr::select(x, y)

    boyce.qa.age[[j]][[i]] <- boyce(fit = qa.rast.age[[j]][[i]],
                                obs = obs,
                                nbins = 10,
                                bin.method = "seq",
                                PEplot = FALSE,
                                rm.duplicate = TRUE,
                                method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 1 sec



perc.use.qa.age <- boyce.qa.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("perc.use") %>%
      set_names(1:length(.))}
  ) %>%
  map_depth(., 1, ~{.x %>%
      bind_rows() %>%
      janitor::remove_empty(which = "rows") %>%
      dplyr::select(10:1) %>%
      apply(., 1, function(x) cumsum(x)) %>%
      data.frame()}
  )

# check fewest bins that contain >=90% of all obs
map(perc.use.qa.age, ~{apply(.x, 2, function(x) which(x >= 0.9)[1]) %>%
    mean()}
)  #4.75 bins; no adults, so remove

perc.use.qa.age <- perc.use.qa.age$Juv

# Viz plot of cumulative percentage of obs per bin (highest to lowest)
perc.use.qa.age %>%
  data.frame() %>%
  mutate(bin = factor(10:1, levels = 10:1),
         Age = factor("Juv", levels = age.class)) %>%
  pivot_longer(cols = -c(Age, bin), names_to = 'month.year', values_to = "cum.perc") %>%
  ggplot(aes(bin, cum.perc)) +
  geom_hline(yintercept = 0.9, linewidth = 0.75, linetype = "dashed", color = "red") +
  geom_line(aes(group = month.year, color = month.year)) +
  theme_bw() +
  labs(x = "Bin", y = "Cumulative Percentage of Total Observations", title = "Qatar") +
  facet_wrap(~Age)


boyce.qa.age2 <- boyce.qa.age %>%
  map_depth(., 2, ~{.x %>%
      pluck("cor")}
  ) %>%
  map(., ~{.x %>%
      unlist() %>%
      data.frame(cor = .,
                 Region = "Qatar",
                 Method = "Age")}
  ) %>%
  bind_rows(.id = "Age.Class")








####################################
### Summarize Validation Results ###
####################################

boyce.fit <- rbind(boyce.br.full.no_age2, boyce.br.sub.no_age2, boyce.qa.no_age2,
                   boyce.br.full.age2, boyce.br.sub.age2, boyce.qa.age2) %>%
  mutate(across(Age.Class, factor, levels = c(age.class, "Pop")))

boyce.mean <- boyce.fit %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 5, position = position_dodge(width = 0.75)) +
  geom_violin(aes(color = Method), fill = "transparent", position = position_dodge(width = 0.75)) +
  geom_point(data = boyce.mean, aes(x = Region, y = mean, group = Method),
             size = 6, position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  facet_wrap(~Age.Class)
