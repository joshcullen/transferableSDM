
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
library(patchwork)
library(MetBrewer)
library(ggh4x)

source('Scripts/helper functions.R')


##########################
### Load fitted models ###
##########################

hgpr.fit <- readRDS("Data_products/HGPR_corr_model_fit.rds")
hgpr.age <- readRDS("Data_products/HGPR_corr_model_fit_scale_age.rds")




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







#########################################
### Validate correlative HGPR w/o age ###
#########################################

# Define vector of covar names
covars <- c("log.bathym","log.npp","log.sst")

# Define 1D meshes to be used for prediction across sites
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)




### Brazil ###

# Create vector storing all month-year values
my.ind.br <- names(cov_list_br$npp)

# Make predictions from model
tic()
br.rast.no_age <- predict.hgpr(cov_list = cov_list_br, model_fit = hgpr.fit, covars = covars,
                                   mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE,
                               method = "corr")
skrrrahh('khaled2')
toc()  #took 20 sec to run


# Normalize predictions on 0-1 scale
br.rast.no_age2 <- normalize(br.rast.no_age)


## Assess model performance via Continuous Boyce Index ##

boyce.br.full.no_age <- vector("list", nlyr(br.rast.no_age2))  #all locations from Brazil
boyce.br.sub.no_age <- vector("list", nlyr(br.rast.no_age2))  #subset from Brazil mainland

# Calculate by month.year
tic()
for (i in 1:nlyr(br.rast.no_age2)) {

  # Full set of Brazil locations
  obs_full <- dat.br %>%
    filter(month.year == my.ind.br[i]) %>%
    dplyr::select(x, y)

  # Subset of Brazil locations from mainland
  obs_sub <- dat.br %>%
    filter(month.year == my.ind.br[i], x < -3800000) %>%  #filter out locations east of this longitude
    dplyr::select(x, y)


  # Calculate Boyce Index per dataset
  boyce.br.full.no_age[[i]] <- boyce(fit = br.rast.no_age2[[i]],
                                   obs = obs_full,
                                   nbins = 10,
                                   bin.method = "seq",
                                   PEplot = FALSE,
                                   rm.duplicate = FALSE,
                                   method = "spearman")

  boyce.br.sub.no_age[[i]] <- boyce(fit = br.rast.no_age2[[i]],
                                  obs = obs_sub,
                                  nbins = 10,
                                  bin.method = "seq",
                                  PEplot = FALSE,
                                  rm.duplicate = FALSE,
                                  method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 2 sec


# Summarize results per dataset
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
             Region = "Brazil_main",
             Method = "No_Age")





### Qatar ###

# Create vector storing all month-year values
my.ind.qa <- names(cov_list_qa$npp)

# Make predictions from model
tic()
qa.rast.no_age <- predict.hgpr(cov_list = cov_list_qa, model_fit = hgpr.fit, covars = covars,
                               mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = FALSE,
                               method = "corr")
skrrrahh('khaled2')
toc()  #took 1 sec to run


# Normalize predictions on 0-1 scale
qa.rast.no_age2 <- normalize(qa.rast.no_age)


## Assess model performance via Continuous Boyce Index ##
boyce.qa.no_age <- vector("list", nlyr(qa.rast.no_age2))

# Make spatial predictions per month.year
tic()
for (i in 1:nlyr(qa.rast.no_age2)) {

  # Subset tracks by month.year
  obs <- dat.qa %>%
    filter(month.year == my.ind.qa[i]) %>%
    dplyr::select(x, y)

  boyce.qa.no_age[[i]] <- boyce(fit = qa.rast.no_age2[[i]],
                                     obs = obs,
                                     nbins = 10,
                                     bin.method = "seq",
                                     PEplot = FALSE,
                                     rm.duplicate = FALSE,
                                     method = "spearman")
}
skrrrahh("khaled3")
toc()  #took 1 sec


# Summarize results
boyce.qa.no_age2 <- boyce.qa.no_age %>%
  map(., pluck, "cor") %>%
  unlist() %>%
  data.frame(Age.Class = "Pop",
             cor = .,
             Region = "Qatar",
             Method = "No_Age")







########################################
### Validate correlative HGPR w/ age ###
########################################

# Define vector of covar names
covars <- c("log.bathym","log.npp","log.sst")

# Define 1D meshes to be used for prediction across sites
mesh.seq <- list(log.bathym = c(0.001, 5500),
                 log.npp = c(20, 200000),
                 log.sst = c(12,35)) %>%
  map(log)




### Brazil ###

# Make predictions from model
tic()
br.rast.age <- predict.hgpr(cov_list = cov_list_br, model_fit = hgpr.age, covars = covars,
                               mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = TRUE,
                            method = "corr")
skrrrahh('khaled2')
toc()  #took 30 sec to run


# Normalize predictions on 0-1 scale
br.rast.age2 <- map(br.rast.age, normalize)


## Assess model performance via Continuous Boyce Index ##
boyce.br.full.age <- boyce.br.sub.age <- vector("list", length(br.rast.age2)) %>%
  set_names(names(br.rast.age2))

# Create vectors storing all month-year values and life stages
my.ind.br <- names(cov_list_br$npp)
age.class <- c("Juv","Adult")

# Calculate by month.year
tic()
for (j in 1:length(br.rast.age2)) {
  for (i in 1:nlyr(br.rast.age2[[j]])) {

    # Full set of Brazil locations
    obs_full <- dat.br %>%
      filter(month.year == my.ind.br[i], Age == age.class[j]) %>%
      dplyr::select(x, y)

    # Subset of Brazil locations from mainland
    obs_sub <- dat.br %>%
      filter(month.year == my.ind.br[i], x < -3800000, Age == age.class[j]) %>%  #filter out locations east of this longitude
      dplyr::select(x, y)


    # Calculate Boyce Index per dataset
    boyce.br.full.age[[j]][[i]] <- boyce(fit = br.rast.age2[[j]][[i]],
                                     obs = obs_full,
                                     nbins = 10,
                                     bin.method = "seq",
                                     PEplot = FALSE,
                                     rm.duplicate = FALSE,
                                     method = "spearman")

    boyce.br.sub.age[[j]][[i]] <- boyce(fit = br.rast.age2[[j]][[i]],
                                    obs = obs_sub,
                                    nbins = 10,
                                    bin.method = "seq",
                                    PEplot = FALSE,
                                    rm.duplicate = FALSE,
                                    method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 5 sec


# Summarize results per dataset
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
                 Region = "Brazil_main",
                 Method = "Age")}
  ) %>%
  bind_rows(.id = "Age.Class")





### Qatar ###

# Make predictions from model
tic()
qa.rast.age <- predict.hgpr(cov_list = cov_list_qa, model_fit = hgpr.age, covars = covars,
                            mesh.seq = mesh.seq, nbasis = 5, degree = 2, age.class = TRUE,
                            method = "corr")
skrrrahh('khaled2')
toc()  #took 1 sec to run


# Normalize predictions on 0-1 scale
qa.rast.age2 <- map(qa.rast.age, normalize)


## Assess model performance via Continuous Boyce Index ##
boyce.qa.age <- vector("list", length(qa.rast.age2)) %>%
  set_names(names(qa.rast.age2))

# Create vector storing all month-year values
my.ind.qa <- names(cov_list_qa$npp)

# Calculate by month.year
tic()
for (j in 1:length(qa.rast.age2)) {
  for (i in 1:nlyr(qa.rast.age2[[j]])) {

    # Subset tracks by month.year
    obs <- dat.qa %>%
      filter(month.year == my.ind.qa[i], Age == age.class[j]) %>%
      dplyr::select(x, y)


    # Calculate Boyce Index
    boyce.qa.age[[j]][[i]] <- boyce(fit = qa.rast.age2[[j]][[i]],
                                obs = obs,
                                nbins = 10,
                                bin.method = "seq",
                                PEplot = FALSE,
                                rm.duplicate = FALSE,
                                method = "spearman")
  }
}
skrrrahh("khaled3")
toc()  #took 1 sec


# Summarize results
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

# Merge all results together
boyce.fit <- rbind(boyce.br.full.no_age2, boyce.br.sub.no_age2, boyce.qa.no_age2,
                   boyce.br.full.age2, boyce.br.sub.age2, boyce.qa.age2) %>%
  mutate(across(Age.Class, \(x) factor(x, levels = c(age.class, "Pop"))))

# Calculate avg Boyce Index for Method x Region
boyce.mean <- boyce.fit %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()


ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(Method, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.mean, aes(x = Region, y = mean, group = Method),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = met.brewer("Isfahan1")[c(3,6)], guide = "none") +
  scale_fill_manual(values = met.brewer("Isfahan1")[c(3,6)], labels = c("Age", "No Age")) +
  scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))





# Create plot of Boyce Index broken out by age class (or pop.) for Supplement
boyce.mean2 <- boyce.fit %>%
  group_by(Method, Region, Age.Class) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()


ggplot(data = boyce.fit, aes(Region, cor)) +
  geom_point(aes(fill = Age.Class), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(Age.Class, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.mean2, aes(x = Region, y = mean, group = Age.Class),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  scale_fill_brewer("Group", palette = "Dark2", labels = c("Juvenile","Adult","Population")) +
  scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))

# ggsave("Tables_Figs/Figure S8.png", width = 7, height = 5, units = "in", dpi = 400)





### Create example prediction maps per method ###

# Define spatial extent for Brazil and Qatar
bbox.br <- ext(br.rast.no_age$`2016-05-01`)
bbox.qa <- ext(qa.rast.no_age$`2014-02-01`)

# Break rasters into bins used for Boyce Index
br.rast.no_age3 <- classify(br.rast.no_age2, seq(0, 1, by = 0.1))
br.rast.no_age3 <- br.rast.no_age3 + 1
br.rast.no_age3.df <- as.data.frame(br.rast.no_age3, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1)))

br.rast.age3 <- map(br.rast.age2, classify, seq(0, 1, by = 0.1))
br.rast.age3 <- map(br.rast.age3, function(x) x + 1)
br.rast.age3.df <- map(br.rast.age3, ~{.x %>%
    as.data.frame(., xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1)))
}) %>%
  bind_rows(.id = "Age")


qa.rast.no_age3 <- classify(qa.rast.no_age2, seq(0, 1, by = 0.1))
qa.rast.no_age3 <- qa.rast.no_age3 + 1
qa.rast.no_age3.df <- as.data.frame(qa.rast.no_age3, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1)))

qa.rast.age3 <- classify(qa.rast.age2$Juv, seq(0, 1, by = 0.1))
qa.rast.age3 <- qa.rast.age3 + 1
qa.rast.age3.df <- as.data.frame(qa.rast.age3, xy = TRUE) %>%
  mutate(across(3:ncol(.), \(x) factor(x, levels = 10:1)))




p.pop.br <- ggplot() +
  geom_raster(data = br.rast.no_age3.df, aes(x, y, fill = `2022-02-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = br.sf) +
  geom_text(aes(x = -4786738, y = -887521.1, label = "Brazil"), size = 5, fontface = "italic") +
  labs(x="",y="", title = "Population") +
  theme_bw() +
  coord_sf(xlim = c(bbox.br[1], bbox.br[2]),
           ylim = c(bbox.br[3], bbox.br[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = "bold"))

p.juv.br <- ggplot() +
  geom_raster(data = br.rast.age3.df %>%
                filter(Age == "Juv"), aes(x, y, fill = `2022-02-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = br.sf) +
  labs(x="",y="", title = "Juvenile") +
  theme_bw() +
  coord_sf(xlim = c(bbox.br[1], bbox.br[2]),
           ylim = c(bbox.br[3], bbox.br[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = "bold"))

p.adult.br <- ggplot() +
  geom_raster(data = br.rast.age3.df %>%
                filter(Age == "Adult"), aes(x, y, fill = `2022-02-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = br.sf) +
  labs(x="",y="", title = "Adult") +
  theme_bw() +
  coord_sf(xlim = c(bbox.br[1], bbox.br[2]),
           ylim = c(bbox.br[3], bbox.br[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = "bold"))



p.pop.qa <- ggplot() +
  geom_raster(data = qa.rast.no_age3.df, aes(x, y, fill = `2014-03-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = qa.sf) +
  geom_text(aes(x = 5700000, y = 2870000, label = "Qatar"), size = 5, fontface = "italic") +
  labs(x="",y="", title = "Population") +
  theme_bw() +
  coord_sf(xlim = c(bbox.qa[1], bbox.qa[2]),
           ylim = c(bbox.qa[3], bbox.qa[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = "bold"))


p.juv.qa <- ggplot() +
  geom_raster(data = qa.rast.age3.df, aes(x, y, fill = `2014-03-01`)) +
  scale_fill_viridis_d("HS Bins", option = 'inferno', direction = -1, drop = FALSE) +
  geom_sf(data = qa.sf) +
  labs(x="",y="", title = "Juvenile") +
  theme_bw() +
  coord_sf(xlim = c(bbox.qa[1], bbox.qa[2]),
           ylim = c(bbox.qa[3], bbox.qa[4]),
           expand = FALSE) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = "bold"))



# Make composite plot
p.pop.br + p.juv.br + p.adult.br + p.pop.qa + p.juv.qa + guide_area() +
  plot_layout(nrow = 2, guides = "collect")

# ggsave("Tables_Figs/Figure 6.png", width = 5, height = 5.5, units = "in", dpi = 400)






###############################################################
### Create composite plot of marginal effects across models ###
###############################################################


# Define 1D mesh per covar
mesh.list <- vector("list", length(covars))
for (i in 1:length(covars)) {
  mesh.list[[i]] <- inla.mesh.1d(seq(mesh.seq[[i]][1], mesh.seq[[i]][2],
                                     length.out = 5),
                                 degree = 2,
                                 boundary = 'free')
}



# Define predictive seq of covars
gp_vars <- data.frame(log.bathym = seq(mesh.seq$log.bathym[1], mesh.seq$log.bathym[2], length.out = 500),
                   log.npp = seq(mesh.seq$log.npp[1], mesh.seq$log.npp[2], length.out = 500),
                   log.sst = seq(mesh.seq$log.sst[1], mesh.seq$log.sst[2], length.out = 500))


# Generate A matrices for prediction
A.mat <- vector("list", length(covars))
for (j in 1:length(covars)) { #one matrix for model estimation and another for generating predictions for plotting
  A.mat[[j]] <- inla.spde.make.A(mesh.list[[j]], loc = gp_vars[[covars[[j]]]])
}
names(A.mat) <- covars

# Replicate list elements for mapping over lists
A.mat2 <- rep(A.mat, each = 2)






### Predict across models ###

# Define coeff values from HGPR
coeff1 <- hgpr.age$summary.random[1:3] %>%
  map(., ~dplyr::select(.x, mean)) %>%
  map(., ~mutate(.x, age = rep(1:2, each = 6))) %>%
  map(., ~split(.x, .x$age)) %>%
  flatten() %>%
  map(., pull, mean) %>%
  set_names(paste(rep(covars, each = 2), rep(1:2, length(covars)), sep = "_"))


# Make predictions via linear algebra
marg.eff <- A.mat2 %>%
  map2(.x = ., .y = coeff1,
       ~{.x %*% .y %>%
           as.vector()}
  ) %>%
  set_names(names(coeff1)) %>%
  bind_cols(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
  split.default(., str_extract(names(.), "[0-9]$")) %>%  #split preds into list by age.class
  set_names(c("Juv","Adult")) %>%
  map(~{.x %>%
      pivot_longer(cols = everything(), names_to = "covar", values_to = "mean") %>%
      arrange(covar) %>%
      mutate(x = c(gp_vars$log.bathym, gp_vars$log.npp, gp_vars$log.sst)) %>%
      mutate(across(c(mean, x), \(z) exp(z))) %>%
      mutate(covar = case_when(str_detect(covar, "log.bathym_.") ~ "depth",
                               str_detect(covar, "log.npp_.") ~ "npp",
                               str_detect(covar, "log.sst_.") ~ "sst")) %>%
      filter(!(covar == "depth" & x > 300))  #to constrain plot to only show depth up to 300 m
    })




# Merge predictions for both life stages
marg.eff2 <- bind_rows(marg.eff, .id = "Age Class") %>%
  mutate(x = case_when(covar == "npp" ~ x / 1000,
                       TRUE ~ x),
         covar = case_when(covar == "depth" ~ "Depth (m)",
                           covar == "npp" ~ "NPP (g C m^-2 d^-1)",
                           covar == "sst" ~ "SST (°C)")) %>%
  mutate(`Age Class` = factor(`Age Class`, levels = c("Juv","Adult")))





# Plot pop-level marginal effects
ggplot() +
  geom_line(data = marg.eff2, aes(x = x, y = mean, color = `Age Class`), linewidth = 1) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  theme_bw(base_size = 14) +
  labs(x = "", y = "Relative Intensity of Use") +
  theme(strip.background = element_rect(fill = NA, color = NA),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank()
  ) +
  ggh4x::facet_grid2(`Age Class` ~ covar,
                     independent = "y",
                     scales = "free",
                     switch = "x")

# ggsave("Tables_Figs/Figure S7.png", width = 11, height = 9, units = "in", dpi = 400)



##################################
### Export Boyce Index results ###
##################################

write_csv(boyce.fit, "Data_products/boyce_age_results.csv")

