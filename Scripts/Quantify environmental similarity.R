
### Evaluate environmental similarity between training and testing regions ###
### Using Shape method from {flexsdm} ###

library(tidyverse)
library(terra)
library(flexsdm)
library(tictoc)
library(BRRR)
library(tidyterra)
library(patchwork)


#################
### Load data ###
#################

rsf.pts_10 <- read_csv("Processed_data/GoM_Cm_RSFprep_10x.csv") |>
  drop_na(bathym, npp, sst)



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






############################################
### Evaluate environmental dissimilarity ###
############################################


### Qatar ###
qa_sim <- vector("list", length = nlyr(cov_list_qa$sst)) |>
  set_names(names(cov_list_qa$sst))

tic()
for (i in seq_along(qa_sim)) {
  qa_sim[[i]] <- extra_eval(training_data = rsf.pts_10,
                            projection_data = rast(c(bathym = cov_list_qa[[1]],
                                                     npp = cov_list_qa[[2]][[i]],
                                                     sst = cov_list_qa[[3]][[i]])),
                            pr_ab = "obs",
                            n_cores = 20,
                            aggreg_factor = 1,
                            metric = "mahalanobis")

  cat("\nLayer", i, "of", length(qa_sim))
}
skrrrahh('khaled2')
toc()  #took 3.5 min w/ 20 cores

plot(qa_sim[[1]], main = "Qatar extrapolation pattern")



### Brazil ###
br_sim <- vector("list", length = nlyr(cov_list_br$sst)) |>
  set_names(names(cov_list_br$sst))

tic()
for (i in seq_along(br_sim)) {
  br_sim[[i]] <- extra_eval(training_data = rsf.pts_10,
                            projection_data = rast(c(bathym = cov_list_br[[1]],
                                                     npp = cov_list_br[[2]][[i]],
                                                     sst = cov_list_br[[3]][[i]])),
                            pr_ab = "obs",
                            n_cores = 20,
                            aggreg_factor = 1,
                            metric = "mahalanobis")

  cat("\nLayer", i, "of", length(br_sim))
}
skrrrahh('khaled2')
toc()  #took 5.75 hrs w/ 20 cores

plot(br_sim[[50]], main = "Brazil extrapolation pattern")
plot(rast(br_sim[1:12]), main = "Brazil extrapolation pattern")



##########################################
### Summarize similarity across months ###
##########################################

# Define index to subset cells that fall w/in neritic habitat (for Brazil only)
br_neritic_ind <- which(values(cov_list_br$bathym) > -200)


# Brazil
br_summ <- br_sim |>
  map(~values(.x) |>
        data.frame()) |>
  map(~summarize(.x,
                 # mean = mean(extrapolation, na.rm = TRUE),
                 # sd = sd(extrapolation, na.rm = TRUE),
                 mean = mean(extrapolation[br_neritic_ind], na.rm = TRUE),
                 sd = sd(extrapolation[br_neritic_ind], na.rm = TRUE))) |>
  bind_rows(.id = "month-year") |>
  mutate(Region = "Brazil")


# Qatar
qa_summ <- qa_sim |>
  map(~values(.x) |>
        data.frame()) |>
  map(~summarize(.x,
                 mean = mean(extrapolation, na.rm = TRUE),
                 sd = sd(extrapolation, na.rm = TRUE))) |>
  bind_rows(.id = "month-year") |>
  mutate(Region = "Qatar")



####################
### Plot results ###
####################

# Join both summaries and create labels for plot
valid_summ <- rbind(br_summ, qa_summ) |>
  mutate(month = month(as_date(`month-year`)))
month.year <- sort(unique(valid_summ$`month-year`))
y_labs <- ifelse(seq_along(month.year) %% 3 == 0, month.year, "")

# DF to highlight points that are mapped
valid_highlight <- valid_summ |>
  filter(Region == "Brazil" & `month-year` == "2020-02-01" |
         Region == "Qatar" & `month-year`== "2014-09-01")


# Plot values of similarity by month-year
p.summ <- ggplot(data = valid_summ) +
  geom_point(aes(mean, `month-year`, fill = month), shape = 21, size = 2) +
  geom_point(data = valid_highlight, aes(mean, `month-year`), shape = 21,
             fill = "transparent", size = 4, color = "blue", stroke = 2) +
  scale_fill_gradientn("Month",
                        colors = c(viridis::magma(n=6), rev(viridis::magma(n=6))),
                        breaks = seq(1, 12, by = 2),
                        labels = month.abb[seq(1, 12, by = 2)]) +
  scale_y_discrete(labels = y_labs) +
  xlim(-2, 32) +
  labs(x = "Mean Environmental Similarity",
       y = "Month-Year") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "grey")) +
  facet_wrap(~ Region)


# Map environ sim for Brazil
p.br <- ggplot() +
  geom_spatraster(data = br_sim[["2020-02-01"]]) +
  scale_fill_viridis_c("Environmental\nSimilarity\n",
                       option = "rocket", limits = c(0,50),
                       oob = scales::squish, labels = c(seq(0, 40, by = 10), "> 50")) +
  coord_sf(expand = FALSE) +
  labs(title = "Brazil: Feb 2020") +
  theme_void() +
  theme(plot.title = element_text(size = 16, face = "bold"))


# Map environ sim for Qatar
p.qa <- ggplot() +
  geom_spatraster(data = qa_sim[["2014-09-01"]]) +
  scale_fill_viridis_c("Environmental\nSimilarity\n",
                       option = "rocket") +
  coord_sf(expand = FALSE) +
  labs(title = "Qatar: Sep 2014") +
  theme_void() +
  theme(plot.title = element_text(size = 16, face = "bold"))


free(p.summ) | (p.br / p.qa)

ggsave("Tables_Figs/Figure S9.png", width = 9, height = 5, units = "in", dpi = 400)


######################
### Export results ###
######################

# Need to first "wrap" list of SpatRasters
qa_sim.w <- map(qa_sim, wrap)
br_sim.w <- map(br_sim, wrap)

saveRDS(qa_sim.w, "Data_products/Qatar_EnvironSim.rds")
saveRDS(br_sim.w, "Data_products/Brazil_EnvironSim.rds")
