

### Fit M4 to estimate resident and migratory phases ###

library(tidyverse)
library(lubridate)
library(bayesmove)  #v2.0.0
library(tictoc)
library(future)
library(furrr)
library(progressr)
library(trelliscopejs)
library(MetBrewer)
library(rnaturalearth)
library(sf)
library(sfarrow)
library(plotly)
library(patchwork)

source('Scripts/helper functions.R')


#################
### Load data ###
#################


dat_gom <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_4hr_aniMotum.csv")
dat_br <- read_csv("Processed_data/Processed_Brazil_Cm_Tracks_SSM_4hr_aniMotum.csv")
dat_qa <- read_csv("Processed_data/Processed_Qatar_Cm_Tracks_SSM_4hr_aniMotum.csv")

glimpse(dat_gom)
glimpse(dat_br)
glimpse(dat_qa)


# Need to change class of 'id' col to character for GoM and Qatar data (to combine w/ Brazil dataset)
dat_gom$id <- as.character(dat_gom$id)
dat_qa$id <- as.character(dat_qa$id)


# Merge all datasets together with column to index their Region
dat <- list(GoM = dat_gom, Brazil = dat_br, Qatar = dat_qa) %>%
  bind_rows(.id = "Region")


# Map all tracks
world <- ne_countries(scale = 50, returnclass = 'sf') %>%
  filter(continent != "Antarctica")
dat.sf <- dat %>%
  drop_na(x, y) %>%
  st_as_sf(., coords = c('x','y'), crs = 3395) %>%
  st_transform(4326) %>%
  mutate(lon = st_coordinates(.)[,1],
         lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry()

# Generate color palette for tracks
set.seed(123)
col.pal <- sapply(names(MetPalettes), function(x) met.brewer(palette_name = x, n = 5, type = "continuous")) %>%
  as.vector() %>%
  sample(., size = n_distinct(dat.sf$id))

p.world <- ggplot() +
  geom_sf(data = world) +
  geom_path(data = dat.sf, aes(lon, lat, group = id, color = id)) +
  geom_rect(aes(xmin = -98, xmax = -77, ymin = 18, ymax = 32), color = "#009ACD", fill = NA, linewidth = 1) +  #inset rect for GoM
  geom_rect(aes(xmin = -49, xmax = -31, ymin = -26, ymax = -2), color = "#00CD00", fill = NA, linewidth = 1) +  #inset rect for Brazil
  geom_rect(aes(xmin = 50.25, xmax = 53, ymin = 23.8, ymax = 27.2), color = "#CD4F39", fill = NA, linewidth = 1) +  #inset rect for Qatar
  scale_color_manual(values = col.pal) +
  labs(x="",y="") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 16)) +
  coord_sf(xlim = c(-100,60), ylim = c(-50,50))

# ggsave("../../Conference Presentations/SERSTM 2023/world_map.png", width = 8, height = 6,
#        units = "in", dpi = 400)


# Map tracks by region
gom.sf <- st_read_parquet("Environ_data/GoM_land.parquet") %>%
  st_transform(4326)
br.sf <- st_read_parquet("Environ_data/Brazil_land.parquet") %>%
  st_transform(4326)
qa.sf <- st_read_parquet("Environ_data/Qatar_land.parquet") %>%
  st_transform(4326)

p.gom <- ggplot() +
  geom_sf(data = gom.sf) +
  geom_path(data = dat.sf, aes(lon, lat, group = id, color = id), linewidth = 0.75) +
  scale_color_manual(values = col.pal) +
  annotate(geom = "text", label = "Gulf of\nMexico", fontface = "italic", size = 8, x = -92, y = 25) +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "#009ACD", fill = NA, linewidth = 2),
        panel.background = element_rect(fill = "white")
        ) +
  coord_sf(xlim = c(-98,-77), ylim = c(18,32))
# ggsave("../../Conference Presentations/SERSTM 2023/gom_map.png", width = 8, height = 6,
#        units = "in", dpi = 400)


p.br <- ggplot() +
  geom_sf(data = br.sf) +
  geom_path(data = dat.sf, aes(lon, lat, group = id, color = id), linewidth = 0.75) +
  scale_color_manual(values = col.pal) +
  annotate(geom = "text", label = "Brazil", fontface = "italic", size = 8, x = -43, y = -8) +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "#00CD00", fill = NA, linewidth = 2),
        panel.background = element_rect(fill = "white")
        ) +
  coord_sf(xlim = c(-49,-31), ylim = c(-26,-2))
# ggsave("../../Conference Presentations/SERSTM 2023/br_map.png", width = 8, height = 6,
#        units = "in", dpi = 400)


p.qa <- ggplot() +
  geom_sf(data = qa.sf) +
  geom_path(data = dat.sf, aes(lon, lat, group = id, color = id), linewidth = 0.75) +
  scale_color_manual(values = col.pal) +
  annotate(geom = "text", label = "Qatar", fontface = "italic", size = 8, x = 52.25, y = 26.5) +
  labs(x="",y="") +
  theme_void() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "#CD4F39", fill = NA, linewidth = 2),
        panel.background = element_rect(fill = "white")
  ) +
  coord_sf(xlim = c(50.25,53), ylim = c(23.8,27.2))
# ggsave("../../Conference Presentations/SERSTM 2023/qa_map.png", width = 4, height = 6,
#        units = "in", dpi = 400)




### Create composite map ###

p.world +
  inset_element(p.gom, left = 0.23, bottom = 0.64, right = 0.65, top = 0.97, align_to = "plot") +
  inset_element(p.br, left = 0.1, bottom = 0.1, right = 0.35, top = 0.6, align_to = "plot") +
  inset_element(p.qa, left = 0.75, bottom = 0.15, right = 0.95, top = 0.65, align_to = "plot")

# ggsave("Tables_Figs/Figure_1.png", width = 10, height = 6, units = "in", dpi = 400)




##############################
### Prepare data for model ###
##############################

# Wrangle data into proper format for {momentuHMM}
dat2 <- dat %>%
  prep_data(dat = ., coord.names = c('x','y'), id = "id")


# Remove any bouts that have large (i.e., 3-day) gaps
dat3 <- dat2 %>%
  drop_na(x, y)


# Since we don't need to filter out obs at other time intervals, we still need to add required variables to data.frame
dat3 <- dat3 %>%
  group_by(id) %>%  #need to number rows separately for each ID
  mutate(time1 = 1:n(),
         obs = 1:n()) %>%
  ungroup()

#verify that it worked properly
dat3 %>%
  dplyr::select(id, date, time1, obs) %>%   #select only a few cols since tibble hides time1 and obs
  split(.$id) %>%
  head()


# Create displacement variable from NSD and convert to km (from m)
dat3$disp <- sqrt(dat3$NSD) / 1000

# Convert step length from m to km
dat3$step <- dat3$step / 1000



### Discretize data stream for models

# Viz density plots of step length and displacement
ggplot(dat3) +
  geom_density(aes(step), fill = "cadetblue") +
  labs(x = "Step Length (km)") +
  theme_bw()

ggplot(dat3) +
  geom_density(aes(disp), fill = "goldenrod") +
  labs(x = "Displacement (km)") +
  theme_bw()


# Define bin limits and number of bins (must be positive, but no upper bound)
step.bin.lims <- c(0, 5, 10, 20, 40, max(dat3$step, na.rm = TRUE))  #5 bins

disp.bin.lims <- c(0, 50, 100, max(dat3$disp, na.rm = TRUE))  #3 bins

step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc <- discrete_move_var(dat3,
                              lims = list(step.bin.lims, disp.bin.lims),
                              varIn = c("step","disp"),
                              varOut = c("SL","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()






##########################
### Segment the tracks ###
##########################

# Convert data to list by ID
dat.list <- dat.disc %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub <- map(dat.list,
                   subset,
                   select = c(id, SL, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha <- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs <- 80000  # number of iterations for Gibbs sampler
nbins <- c(max(dat.disc$SL, na.rm = TRUE), max(dat.disc$Disp))  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg <- segment_behavior(data = dat.list.sub, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 3 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg, type = "LML")  #appears to have converged for each track



# Determine MAP for selecting breakpoints
MAP.est <- get_MAP(dat = dat.res.seg$LML, nburn = ngibbs/2)
brkpts <- get_breakpts(dat = dat.res.seg$brkpts, MAP.est = MAP.est)

# How many breakpoints estimated per ID?
apply(brkpts[,-1], 1, function(x) length(purrr::discard(x, is.na)))
brkpts


# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list, as_date = TRUE, var_names = c("step","disp"),
                 var_labels = c("Step Length (km)","Displacement (km)"),
                 brkpts = brkpts)

plot_breakpoints(data = dat.list, as_date = FALSE, var_names = c("SL","Disp"),
                 var_labels = c("Step Length (km)","Displacement (km)"),
                 brkpts = brkpts)


# Assign track segments to each ID
dat.seg <- assign_tseg(dat = dat.list, brkpts = brkpts)

head(dat.seg)





###############################################
### Cluster segments into behavioral states ###
###############################################

#Select only id, tseg, and discretized data streams
dat.seg2 <- dat.seg[,c("id","tseg","SL","Disp")]

#Summarize observations by track segment
obs <- summarize_tsegs(dat = dat.seg2, nbins = nbins)
obs


set.seed(123)

# Prepare for Gibbs sampler
ngibbs <- 10000  #number of MCMC iterations for Gibbs sampler
nburn <- ngibbs/2  #number of iterations for burn-in
nmaxclust <- 4  #max number of groups to test for
ndata.types <- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1 <- 0.1
alpha <- 0.1

# Run LDA model
dat.res.segclust<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                    ngibbs = ngibbs, nmaxclust = nmaxclust,
                                    nburn = nburn, ndata.types = ndata.types)
# takes 2.5 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim <- extract_prop(res = dat.res.segclust, ngibbs = ngibbs, nburn = nburn,
                           nmaxclust = nmaxclust)

theta.estim_df <- theta.estim %>%
  as.data.frame() %>%
  pivot_longer(., cols = everything(), names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior) <- 1:nmaxclust

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


#Calculate mean proportions per behavior
(theta.means <- round(colMeans(theta.estim), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res.seg <- get_behav_hist(dat = dat.res.segclust, nburn = nburn, ngibbs = ngibbs,
                               nmaxclust = nmaxclust,
                               var.names = c("Step Length","Displacement"))

# Add bin lim range to each label
step.lims <- data.frame(bin.vals = cut(dat3$step, round(step.bin.lims, 2)) %>%
                          levels(),
                        bin = 1:(length(step.bin.lims) - 1),
                        var = "Step Length")

disp.lims <- data.frame(bin.vals = cut(dat3$disp, round(disp.bin.lims, 2)) %>%
                          levels(),
                        bin = 1:(length(disp.bin.lims) - 1),
                        var = "Displacement")
lims <- rbind(step.lims, disp.lims)

behav.res.seg <- left_join(behav.res.seg, lims, by = c('var','bin'))
behav.res.seg$bin.vals <- factor(behav.res.seg$bin.vals, levels = unique(behav.res.seg$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.seg, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = MetPalettes$Egypt[[1]], guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")


# Merge Resident states together into one
theta.estim[,1] <- theta.estim[,1] + theta.estim[,3] + theta.estim[,4]
theta.estim <- theta.estim[,1:2]



#Reformat proportion estimates for all track segments
theta.estim.long <- expand_behavior(dat = dat.seg, theta.estim = theta.estim, obs = obs, nbehav = 2,
                                   behav.names = c("Resident","Migratory"),
                                   behav.order = 1:2)

#Plot results
ggplot(theta.estim.long %>%
         dplyr::select(id, date, prop, behavior), aes(x=date, y=prop, fill = behavior)) +
  geom_area(color = "black", linewidth = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_trelliscope(~id, scales = "free_x", ncol = 3, nrow = 3)




#### Assign states to segments and map ####
world <- ne_countries(scale = 50, returnclass = 'sf') %>%
  filter(continent != 'Antarctica') %>%
  st_transform(3395)

# Convert segmented dataset into list
dat.seg.list <- df_to_list(dat = dat.seg, ind = "id")

# Merge results with original data
dat.out <- assign_behavior(dat.orig = dat3,
                          dat.seg.list = dat.seg.list,
                          theta.estim.long = theta.estim.long,
                          behav.names = levels(theta.estim.long$behavior))

# Map dominant behavior for all IDs
plotly::ggplotly(
  ggplot() +
    geom_sf(data = world) +
    geom_path(data = dat.out, aes(x=x, y=y, group = id), color="grey60", linewidth=0.25) +
    geom_point(data = dat.out, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.out %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.out %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_viridis_d("Behavior", na.value = "grey50") +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14)))
)



# Split by Region for export
gom.dat.out <- dat.out %>%
  filter(Region == "GoM")
br.dat.out <- dat.out %>%
  filter(Region == "Brazil")
qa.dat.out <- dat.out %>%
  filter(Region == "Qatar")


###############################
### Export annotated tracks ###
###############################

write_csv(gom.dat.out, "Processed_data/GoM_Cm_Tracks_behav.csv")
write_csv(br.dat.out, "Processed_data/Brazil_Cm_Tracks_behav.csv")
write_csv(qa.dat.out, "Processed_data/Qatar_Cm_Tracks_behav.csv")
