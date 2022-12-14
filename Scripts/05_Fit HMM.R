

### Fit HMM to estimate behavioral states ###

library(tidyverse)
library(lubridate)
library(momentuHMM)
library(tictoc)
library(future)
library(furrr)
library(progressr)
library(trelliscopejs)

source('Scripts/helper functions.R')


#################
### Load data ###
#################


dat_gom <- read_csv("Processed_data/Processed_GoM_Cm_Tracks_SSM_2hr_foiegras.csv")
dat_br <- read_csv("Processed_data/Processed_Brazil_Cm_Tracks_SSM_2hr_foiegras.csv")
dat_qa <- read_csv("Processed_data/Processed_Qatar_Cm_Tracks_SSM_2hr_foiegras.csv")

glimpse(dat_gom)
glimpse(dat_br)
glimpse(dat_qa)


# Need to change class of 'id' col to character for GoM and Qatar data (to combine w/ Brazil dataset)
dat_gom$id <- as.character(dat_gom$id)
dat_qa$id <- as.character(dat_qa$id)


# Merge all datasets together with column to index their Region
dat <- list(GoM = dat_gom, Brazil = dat_br, Qatar = dat_qa) %>%
  bind_rows(.id = "Region")


# Wrangle data into proper format for {momentuHMM}
dat <- dat %>%
  rename(ID = id) %>%
  data.frame() %>%
  prepData(., type = "UTM", coordNames = c('x','y'))


# Remove any bouts that have large (i.e., 7-day) gaps
dat2 <- dat %>%
  filter(!is.na(bout))




#######################
### Fit 2-state HMM ###
#######################

# Viz SL and TA over time by ID
ggplot(dat, aes(date, step)) +
  geom_path(aes(group = ID, color = ID)) +
  theme_bw() +
  facet_trelliscope(~ID, nrow = 5, ncol = 5, scales = "free")


# Pre-define possible states to determine 'good' initial values
dat <- dat %>%
  mutate(phase = case_when(step > 1000 ~ 'Migratory',
                           TRUE ~ 'ARS'))


ggplot(dat, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


# Summarize mean and SD per phase
dat %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T))



### Set initial values for HMM

# initial step length distribution (natural scale parameters); gamma distribution
sum(dat$step == 0)  #Check to see if I need to account for zero mass
stepPar0 <- c(50, 1000, 100, 1000)  #(mu_1, mu_2, sd_1, sd_2)

# initial turning angle distribution (natural scale parameters); wrapped Cauchy distribution
anglePar0 <- c(-3.1, 0, 0.5, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)


Par0 <- list(step = stepPar0, angle = anglePar0)

# Remove 'phase' col
dat <- dat %>%
  dplyr::select(-phase)


### Fit HMM w/ random perturbations to determine global ML
set.seed(2022)
tic()
# fit_hmm <- run.HMMs(data = dat, K = 2, Par0 = Par0, state.names = c('ARS','Migratory'),
#                             niter = 20)
hmm.res <- fitHMM(data = dat, nbStates = 2,
                  Par0 = list(step = stepPar0, angle = anglePar0),
                  dist = list(step = "gamma", angle = "wrpcauchy"),
                  formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                  estAngleMean = list(angle=TRUE),
                  stateNames = c('ARS','Migratory'),
                  retryFits = 5
)
toc()  #took 10 min to run


fit_hmm

plot(fit_hmm)
plotStates(fit_hmm)
timeInStates(fit_hmm)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_hmm, ncores = 5)  #look decent for SL and TA, but Disp could be improved
