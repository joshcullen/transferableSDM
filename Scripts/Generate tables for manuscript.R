
#### Create summary tables for manuscript ####

library(tidyverse)
library(lubridate)


#################
### Load data ###
#################

dat <- read_csv("Raw_data/Master Sat Tag Dataset.csv")





### Summarize number of obs post-SSM that includes Region, ID, Age Class, Start Date, End Date, and Duration per turtle
dat.summary <- dat %>%
  group_by(Region, Ptt, Age) %>%
  mutate(Start = as_date(first(Date)),
         End = as_date(last(Date)),
         Duration = End - Start,
         N = n()) %>%
  dplyr::select(Region, Ptt, Age, Start, End, Duration, N) %>%
  distinct() %>%
  ungroup() %>%
  # mutate(across(Start:End, \(x) as_date(x))) %>%
  mutate(across(Duration, \(x) as.numeric(round(x,0)))) %>%
  mutate(across(Region, \(x) factor(x, levels = c('GoM','Brazil','Qatar')))) %>%
  arrange(Region, Start)

write_csv(dat.summary, "Tables_Figs/Table S1.csv")




### Summarize number of turtles tracked, number of observations, breakdown by age class, and range of study years per region
dat.summary.general <- dat %>%
  group_by(Region) %>%
  arrange(Date) %>%
  summarize(N_id = n_distinct(Ptt),
          N_juv = n_distinct(Ptt[Age == 'Juv']),
          N_adult = n_distinct(Ptt[Age == 'Adult']),
          N_obs = n(),
          Start = year(first(Date)),
          End = year(last(Date))
  ) %>%
  arrange(desc(N_id))

write_csv(dat.summary.general, "Tables_Figs/Table 1.csv")
