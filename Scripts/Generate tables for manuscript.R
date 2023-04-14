
#### Create summary tables for manuscript ####

library(tidyverse)


#################
### Load data ###
#################

dat <- read_csv("Processed_data/Prefiltered_Cm_Tracks.csv")





### Summarize number of obs post-SSM that includes Region, ID, Age Class, Start Date, End Date, and Duration per turtle
dat.summary <- dat %>%
  group_by(Region, Ptt, Age) %>%
  mutate(Start = first(date),
         End = last(date),
         Duration = End - Start,
         N = n()) %>%
  dplyr::select(Region, Ptt, Age, Start, End, Duration, N) %>%
  distinct() %>%
  ungroup() %>%
  mutate(across(Start:End, \(x) as_date(x))) %>%
  mutate(across(Duration, \(x) as.numeric(round(x,0)))) %>%
  mutate(across(Region, \(x) factor(x, levels = c('GoM','Brazil','Qatar')))) %>%
  arrange(Region, Start) %>%
  print(n=98)

write_csv(dat.summary, "Tables_Figs/Table S1.csv")




### Summarize number of turtles tracked, number of observations, breakdown by age class, and range of study years per region
dat.summary.general <- dat %>%
  group_by(Region) %>%
  summarize(N_id = n_distinct(Ptt),
          N_obs = n(),
          N_juv = n_distinct(Ptt[Age == 'Juv']),
          N_adult = n_distinct(Ptt[Age == 'Adult']),
          Start = year(first(date)),
          End = year(last(date))
  ) %>%
  arrange(desc(N_id))

write_csv(dat.summary.general, "Tables_Figs/Table 1.csv")
