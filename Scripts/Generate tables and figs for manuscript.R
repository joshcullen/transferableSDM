
#### Create summary tables and extra figs for manuscript ####

library(tidyverse)
library(lubridate)
library(MetBrewer)
library(patchwork)


#################
### Load data ###
#################

dat <- read_csv("Raw_data/Master Sat Tag Dataset.csv")
dat.meta <- read_csv("Raw_data/Turtle tag metadata.csv") %>%
  dplyr::select(Ptt, CCL_cm)


boyce.alg <- read_csv("Data_products/boyce_alg_results.csv")
boyce.scale <- read_csv("Data_products/boyce_scale_results.csv") |>
  mutate(across(scale, \(x) factor(x, levels = c(5, 10, 20, 40))))
boyce.age <- read_csv("Data_products/boyce_age_results.csv")




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
  arrange(Region, Start) %>%
  left_join(., dat.meta, by = 'Ptt') %>%
  rename(CCL = CCL_cm, ID = Ptt, `Life stage` = Age) %>%
  dplyr::relocate(CCL, .after = ID) %>%
  mutate(CCL = round(CCL, 1))

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
  mutate(N_adult = case_when(Region == "Brazil" ~ N_adult - 1,  #one of adults tracked twice in Brazil
                             TRUE ~ N_adult),
         N_id = N_juv + N_adult) %>%
  arrange(desc(N_id))

write_csv(dat.summary.general, "Tables_Figs/Table 1.csv")




### Create composite figure of Boyce Index plots

boyce.alg.mean <- boyce.alg %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

boyce.scale.mean <- boyce.scale %>%
  group_by(scale, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()

boyce.age.mean <- boyce.age %>%
  group_by(Method, Region) %>%
  summarize(mean = mean(cor, na.rm = TRUE)) %>%
  ungroup()


p.alg <- ggplot(data = boyce.alg, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(Method, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.alg.mean, aes(x = Region, y = mean, group = Method),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_met_d('Egypt') +
  scale_fill_met_d('Egypt') +
  # scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank()) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))


p.scale <- ggplot(data = boyce.scale, aes(Region, cor)) +
  geom_point(aes(fill = scale), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(scale, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.scale.mean, aes(x = Region, y = mean, group = scale),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_viridis_d(option = "mako", guide = "none", begin = 0.5, end = 0.95, direction = -1) +
  scale_fill_viridis_d("Scale (km)", option = "mako", begin = 0.5, end = 0.95, direction = -1) +
  # scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "Boyce Index") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank()) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))


p.age <- ggplot(data = boyce.age, aes(Region, cor)) +
  geom_point(aes(fill = Method), pch = 21, alpha = 0.7, size = 4,
             position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(group = interaction(Method, Region)), fill = "transparent",
               position = position_dodge(width = 0.75),
               outlier.shape = NA, width = 0.6, size = 0.75) +
  geom_point(data = boyce.age.mean, aes(x = Region, y = mean, group = Method),
             size = 4, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = met.brewer("Isfahan1")[c(3,6)], guide = "none") +
  scale_fill_manual(values = met.brewer("Isfahan1")[c(3,6)], labels = c("Age", "No Age")) +
  scale_x_discrete(labels = c("Brazil (all)", "Brazil (main)", "Qatar")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  lims(y = c(-1,1)) +
  labs(x="", y = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 24)) +
  guides(color = "none",
         fill = guide_legend(override.aes = list(alpha = 1)))


plot_spacer() + p.alg + p.scale + p.age +
  plot_layout(nrow = 4, heights = c(0.1, 1, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "(",
                  tag_suffix = ")") &
  theme(plot.tag.position = c(0.05, 1),
        plot.tag = element_text(size = 18, face = "bold", hjust = 0.8, vjust = -0.5))

ggsave("Tables_Figs/boyce.png", width = 7, height = 9, units = "in", dpi = 400)
