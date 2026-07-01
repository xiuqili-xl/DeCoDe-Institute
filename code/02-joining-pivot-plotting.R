# Overview -----
# DeCoDe Institute 2026-07-01 Training
# Joining, pivoting, and plotting in R

# Continuing from `01-intro-R-tidyverse.R`


# Load packages ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(ggbeeswarm)
library(here)


# Load data ----
data(PHENO)

group_assignment <- read_csv(here("outputs_csv", "group_assignment.csv"))
pheno_bw_key_timepoint <- read_csv(here("outputs_csv", "pheno_bw_key_timepoints.csv")) %>%
  mutate(pid = as.character(pid))
pheno_weights_terminal <- read_csv(here("outputs_csv", "pheno_weights_terminal.csv")) %>%
  mutate(pid = as.character(pid))
pheno_vo2 <- read_csv(here("outputs_csv", "pheno_vo2.csv")) %>%
  mutate(pid = as.character(pid))


# Distribution of VO2 max1  ----
glimpse(pheno_vo2)

ggplot(data = pheno_vo2, 
       mapping = aes(x = visit1_vo2_max)) +
  geom_histogram(binwidth = 2) +
  theme_bw() +
  labs(title = "Histogram of VO2 max from visit 1")


ggplot(data = pheno_vo2, 
       mapping = aes(x = visit1_vo2_max)) +
  geom_density() +
  theme_bw() +
  labs(title = "Density plot of VO2 max from visit 1")


# Box plot of sex vs VO2 max 1 ----
ggplot(data = pheno_vo2,
       mapping = aes(x = sex, y = visit1_vo2_max, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Boxplot of VO2 max 1")


# Scatterplot of VO max 1 and VO max 2 ----

ggplot(data = pheno_vo2,
       mapping = aes(x = visit1_vo2_max, y = visit2_vo2_max, 
                     shape = sex, color = group)) +
  geom_point() +
  theme_bw() +
  labs(title = "Scatterplot of VO2 max 1 and VO2 max 2")



# Change in VO2 max (Schenk 2024, Figure 3) ----
PHENO %>% 
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, 
         vo2_day1 = vo2.max.test.d_vo2_1, vo2_day2 = vo2.max.test.d_vo2_2, 
         familiarization = familiarization.d_visit, training_day1 = training.day1date, 
         training_day40 = training.day40date, sacrifice = key.d_sacrifice) %>%
  mutate(across(c(vo2_day1, vo2_day2, familiarization, training_day1, training_day40, sacrifice), dmy)) %>%
  mutate(training_from_fam = training_day1 - familiarization, 
         vo2_day1_from_fam = vo2_day1 - familiarization,
         vo2_day1_from_training_day1 = vo2_day1 - training_day1,
         vo2_day2_from_training_day1 = vo2_day2 - training_day1,
         sacrifice_from_training_day40 = sacrifice - training_day40,
         sacrifice_from_vo2_day2 = sacrifice - vo2_day2) %>%
  distinct() %>%
  view()
# so.... vo2_day 1 is ~ a week after familiarization and ~ a week before first day of training
# sacrifice is two days after training_day40, vo2_day2 is 6 days before sacrifice
# matches what's report in Schenk 2024 (Figure 1)

## Graph VO2 weight change ----
pheno_vo2_change <- pheno_vo2 %>%
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, ends_with("vo2_max")) %>%
  arrange(sex, group, visit1_vo2_max) %>%
  group_by(group, sex) %>%
  mutate(vo2_change = visit2_vo2_max - visit1_vo2_max, 
         vo2_change_direction = if_else(vo2_change > 0, "increase", "decrease"),
         graph_order = row_number(),
         group = factor(group, levels = c("control", "8w"))) %>%
  ungroup()


ggplot(data = pheno_vo2_change) +
  geom_segment(mapping = aes(x = graph_order, y = visit1_vo2_max, 
                             xend = graph_order, yend = visit2_vo2_max,
                             color = vo2_change_direction),
               arrow = arrow(length = unit(3, "pt"))) +
  scale_color_manual(values = c("increase" = "darkred", "decrease" = "blue4")) +
  facet_grid(sex ~ group, scales = "free", switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside") +
  labs(title = "Change in VO2 max before and after training", 
       subtitle = "(`vo2.max.test.vo2_max_1` vs `vo2.max.test.vo2_max_2`)", 
       x = "", y = "VO2(mL/kg/min)") 

# based on the number and comparing it to Schnek et al 2024 (Figure 3), 
# it looks like these have already be normalize by animal body weight

ggsave(path = here("graphs"), filename = "pheno_vo2_pre_vs_post.png",
       width = 5, height = 7, dpi = 300, units = "in")




# Body Composition (Schenk 2024, Figure 5) ----
pheno_nmr <- PHENO %>%
  select(pid, sex, group, 
         visit1_date = nmr.testing.d_nmr_1, visit2_date = nmr.testing.d_nmr_2, 
         visit1_weight = nmr.testing.nmr_weight_1, visit2_weight = nmr.testing.nmr_weight_2, 
         visit1_fat = nmr.testing.nmr_fat_1, visit2_fat = nmr.testing.nmr_fat_2, 
         visit1_lean = nmr.testing.nmr_lean_1, visit2_lean = nmr.testing.nmr_lean_2, 
         visit1_fluid = nmr.testing.nmr_fluid_1, visit2_fluid = nmr.testing.nmr_fluid_2,
         visit2_comments = nmr.testing.nmr_comments_2) %>%
  distinct() %>%
  arrange(sex, group, pid)

write_csv(pheno_nmr, here("outputs_csv", "pheno_nmr.csv"))


## sanity check: when where nmr testing carried out?
PHENO %>% 
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, 
         nmr_day1 = nmr.testing.d_nmr_1, nmr_day2 = nmr.testing.d_nmr_2, 
         familiarization = familiarization.d_visit, training_day1 = training.day1date, 
         training_day40 = training.day40date, sacrifice = key.d_sacrifice) %>%
  mutate(across(c(nmr_day1, nmr_day2, familiarization, training_day1, training_day40, sacrifice), dmy)) %>%
  mutate(nmr_day1_from_fam = nmr_day1 - familiarization,
         nmr_day1_from_training_day1 = nmr_day1 - training_day1,
         nmr_day2_from_training_day1 = nmr_day2 - training_day1,
         sacrifice_from_nmr_day2 = sacrifice - nmr_day2) %>%
  distinct() %>%
  view()

# nmr day1 was taken 1 day after familiarization, about 2 weeks before start of training
# nmr day2 was taken 5 days before training wraps, so in the middle of last week of training
# matches reporting in Schenk et al 2024


## graph NMR weight change ----
pheno_nmr_weight_change <- pheno_nmr %>% 
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, ends_with("weight")) %>%
  arrange(sex, group, visit1_weight) %>%
  group_by(group, sex) %>%
  mutate(wt_change = visit2_weight - visit1_weight, 
         wt_change_direction = if_else(wt_change > 0, "increase", "decrease"),
         graph_order = row_number(),
         group = factor(group, levels = c("control", "8w"))) %>%
  ungroup() 


ggplot(data = pheno_nmr_weight_change) +
  geom_segment(mapping = aes(x = graph_order, y = visit1_weight, 
                             xend = graph_order, yend = visit2_weight,
                             color = wt_change_direction),
               arrow = arrow(length = unit(3, "pt"))) +
  scale_color_manual(values = c("increase" = "darkred", "decrease" = "blue4")) +
  facet_grid(sex ~ group, scales = "free", switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside") +
  labs(title = "Change in body weight before and after training", 
       subtitle = "(`nmr.testing.nmr_weight_1` vs `nmr.testing.nmr_weight_2`)", 
       x = "", y = "Body mass(g)") 

ggsave(path = here("graphs"), filename = "pheno_nmr_weight_pre_vs_post.png",
       width = 5, height = 7, dpi = 300, units = "in")



## Graph LEAN percentage change ----
pheno_nmr_lean_change <- pheno_nmr %>% 
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, ends_with("lean")) %>%
  arrange(sex, group, visit1_lean) %>%
  group_by(group, sex) %>%
  mutate(lean_change = visit2_lean - visit1_lean, 
         lean_change_direction = if_else(lean_change > 0, "increase", "decrease"),
         graph_order = row_number(),
         group = factor(group, levels = c("control", "8w"))) %>%
  ungroup() 


ggplot(data = pheno_nmr_lean_change) +
  geom_segment(mapping = aes(x = graph_order, y = visit1_lean, 
                             xend = graph_order, yend = visit2_lean,
                             color = lean_change_direction),
               arrow = arrow(length = unit(3, "pt"))) +
  scale_color_manual(values = c("increase" = "darkred", "decrease" = "blue4")) +
  facet_grid(sex ~ group, scales = "free_x", switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside") +
  labs(title = "Change in percent lean mass before and after training", 
       subtitle = "(`nmr.testing.nmr_lean_1` vs `nmr.testing.nmr_lean_2`)", 
       x = "", y = "NMR % lean mass") 

ggsave(path = here("graphs"), filename = "pheno_nmr_lean_percent_pre_vs_post.png",
       width = 5, height = 7, dpi = 300, units = "in")

## overall matches Schenk Figure S5



## Graph LEAN percentage change ----
pheno_nmr_fat_change <- pheno_nmr %>% 
  filter(group %in% c("control", "8w")) %>%
  select(pid, group, sex, ends_with("fat")) %>%
  arrange(sex, group, visit1_fat) %>%
  group_by(group, sex) %>%
  mutate(fat_change = visit2_fat - visit1_fat, 
         fat_change_direction = if_else(fat_change > 0, "increase", "decrease"),
         graph_order = row_number(),
         group = factor(group, levels = c("control", "8w"))) %>%
  ungroup() 


ggplot(data = pheno_nmr_fat_change) +
  geom_segment(mapping = aes(x = graph_order, y = visit1_fat, 
                             xend = graph_order, yend = visit2_fat,
                             color = fat_change_direction),
               arrow = arrow(length = unit(3, "pt"))) +
  scale_color_manual(values = c("increase" = "darkred", "decrease" = "blue4")) +
  facet_grid(sex ~ group, scales = "free_x", switch = "y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside") +
  labs(title = "Change in percent fat mass before and after training", 
       subtitle = "(`nmr.testing.nmr_fat_1` vs `nmr.testing.nmr_fat_2`)", 
       x = "", y = "NMR % fat mass") 

ggsave(path = here("graphs"), filename = "pheno_nmr_fat_percent_pre_vs_post.png",
       width = 5, height = 7, dpi = 300, units = "in")
## overall matches Schenk Figure S5



## Sanity chack: What does fat, lean, and fluid add up to?
pheno_nmr %>%
  filter(group %in% c("control", "8w")) %>%
  mutate(visit1_sum = visit1_fat + visit1_lean + visit1_fluid,
         visit2_sum = visit2_fat + visit2_lean + visit2_fluid) %>%
  view()
# the sums are coming up to around 80%, which is expected
# the result are undetectable substances + systematic errors
# undetectable substances include bone mineral content, claws, hair, etc 
# that do not contribute to NMR signal detection



# Training data -----
pheno_training_info <- PHENO %>%
  select(pid, sex, group, starts_with("training.day")) %>%
  select(pid, sex, group, ends_with("date"), ends_with("time"), 
         ends_with("treadmillspeed"), ends_with("treadmillincline"), ends_with("treadmill"),
         ends_with("score"), ends_with("comments")) %>%
  distinct() %>%
  arrange(sex, group, pid)

write_csv(pheno_training_info, here("outputs_csv", "pheno_training_info.csv"))


pheno_training_data <- PHENO %>%
  select(pid, sex, group, starts_with("training.day")) %>%
  select(pid, sex, group, ends_with("weight"), ends_with("posttrainlact")) %>%
  distinct() %>%
  arrange(sex, group, pid)

write_csv(pheno_training_data, here("outputs_csv", "pheno_training_data.csv"))


# Graph weight change during training -----
pheno_training_wt_summary <- pheno_training_data %>%
  select(pid, sex, group, ends_with("weight")) %>%
  pivot_longer(cols = -c(pid, sex, group), names_to = "training_day", values_to = "weight",
               values_drop_na = TRUE) %>%
  mutate(training_day = str_remove(training_day, "training.day"),
         training_day = str_remove(training_day, "_weight"),
         training_day = as.numeric(training_day)) %>%
  group_by(sex, group, training_day) %>%
  summarise(weight_mean = mean(weight),
            weight_sd = sd(weight),
            weight_se = weight_sd / sqrt(n())) %>%
  ungroup()


ggplot(data = pheno_training_wt_summary %>% 
         filter(group %in% c("control", "8w")) %>%
         mutate(group = forcats::fct_rev(group)),
       mapping = aes(x = training_day, y = weight_mean, color = group, group = group)) +
  geom_errorbar(aes(ymin = weight_mean - weight_se, ymax = weight_mean + weight_se), 
                width = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(0, 40), expand = 0.01) +
  facet_wrap(~sex, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Body weight of mice during 40-day training",
       subtitle = "(mean +/- standard error)", 
       x = "Days of training", y = "Body mass (g)")

ggsave(path = here("graphs"), filename = "pheno_training_weight.png",
       width = 6, height = 6, dpi = 300, units = "in")



# Graph blood lactate change during training -----
pheno_training_lact_summary <- pheno_training_data %>%
  select(pid, sex, group, ends_with("lact")) %>%
  pivot_longer(cols = -c(pid, sex, group), names_to = "training_day", values_to = "lactate",
               values_drop_na = TRUE) %>%
  mutate(training_day = str_remove(training_day, "training.day"),
         training_day = str_remove(training_day, "_posttrainlact"),
         training_day = as.numeric(training_day)) %>%
  group_by(sex, group, training_day) %>%
  summarise(lactate_mean = mean(lactate),
            lactate_sd = sd(lactate),
            lactate_se = lactate_sd / sqrt(n())) %>%
  ungroup()


ggplot(data = pheno_training_lact_summary %>% 
         filter(group %in% c("control", "8w")) %>%
         mutate(group = forcats::fct_rev(group)),
       mapping = aes(x = training_day, y = lactate_mean, color = sex, group = sex)) +
  geom_errorbar(aes(ymin = lactate_mean - lactate_se, ymax = lactate_mean + lactate_se), 
                width = 0.5) +
  geom_point() +
  geom_line() +
  scale_x_continuous(limits = c(0, 40), expand = 0.01) +
  facet_wrap(~sex, ncol = 1) +
  theme_bw() +
  labs(title = "Post training lactate of mice during 40-day training",
       subtitle = "(only for rats doing endurance training; mean +/- standard error)", 
       x = "Days of training", y = "Blood lactate (?unit)")


ggsave(path = here("graphs"), filename = "pheno_training_lactate.png",
       width = 6, height = 6, dpi = 300, units = "in")
