# Overview -----
# DeCoDe Institute 2026-06-30 Training
# An Introduction to R and Data Wrangling with the Tidyverse

# Course materials https://github.com/CFDETrainingCenter/decode-institute-2026 (CC0 license)
# Downloaded the Rmd file from `06_intro_to_R` and stored in `code_markdown`, with a couple of modifications
# (1) changed file path for images
# (2) added an initial code chunk so no future chunks will excut
# (3) turned output to markdown notebook for easy viewing (no need to knit)


# MoTrPAC dataset ----
# pronounced: motor-pack
# Project website: https://motrpac.org/
# Data hub: https://motrpac-data.org/

# MotrpacRatTraining6moData package: https://motrpac.github.io/MotrpacRatTraining6moData/

# install this package
#options(timeout=1e5)                     # extend the timeout
#devtools::install_github("MoTrPAC/MotrpacRatTraining6moData")



# Load packages ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(ggbeeswarm)
library(here)



# Explore data available through MotrpacRatTraining6moData ----
# list available objects
data(package = "MotrpacRatTraining6moData")

# get the documentation of a give data object
?PHENO         # Phenotypic data

# load a data object into the environment
data(PHENO)


# Explore the PHENO data ----
glimpse(PHENO)

nrow(PHENO)      # 6156 rows

PHENO %>%
  select(pid, labelid, viallabel) %>%
  distinct() %>%
  nrow()
# so, each entry is uniquly identified by pid, labelid, and viallabel


length(unique(PHENO$pid))              # a total of 147 animals
unique(PHENO$group)
unique(PHENO$sex)

# figure out group assignment
group_assignment <- PHENO %>%
  count(pid, key.d_arrive, key.d_sacrifice, key.intervention, group, sex, name = "sample_count") %>%
  mutate(group = factor(group, levels = c("1w", "2w", "4w", "8w", "control"))) %>%
  arrange(sex, group, pid)

write_csv(group_assignment, here("outputs_csv", "group_assignment.csv"))



# Explore terminal weight in PHENO ----
pheno_terminal_weight <- PHENO %>%
  select(pid, sex, group, starts_with("terminal.weight")) %>%
  distinct() %>%
  arrange(sex, group, pid)

write_csv(pheno_terminal_weight, here("outputs_csv", "pheno_weights_terminal.csv"))

pheno_terminal_weight %>%
  group_by(sex, group) %>%
  summarise(terminal_bw_mean = mean(terminal.weight.bw),
            terminal_bw_sd = sd(terminal.weight.bw)) %>%
  ungroup()

# only drawing 8 wks vs control so rats are age-matched
pheno_terminal_weight %>% 
  filter(group %in% c("control", "8w")) %>%
  mutate(group = factor(group, levels = c("control", "8w"))) %>%
  ggplot(mapping = aes(x = group, y = terminal.weight.bw)) +
  geom_violin(mapping = aes(color = group),
              quantiles = 0.5, quantile.linetype = "solid" ) +
  geom_beeswarm(mapping = aes(color = group), shape = 21, size = 2, cex = 4) +
  scale_y_continuous(limits = c(150, 400), expand = 0.01) +
  facet_wrap(~sex) +
  theme_bw() +
  labs(title = "Terminal body weight of mice (8wks vs control)",
       subtitle = "(`terminal.weight.bw`)") +
  theme(legend.position = "none")


ggsave(path = here("graphs"), filename = "pheno_bw_terminal.png",
       width = 6, height = 4, dpi = 300, units = "in")



# Weight change (Figure 5 of Schneck) -----
pheno_bw_key_timepoints <- PHENO %>%
  select(pid, sex, group, 
         registration.d_arrive, registration.weight, 
         familiarization.d_visit, familiarization.weight, 
         key.d_sacrifice, terminal.weight.bw) %>% 
  distinct() %>%
  arrange(sex, group, pid)

write_csv(pheno_bw_key_timepoints, here("outputs_csv", "pheno_bw_key_timepoints.csv"))


pheno_bw_change <- pheno_bw_key_timepoints %>%
  select(pid, group, sex, start.bw = familiarization.weight, end.bw = terminal.weight.bw) %>%
  filter(group %in% c("control", "8w")) %>%
  arrange(sex, group, pid) %>%
  mutate(change.bw = end.bw - start.bw,
         change.label = if_else(change.bw > 0, "increase", "decrease")) %>%
  arrange(sex, start.bw, change.bw) %>%
  group_by(sex, group) %>%
  mutate(graph_order = row_number(),
         group = factor(group, levels = c("control", "8w"))) %>%
  ungroup()

ggplot(data = pheno_bw_change) +
  geom_segment(mapping = aes(x = graph_order, y = start.bw, 
                             xend = graph_order, yend = end.bw,
                             color = change.label),
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
  labs(title = "Change in bw pre- and post-training", 
       subtitle = "(`familiarization.weight` vs `terminal.weight.bw`)", 
       x = "", y = "body weight(g)") 

ggsave(path = here("graphs"), filename = "pheno_bw_familiarization_vs_terminal.png",
       width = 4, height = 6, dpi = 300, units = "in")



# Plasma samples -----
unique(PHENO$tissue)         # note, there is an NA there

# examine entries where tissue is NA
pheno_tissue_na <- PHENO %>%
  filter(is.na(tissue)) %>%
  select(pid, labelid, viallabel, tissue, tissue_description, tissue_code_no) 

pheno_tissue_na %>%
  count(tissue_description)
# so... all NA have tissue description of Lateral Gastrocnemius (Histology)

unique(PHENO$tissue_description)       # all options for tissue_description (controlled vocab)

PHENO %>%
  # any entries with lateral gastrocnemius as description, but also has a tissue assignment? 
  filter(!is.na(tissue), tissue_description == "Lateral Gastrocnemius (Histology)") %>%
  nrow()
# so... all is.na(tissue) are for Lateral Gastrocnemius (Histology)

sample_plasma <- PHENO %>% 
  filter(tissue == "PLASMA") %>%
  select(pid, labelid, viallabel, tissue, tissue_description, tissue_code_no)
  
nrow(sample_plasma)
nrow(sample_plasma %>% distinct())



# VO2 max -----
pheno_vo2 <- PHENO %>% 
  select(pid, sex, group, 
         visit1_date = vo2.max.test.d_vo2_1, visit2_date = vo2.max.test.d_vo2_2,
         visit1_blactate_begin = vo2.max.test.blactate_begin_1, visit2_blactate_begin = vo2.max.test.blactate_begin_2,
         visit1_vo2_max = vo2.max.test.vo2_max_1, visit2_vo2_max = vo2.max.test.vo2_max_2,
         visit1_vco2_max = vo2.max.test.vco2_max_1, visit2_vco2_max = vo2.max.test.vco2_max_2,
         visit1_rer_max = vo2.max.test.rer_max_1, visit2_rer_max = vo2.max.test.rer_max_2,
         visit1_speed_max = vo2.max.test.speed_max_1, visit2_speed_max = vo2.max.test.speed_max_2,
         visit1_blactate_end = vo2.max.test.blactate_end_1, visit2_blactate_end = vo2.max.test.blactate_end_2,
         visit1_comments = vo2.max.test.vo2_comments_1, visit2_comments = vo2.max.test.vo2_comments_2) %>%
  distinct() %>%
  arrange(sex, group, pid)

nrow(pheno_vo2)

write_csv(pheno_vo2, here("outputs_csv", "pheno_vo2.csv"))


pheno_vo2 %>% 
  filter(group %in% c("control", "8w")) %>%
  mutate(group = factor(group, levels = c("control", "8w"))) %>%
  ggplot(mapping = aes(x = group, y = visit1_vo2_max)) +
  geom_violin(mapping = aes(color = group),
              quantiles = 0.5, quantile.linetype = "solid" ) +
  geom_beeswarm(mapping = aes(color = group), shape = 21, size = 2, cex = 4) +
  scale_y_continuous(limits = c(50, 90), expand = 0.01) +
  facet_wrap(~sex) +
  theme_bw() +
  labs(title = "VO2 of mice (8wks vs control) during first visit",
       subtitle = "(`vo2.max.test.vo2_max_1`)") +
  theme(legend.position = "none")

ggsave(path = here("graphs"), filename = "pheno_vo2_visit1_vo2max.png",
       width = 6, height = 4, dpi = 300, units = "in")




# could investigate: 
# weight change over time
# post-training lactate
# body mass and body composition of 8wk vs control
# VO2 max and max run speed (MRS)


# Explore features ----
?FEATURE_TO_GENE                # 4044034 rows and 9 variables
data(FEATURE_TO_GENE)           # this takes a while, consider running the one below

feature_to_gene <- head(FEATURE_TO_GENE, 10000)
glimpse(feature_to_gene)



