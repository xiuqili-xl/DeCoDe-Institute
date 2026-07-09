# Overview -----
# DeCoDe Institute post session exploration
# Following instructions for Omics Data and subsequent sections from `06_R_Intro.Rmd`


# Load libraries ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(here)
library(ggfortify)
library(ggrepel)


# Load data ----
rm(list = ls())

## to look at all available objects
data(package = "MotrpacRatTraining6moData")

## subcutaneous white adipose RNA data sets include
## TRNSCRPT_WATSC_DA                Differential analysis of RNA-seq datasets
## TRNSCRPT_WATSC_NORM_DATA         Normalized RNA-seq data
## TRNSCRPT_WATSC_RAW_COUNTS        RNA-seq raw counts

watsc_rna <- TRNSCRPT_WATSC_NORM_DATA          
watsc_rna_da <- TRNSCRPT_WATSC_DA
watsc_rna_raw <- TRNSCRPT_WATSC_RAW_COUNTS

feature_to_gene <- FEATURE_TO_GENE
pheno <- PHENO


# Explore NORM_DATA----
?TRNSCRPT_WATSC_NORM_DATA
## Normalized sample-level RNA-seq (TRNSCRPT) data used for visualization
## A data frame with genes in rows (feature_ID) and samples in columns (viallabel)
## Filtering of lowly expressed genes and normalization were performed separately in each tissue.
## Aka, lowly expressed genes has already been removed from NORM_DATA
glimpse(watsc_rna)

?FEATURE_TO_GENE
glimpse(feature_to_gene)
# note the last three columns are not empty, though most of them are
unique(feature_to_gene$relationship_to_gene)
unique(feature_to_gene$custom_annotation)
unique(feature_to_gene$kegg_id)



# Combine with gene info  ----
# inner join, retain only rows in both sets
watsc_rna_gene <- inner_join(watsc_rna, feature_to_gene, by = "feature_ID")

glimpse(watsc_rna_gene)
# note, there are more rows in watsc_rna_gene than in watsc_rna

watsc_rna %>%
  count(feature_ID) %>%
  filter(n >= 2)
# feature_ID in watsc_rna are unique

feature_to_gene %>%
  count(feature_ID) %>%
  filter(n >= 2)
# feature_ID in watsc_rna are not unique...

# question: what feature_ID are now repeated
watsc_rna_gene_feature_ID_count <- watsc_rna_gene %>%
  count(feature_ID) %>%
  filter(n >= 2)

feature_to_gene %>%
  filter(feature_ID %in% watsc_rna_gene_feature_ID_count$feature_ID) %>%
  arrange(feature_ID) %>%
  view()


# Rmd suggests that we limit our analysis to feature_id's that correspond to 1 gene
watsc_rna_gene_filtered <- watsc_rna_gene %>%
  group_by(gene_symbol) %>%
  filter(n() == 1) %>%
  ungroup()

glimpse(watsc_rna_gene_filtered)
# 15,963 rows; 700 fewer than watsc_rna



# Combine with sample, animal info ----
## first, pull sample info from pheno
glimpse(pheno) 

sample_info <- pheno %>%
  distinct(viallabel, pid, group, sex, tissue) %>%
  mutate(pid = as.character(pid)) %>%
  arrange(sex, group, pid, tissue, viallabel) 

glimpse(sample_info)
write_csv(sample_info, here("outputs_csv", "exp_sample_info.csv"))


# wrangle watsc_rna_gene_filtered 
watsc_rna_wide <- watsc_rna_gene_filtered %>%
  select(-c(feature, entrez_gene:rgd_gene, old_gene_symbol:kegg_id)) %>%
  pivot_longer(cols = -c(feature_ID, gene_symbol, tissue, assay), 
               names_to = "viallabel", values_to = "expression") %>%
  pivot_wider(id_cols = c("viallabel", "tissue", "assay"), 
              names_from = gene_symbol, values_from = expression)

## each row has a unique viallabel
length(unique(watsc_rna_wide$viallabel)) == nrow(watsc_rna_wide)

# combine with sample_info
watsc_exp <- watsc_rna_wide %>%
  left_join(sample_info, by = join_by(viallabel, tissue)) %>% 
  select(viallabel, pid, group, sex, tissue, assay, everything())


# PCA graph ----
watsc_exp_pca <- prcomp(select(watsc_exp, -c(viallabel:assay)))     
# prcomp only wants the numeric expression matrix

watsc_exp_pca

autoplot(watsc_exp_pca, data = watsc_exp, color = "group", shape = "sex") 
# by default graphs the first two principal components

autoplot(watsc_exp_pca, data = watsc_exp, x = 2, y = 3, color="group", shape="sex")
# graphs the 2nd and 3rd principal component

# Question: why do the feature_id & gene inner_join + filtering when we could do the PCA on feature ID?


## PCA of male rats 8wk vs control
watsc_exp_subset <- watsc_exp %>%
  filter(group %in% c("8w", "control"))

watsc_exp_subset_pca <- prcomp(select(watsc_exp_subset, -c(viallabel:assay)))    

autoplot(watsc_exp_subset_pca, data = watsc_exp_subset, color = "group", shape = "sex", size = 3) +
  scale_x_continuous(limits = c(-0.3, 0.3), expand = 0.005) +
  scale_y_continuous(limits = c(-0.3, 0.3), expand = 0.005) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "PCA of Subcutaneous WAT gene expression",
       subtitle = "(8wk endurance training vs control; female and male)")

ggsave(path = here("graphs"), filename = "watsc_pac_8wk.png",
       width = 6, height = 5, dpi = 300, units = "in")


# Volcano plot ----
glimpse(watsc_rna_da)

unique(watsc_rna_da$comparison_group)

watsc_rna_da %>%
  count(feature_ID, sex) %>%
  filter(n != 4)
# looks like all comparison are made to the control group, 
# thus, for each feature there are 4 rows for male and 4 rows for female

watsc_rna_da_graph_male <- watsc_rna_da %>%
  select(feature_ID, sex, comparison_group, logFC, adj_p_value) %>%
  filter(sex == "male", comparison_group == "8w") %>%
  inner_join(feature_to_gene, by = "feature_ID") %>%
  group_by(gene_symbol) %>%
  filter(n() == 1) %>%
  mutate(point_color = if_else(adj_p_value < 0.05, "coral", "grey")) 

ggplot() +
  geom_point(data = watsc_rna_da_graph_male,
             mapping = aes(x = logFC, y = -log10(adj_p_value), fill = point_color),
             shape = 21) +
  geom_text_repel(data = watsc_rna_da_graph_male %>% filter(-log10(adj_p_value) > 5),
             mapping = aes(x = logFC, y = -log10(adj_p_value), label = gene_symbol),
             size = 2, box.padding = 0.3) +
  scale_x_continuous(limits = c(-60, 60), expand = 0.01) +
  scale_y_continuous(limits = c(0, 28), expand = 0.03) +
  scale_fill_manual(values = c("coral" = "coral", "grey" = "grey")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "PCA of Subcutaneous WAT gene expression",
       subtitle = "(8wk endurance training vs control; male)")

ggsave(path = here("graphs"), filename = "watsc_volcano_8wk_male.png",
       width = 6, height = 5, dpi = 300, units = "in")


