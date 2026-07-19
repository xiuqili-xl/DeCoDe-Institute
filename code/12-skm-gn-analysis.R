# Overview -----
# MoTrPAC data exploration -- focusing on skeletal muscle (gastrocnemius)
# the first STAR alignment on Galaxy to finish is for sample SRR25250934
# which is skm-gn_control_male_8w_biological_rep5_vial_90239015512


# Load libraries ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(here)

library(ComplexHeatmap)
library(ggplotify)


# Clean environment ----
rm(list = ls())


# Look at mete-data ----
## SRA metadata  ----
sra_rna_metadata <- read_csv(here("data_raw", "MoTrPAC", "SRARunSelector_Metadata_RNA-seq.csv"),
                             name_repair = "universal")
sra_rna_metadata_t <- sra_rna_metadata %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = -Run, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(id_cols = Variable, names_from = Run, values_from = Value)

#write_csv(sra_rna_metadata_t,
#          here("outputs_csv", "MoTrPAC_exp", "exp_sample_info_transposed.csv"))

sra_rna_metadata_t %>%
  select(Variable, SRR25250934) %>%
  print(n = 35)
## control male, rep5, viallabel == 90239015512 


## Sample info ----
exp_sample_info <- read_csv(here("data_processed", "MoTrPAC", "exp_sample_info_copy.csv"),
                            col_types = c(rep("c", 5)))
exp_sample_info %>%
  filter(viallabel == "90239015512" )
## pid == 10027599




# Import SKM-GN TRNSCRPT data  ----
## TRNSCRPT_SKMGN_DA                Differential analysis of RNA-seq datasets
## TRNSCRPT_SKMGN_NORM_DATA         Normalized RNA-seq data
## TRNSCRPT_SKMGN_RAW_COUNTS        RNA-seq raw counts

skm_gn_rna_raw <- TRNSCRPT_SKMGN_RAW_COUNTS
skm_gn_rna_norm <- TRNSCRPT_SKMGN_NORM_DATA          
skm_gn_rna_da <- TRNSCRPT_SKMGN_DA

feature_to_gene <- FEATURE_TO_GENE
head(feature_to_gene)


# Quick analysis | DE genes in control vs 8 wk males ----
## use skm_gn_rna_da directly
skm_gn_rna_da_m8w <- skm_gn_rna_da %>%
  filter(sex == "male", comparison_group == "8w") %>% 
  select(feature_ID, ends_with("intensity"), ends_with("logFC"), ends_with("p_value")) %>%
  left_join(feature_to_gene %>% select(feature_ID, gene_symbol), by = "feature_ID") %>%
  arrange(adj_p_value) %>%
  select(feature_ID, gene_symbol, everything())

head(skm_gn_rna_da_8w) 

## visualize using a volcano plot
ggplot(data = skm_gn_rna_da_m8w,
       mapping = aes(x = shrunk_logFC,
                     y = -log10(adj_p_value))) +
  geom_point(shape = 21, color = "grey25", fill = "grey", alpha = 60) +
  scale_x_continuous(limits = c(-3.1, 3.1), breaks = seq(-3, 3, by = 1),
                     expand = 0.005, oob = scales::squish) +
  scale_y_continuous(limits = c(0, 15), expand = 0.02) +
  theme_bw()


## explore genes with abs(shrunk_logFC) > 1 -- doubled or halfed
skm_gn_rna_da_m8w %>%
  filter(abs(shrunk_logFC) >=1) %>%
  select(feature_ID, gene_symbol, ref_int = reference_average_intensity,
         comp_int = comparison_average_intensity, everything())
## 12 genes with Erfe! -- nice to see an old friend

feature_to_gene %>%
  filter(feature_ID == "ENSRNOG00000024688")


## explore genes with adj_p_value < 0.05
skm_gn_rna_da_m8w_sig <- skm_gn_rna_da_m8w %>%
  filter(adj_p_value < 0.05)

head(skm_gn_rna_da_m8w_sig)
nrow(skm_gn_rna_da_m8w_sig)      # 102 features, to too bad ;D


## create a heatmap using skm_gn_rna_norm, so we can visualize representative samples
## ultimate goal: find a SRA to pull from (skm-ng, 8w training, male) for Galaxy
unique(skm_gn_rna_norm$tissue)
unique(skm_gn_rna_norm$assay)

### pull relevant meta data for skm_gn_rna (include both sex and all weeks)
skm_gn_rna_metadata <- sra_rna_metadata %>%
  select(sex, tissue, treatment, viallabel = Vial_Label, Sample.Name) %>%
  filter(str_detect(Sample.Name, "skm-gn"), treatment != "Reference") %>%
  arrange(treatment) %>%
  mutate(rep = str_extract(Sample.Name , "rep\\d+"),
         rep = str_replace(rep, "rep", "Rep"),
         viallabel = as.character(viallabel)) %>%
  select(viallabel, sex, tissue, treatment, rep) %>%
  arrange(treatment, rep)


### reorganize rna_norm data (restrict to control and 8wk, male)
skm_gn_rna_norm_m8w <-  skm_gn_rna_norm %>%
  select(-feature, -tissue, -assay) %>%
  pivot_longer(cols = -feature_ID, names_to = "viallabel", values_to = "norm_exp") %>%
  ## join with metadata and clean up samples (selection, name, etc)
  left_join(skm_gn_rna_metadata, by = "viallabel") %>%
  filter(sex == "male", str_detect(treatment, "8 weeks")) %>%
  mutate(treatment = str_remove(treatment, " - 8 weeks"),
         sample = paste(treatment, rep, sep = "_")) %>%
  select(feature_ID, sample, norm_exp) %>%
  ## pivot wider 
  pivot_wider(id_cols = feature_ID, names_from = sample, values_from = norm_exp) %>%
  ## join with feature_to_gene to get gene names
  left_join(feature_to_gene %>% select(feature_ID, gene_symbol), by = "feature_ID") %>%
  select(feature_ID, gene_symbol, everything())
  

## subset rna_norm data to just the DE genes
skm_gn_rna_norm_m8w_heatmap_data <- skm_gn_rna_norm_m8w %>%
  filter(feature_ID %in% skm_gn_rna_da_m8w_sig$feature_ID) %>%
  select(-feature_ID) %>%
  column_to_rownames(var = "gene_symbol") %>%
  relocate(sort(names(.))) 

skm_gn_rna_norm_m8w_heatmap_scaled <- t(scale(t(skm_gn_rna_norm_m8w_heatmap_data))) %>%
  as.data.frame()


### create metadata for heatmap annotation
skm_gn_rna_norm_m8w_heatmap_metadata <- data.frame(
  sample = colnames(skm_gn_rna_norm_m8w_heatmap_data)
) %>% 
  separate(col = sample, into = c("group", "rep"), remove = FALSE) %>%
  column_to_rownames(var = "sample")


### generate heatmp
Heatmap(
  skm_gn_rna_norm_m8w_heatmap_scaled,
  name = "z-score",
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(
    group = skm_gn_rna_norm_m8w_heatmap_metadata$group,
    col = list(
      group = c(Control = "#56b4e9", Training  = "#0072b2")
    )
  )
) %>%
  as.ggplot() +
  labs(title = "Heatmap showing z-score of differentially expressed genes",
       subtitle = "(MoTrPAC male, control vs 8wk, skm-gn, adj pvalue < 0.05)")

ggsave(path = here("graphs", "MoTrPAC_DE"), 
       filename = "skm_gn_heatmap_8wk_male.png",
       width = 6, height = 8, dpi = 300, unit = "in", bg = "white")


### look up Erfe
skm_gn_rna_norm_m8w_heatmap_scaled ["Erfe", ]    # down regulated in Training

### let's pick Training_Rep5 for Galaxy investigation
sra_rna_metadata %>%
  filter(sex == "male", treatment == "Training - 8 weeks",
         str_detect(Sample.Name, "skm-gn")) %>%
  select(Run, Assay.Type, Sample.Name, treatment, Vial_Label)
### Run == SRR25250938; Vial_Label == 90227015512

sra_rna_metadata %>%
  filter(sex == "male", treatment == "Control - 8 weeks",
         str_detect(Sample.Name, "skm-gn")) %>%
  select(Run, Assay.Type, Sample.Name, treatment, Vial_Label)
