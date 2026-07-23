# Overview -----
# MoTrPAC data exploration: WAT-SC
# One of the tissues with the most differentially expressed genes
# Looking at MoTrPAC Nature 2024: WAT-SC has a lot of differentially expressed genes, but not so many protein
# This is in contrast to the liver, which is the opposite

# (1) Pick two samples for analysis on Galaxy
# (2) DESeq analysis


# Load libraries ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(here)

library(ComplexHeatmap)
library(ggplotify)


# Clean environment ----
rm(list = ls())


# Import data ----
## Sample info ---- 
exp_sample_info <- read_csv(here("data_processed", "MoTrPAC", "exp_sample_info_copy.csv"),
                            col_types = c(rep("c", 5)))


## SRA metadata  ----
sra_rna_metadata <- read_csv(here("data_raw", "MoTrPAC", "SRARunSelector_Metadata_RNA-seq.csv"),
                             name_repair = "universal")

sra_rna_metadata_t <- read_csv(here("data_processed", "MoTrPAC", "exp_sample_info_transposed_copy.csv"))


## WAT-SC data ----
## TRNSCRPT_WATSC_DA                Differential analysis of RNA-seq datasets
## TRNSCRPT_WATSC_NORM_DATA         Normalized RNA-seq data
## TRNSCRPT_WATSC_RAW_COUNTS        RNA-seq raw counts

wat_sc_rna_raw <- TRNSCRPT_WATSC_RAW_COUNTS
wat_sc_rna_norm <- TRNSCRPT_WATSC_NORM_DATA          
wat_sc_rna_da <- TRNSCRPT_WATSC_DA

feature_to_gene <- FEATURE_TO_GENE
head(feature_to_gene)


# Select samples of Galaxy ----
## consider (1) representative sample, (2) match animals already used for liver or skm-gn

## for liver, we have SRR25251380 (female, control) and SRR25251378 (female, 8wk)
## for skm-gn, we have SRR25250934 (male, control) and SRR25250938 (male, 8wk)

## quick question: does the number of training regulated feature differ by sex?
training_regulated_features <- TRAINING_REGULATED_FEATURES
glimpse(training_regulated_features)

training_regulated_features %>%
  filter(assay == "TRNSCRPT", tissue == "WAT-SC") %>%
  count(sex, feature_ID) %>%
  count(sex)             # the same... does this mean the stats was calculate for entire group ~sex + group


## let's pick females for now...
## Invesigate the sample already used ----
## find viallabel based on SRR code
sra_rna_metadata %>%
  filter(Run %in% c("SRR25251380", "SRR25251378")) %>%
  select(Run, BioProject, BioSample, SRA.Study, Vial_Label, Sample.Name, tissue, treatment) 

## find pid based on viallabel
exp_sample_info %>%
  filter(viallabel %in% c("90258016803", "90252016803"))

## see if the same animals were used for wat-sc
wat_sc_samples_selected <- exp_sample_info %>%
  filter(pid %in% c(10045228, 10044337)) %>%
  filter(tissue == "WAT-SC")

## go back to sra_rna_metadata and see if there are any seq data from these animals
sra_rna_metadata %>%
  filter(Vial_Label %in% wat_sc_samples_selected$viallabel) %>%
  select(Run, BioProject, BioSample, SRA.Study, Vial_Label, Sample.Name, tissue, treatment)
## so there are seq data from the wat_sc from these animals
## SRR25251200 (viallabel: 90258017005) and SRR25251202 (viallabel: 90252017005)


## Test to see if these samples are representative ----
## probably easiest by looking at a heatmap of the rna_norm

## pull metadata from sra ----
wat_sc_rna_sample_metadata <- sra_rna_metadata %>%
  filter(tissue == "White Adipose Powder") %>%
  mutate(Vial_Label = as.character(Vial_Label)) %>%
  select(Run, viallabel = Vial_Label, sex, treatment) %>%
  arrange(Run)


## subset norm_data ----
## we'll focus on female, control vs 8w
wat_sc_rna_norm_f8w <- wat_sc_rna_norm %>%
  select(-feature, -tissue, -assay) %>%
  pivot_longer(cols = -feature_ID, names_to = "viallabel", values_to = "value") %>%
  ## join with sample metadata
  left_join(wat_sc_rna_sample_metadata, by = "viallabel") %>%
  ## filter for female, control vs 8 weeks 
  filter(sex == "female", str_detect(treatment, "8 weeks")) %>%
  select(feature_ID, Run, value) %>%
  arrange(Run) %>%
  pivot_wider(id_cols = "feature_ID", names_from = "Run", values_from = "value")
  

## find the differentially expressed features ----
wat_sc_rna_da_f8w_top100 <- wat_sc_rna_da %>%
  filter(sex == "female", comparison_group == "8w") %>%
  filter(selection_fdr < 0.05) %>%             
  # at this point we have 1742 rows, matching the reported no
  # we probably want to only focus on the top DE genes (top 100)
  arrange(adj_p_value) %>%
  head(100)


## prep for graphing heatmap
wat_sc_rna_norm_f8w_heatmap_data <- wat_sc_rna_norm_f8w %>%
  filter(feature_ID %in% wat_sc_rna_da_f8w_top100$feature_ID) %>%
  left_join(feature_to_gene %>% select(feature_ID, gene_symbol), by = "feature_ID") %>%
  select(-feature_ID) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  column_to_rownames(var = "gene_symbol")

wat_sc_rna_norm_f8w_heatmap_scaled <- t(scale(t(wat_sc_rna_norm_f8w_heatmap_data))) %>%
  as.data.frame()

wat_sc_rna_norm_f8w_heatmap_metadata <- wat_sc_rna_sample_metadata %>%
  filter(sex == "female", str_detect(treatment, "8 weeks")) %>%
  separate(col = treatment, into = c("treatment", NA), sep = " - ") %>%
  mutate(treatment = factor(treatment, levels = c("Control", "Training"))) %>%
  column_to_rownames(var = "Run")

all(colnames(wat_sc_rna_norm_f8w_heatmap_data) == rownames(wat_sc_rna_norm_f8w_heatmap_metadata))
# sanity check that the colnames and the rownames match


## graph heatmap! 
Heatmap(
  wat_sc_rna_norm_f8w_heatmap_scaled,
  name = "z-score",
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(
    group = wat_sc_rna_norm_f8w_heatmap_metadata$treatment,
    col = list(
      group = c(Control = "#56b4e9", Training  = "#0072b2")
    )
  )
) %>%
  as.ggplot() +
  labs(title = "Heatmap showing z-score of differentially expressed genes",
       subtitle = "(MoTrPAC female, control vs 8wk, wat-sc, top 100 DE features)")

ggsave(path = here("graphs", "MoTrPAC_DE"), 
       filename = "wat_sc_heatmap_8wk_female.png",
       width = 6, height = 8, dpi = 300, unit = "in", bg = "white")


# Hmmm.... it might be better to also have other animals
# control: SRR25251198, (viallabel: 90265017005) / SRR25251202 (viallabel: 90252017005)
# training: SRR25251199 (viallabel: 90259017005) / SRR25251200 (viallabel: 90258017005)
sra_rna_metadata %>%
  filter(Run %in% c("SRR25251198", "SRR25251199")) %>%
  select(Run, BioProject, BioSample, SRA.Study, Vial_Label, Sample.Name, tissue, treatment) 







