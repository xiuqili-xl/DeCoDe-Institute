# Goal ----



# Dataset ----

## gtex_tissue_train_200x100.tsv ----
# A subset of GTEx transcripts_tpm (genes and people) used during DeCoDe Galaxy Day 2 Training
# tissue with a lot of expression data + genes with bigger variability in expression (ANOVA)
# http://data.schatz-lab.org/cwic_training/
# also avaiable on https://github.com/mschatz/data/tree/main/cwic_training)

# Load libraries ----
library(tidyverse)
library(here)


# Import data ----
train <- read_delim(here("data_raw", "gtex_tissue_train_200x100.tsv"))
test <- read_delim(here("data_raw", "gtex_tissue_test_labeled_100x100.tsv"))


# Explore data ----
glimpse(train)                # columns: SampleID, Tissue, and 100 genes

unique(train$SampleID)        # 200 unique Sample ID
## named GTEX-<DonorID>-<TissueID>

unique(train$Tissue)          # 10 different types of tissues
train %>% count(Tissue)       # 20 donors per tissue
