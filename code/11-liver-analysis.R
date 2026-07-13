# Overview -----
# MoTrPAC data exploration -- focusing on liver data, control vs 8wk, male and female
# Start with TRNSCRPT_LIVER, potentially moving to PROT_LIVER and METAB_LIVER


# Load libraries ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(here)


# Import data ----
## Previously saved animal and sample data 
exp_group_assignment <- read_csv(here("data_processed", "exp_group_assignment_copy.csv"),
                                 col_types = c(rep("c", 6), "n"))
exp_sample_info <- read_csv(here("data_processed", "exp_sample_info_copy.csv"),
                            col_types = c(rep("c", 5)))


## SRA metadata 
## downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA908279&o=acc_s%3Aa
sra_rna_metadata <- read_csv(here("data_raw", "SRARunSelector_Metadata_RNA-seq.csv"),
                             name_repair = "universal")
glimpse(sra_rna_metadata)


## GEO metadata 
geo_rna_metadata <- read_csv(here("data_raw", "GEO_Samples_Metadata_RNA-Seq.csv"),
                             name_repair = "universal")
glimpse(geo_rna_metadata)



# TRNSCRPT_LIVER ----
## Import data ----
rm(list = ls())

## to look at all available objects
data(package = "MotrpacRatTraining6moData")
## TRNSCRPT_LIVER_DA                Differential analysis of RNA-seq datasets
## TRNSCRPT_LIVER_NORM_DATA         Normalized RNA-seq data
## TRNSCRPT_LIVER_RAW_COUNTS        RNA-seq raw counts

liver_rna_raw <- TRNSCRPT_LIVER_RAW_COUNTS
liver_rna_norm <- TRNSCRPT_LIVER_NORM_DATA          
liver_rna_da <- TRNSCRPT_LIVER_DA



## What are the sample in liver RNA dataset? ----
## for examining FASTQ and maybe BAM file (i.e., repeaing earlier parts of the preprocessing)

glimpse(liver_rna_raw)

liver_rna_viallabel <- colnames(liver_rna_raw)[5:ncol(liver_rna_raw)]

liver_rna_sample_info <- exp_sample_info %>%
  filter(viallabel %in% liver_rna_viallabel) %>%
  filter(group %in% c("8w", "control"))

liver_rna_sample_info
## 5 samples per group


## Find corresponding SRA records ----
## viallabel is included in the file name of SRA record
sra_rna_metadata_for_analysis <- sra_rna_metadata %>%
  filter(str_detect(Sample.Name, "liver")) %>%
  filter(treatment %in% c("Reference", "Control - 8 weeks", "Training - 8 weeks")) %>%
  select(Run, sex, strain, tissue, treatment, Vial_Label, starts_with("Library"), 
         Bytes, BioProject, BioSample, Experiment, SRA.Study, Sample.Name) %>%
  arrange(sex, treatment, Vial_Label) 

write_csv(sra_rna_metadata_for_analysis, here("outputs_csv", "SRARunSelector_metadata_for_analysis.csv"))

## pick vial 3 for female, control vs 8 wk --- analysis in Galaxy
## SRR25251380 (female, control, ) and SRR25251378 (female, 8 wk)
sra_rna_metadata_for_analysis %>%
  filter(Run %in% c("SRR25251380", "SRR25251378")) %>%
  view()


## Find corresponding GEO records ----
geo_rna_metadata %>%
  filter(SRA.Accession %in% c("SRX20997600", "SRX20997602")) %>%
  view()
## GEO Accession are GSM7760505 and GSM7760507
## points to SRA for data


