# Goal -----
# MoTrPAC data exploration -- focusing on liver data, control vs 8wk, male and female
# Start with TRNSCRPT_LIVER, potentially moving to PROT_LIVER and METAB_LIVER


# Load libraries ----
library(MotrpacRatTraining6moData)
library(tidyverse)
library(here)


# Import data ----
## Clean environment
rm(list = ls())

## Previously saved animal and sample data 
exp_group_assignment <- read_csv(here("data_processed", "MoTrPAC", "exp_group_assignment_copy.csv"),
                                 col_types = c(rep("c", 6), "n"))
exp_sample_info <- read_csv(here("data_processed", "MoTrPAC", "exp_sample_info_copy.csv"),
                            col_types = c(rep("c", 5)))


## SRA metadata ----
## downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA908279&o=acc_s%3Aa
sra_rna_metadata <- read_csv(here("data_raw", "MoTrPAC", "SRARunSelector_Metadata_RNA-seq.csv"),
                             name_repair = "universal")
glimpse(sra_rna_metadata)


## GEO metadata ----
geo_rna_metadata <- read_csv(here("data_raw", "MoTrPAC", "GEO_Samples_Metadata_RNA-Seq.csv"),
                             name_repair = "universal")
glimpse(geo_rna_metadata)


## Liver RNA Data ----
## to look at all available objects
data(package = "MotrpacRatTraining6moData")
## TRNSCRPT_LIVER_DA                Differential analysis of RNA-seq datasets
## TRNSCRPT_LIVER_NORM_DATA         Normalized RNA-seq data
## TRNSCRPT_LIVER_RAW_COUNTS        RNA-seq raw counts

liver_rna_raw <- TRNSCRPT_LIVER_RAW_COUNTS
liver_rna_norm <- TRNSCRPT_LIVER_NORM_DATA          
liver_rna_da <- TRNSCRPT_LIVER_DA



# What are the sample in liver RNA dataset? ----
## for examining FASTQ and maybe BAM file (i.e., repeating parts of the pre-processing)
glimpse(liver_rna_raw)

liver_rna_viallabel <- colnames(liver_rna_raw)[5:ncol(liver_rna_raw)]

liver_rna_sample_info <- exp_sample_info %>%
  filter(viallabel %in% liver_rna_viallabel) %>%
  filter(group %in% c("8w", "control"))

liver_rna_sample_info
## 5 samples per group


## Find corresponding SRA records ----
## viallabel is included in the file name of SRA record
sra_rna_metadata_liver <- sra_rna_metadata %>%
  filter(str_detect(Sample.Name, "liver")) %>%
  filter(treatment %in% c("Reference", "Control - 8 weeks", "Training - 8 weeks")) %>%
  select(Run, sex, strain, tissue, treatment, Vial_Label, starts_with("Library"), 
         Bytes, BioProject, BioSample, Experiment, SRA.Study, Sample.Name) %>%
  arrange(sex, treatment, Vial_Label) 

write_csv(sra_rna_metadata_liver, 
          here("outputs_csv", "MoTrPAC_exp", "SRARunSelector_metadata_liver.csv"))


## pick vial 3 for female, control vs 8 wk --- analysis in Galaxy
## SRR25251380 (female, control, 90252016803) and SRR25251378 (female, 8 wk, 90258016803)
sra_rna_metadata_liver %>%
  filter(Run %in% c("SRR25251380", "SRR25251378")) %>%
  view()


## Find corresponding GEO records ----
geo_rna_metadata %>%
  filter(SRA.Accession %in% c("SRX20997600", "SRX20997602")) %>%
  view()
## GEO Accession are GSM7760505 and GSM7760507
## points to SRA for data


# Try DESeq2 analysis ----
library(DESeq2)

## DESeq2 is the analysis used in the paper. From MoTrPAC Nature 2024 Suppl:
## "filtered raw counts were used as input for differential analysis with DESeq2"
glimpse(liver_rna_raw)        ## raw counts, presumably what goes into DESeq2
glimpse(liver_rna_norm)       ## normalized count, output from DESeq2?

glimpse(liver_rna_da)         ## includes control vs (1w, 2w, 4w, 8w) for male and female

## investigate the number of rows that below to each comparison in the liver_rna_da df
liver_rna_da %>%
  dplyr::count(sex, comparison_group)     
## matches the numer of rows in liver_rna_norm

## also note, nrows in liver_rna_norm is just little under 50% of liver_rna_raw
## suggest some sort of pre-filtering


## Female rats control vs 8 wk ----
## here is another decision: let's limit first attempt to just this set
liver_rna_da_f8w <- liver_rna_da %>%
  filter(sex == "female", comparison_group == "8w") %>%
  select(feature_ID, p_value, adj_p_value, starts_with("logFC"), starts_with("shrunk"),
         covariates)

head(liver_rna_da_f8w)

unique(liver_rna_da_f8w$covariates)
## "pct_globin,rin,pct_umi_dup,median_5_3_bias"
## wonder what this means... time-series?


### Set up by reshaping df ----
### create a new identifier for DESeq2 (bc it'll be used as names)
liver_rna_sample_f8w <- liver_rna_sample_info %>%
  filter(sex == "female") %>%
  mutate(PID = paste0("R", pid))       # add R for rat ;D


### create count df for DESeq
count_data_f8w <- liver_rna_raw %>%
  select(-feature, -tissue, -assay) %>%
  pivot_longer(cols = -feature_ID, names_to = "viallabel", values_to = "value") %>%
  filter(viallabel %in% liver_rna_sample_f8w$viallabel) %>%
  left_join(liver_rna_sample_f8w %>% select(viallabel, PID), by = "viallabel") %>%
  select(-viallabel) %>%
  pivot_wider(id_cols = feature_ID, names_from = PID, values_from = value) %>%
  column_to_rownames(var = "feature_ID")
  
### create metadata for DESeq 
metadata_f8w <- liver_rna_sample_f8w %>%
  select(PID, group) %>%
  column_to_rownames(var = "PID") %>%
  mutate(group = factor(group, levels = c("control", "8w")))

### rearrange count_data_f8w columns by the order of rows in metadata_f8w
count_data_f8w <- count_data_f8w[ , rownames(metadata_f8w)]

### check the col of count_data is in the same order as rows of metadata
colnames(count_data_f8w) == rownames(metadata_f8w)


### Construct DESeqDataSet ----
dds_f8w <- DESeqDataSetFromMatrix(countData = count_data_f8w,
                                  colData = metadata_f8w,
                                  design = ~ group)
dds_f8w

### pre-filtering to remove rows with low reads
summary(rowSums(counts(dds_f8w)))       # median is 8

dds_f8w <- dds_f8w[rowSums(counts(dds_f8w)) > 10, ]
dds_f8w          # has 16,106 features, could increase the threshold to cut more 


### Differential Expression ----
dds_f8w <- DESeq(dds_f8w)

res_f8w <- results(dds_f8w, contrast = c("group", "control", "8w"), alpha = 0.1)
res_f8w
summary(res_f8w) 
### 43+67 genes with adjusted p-value < 0.05

### pull these out and cross reference these genes with liver_rna_da_f8w
res_f8w_df <- as.data.frame(res_f8w)

res_f8w_sig <- res_f8w_df %>%
  filter(padj < 0.1) %>%
  rownames_to_column(var = "feature_ID")

liver_rna_da_f8w %>%
  filter(adj_p_value < 0.1) %>%      # there are only 36
  select(feature_ID, adj_p_value, starts_with("logFC")) %>%    # for now, just keep the col we'll look at
  full_join(res_f8w_sig %>% select(feature_ID, log2FoldChange, lfcSE, padj), by = "feature_ID")
## 13 features in liver_rna_da_f8w, but not in res_f8w_sig, these tend to have super large logFC
## logFC for features in both datasets are pretty close, except with opposite signs - likely issue with setting control


### MA-plot ----
plotMA(res_f8w)
plotMA(res_f8w, ylim = c(-6, 6), size = 1) 


### Plot dispersion ----
plotDispEsts(dds_f8w)


### PCA ----
vsd_8fw <- vst(dds_f8w, blind = TRUE)           # Variance-stabilizing transformation
plotPCA(vsd_8fw, intgroup = "group") 

rld_f8w <- rlog(dds_f8w, blind = TRUE)          # Regularized log transformation
plotPCA(rld_f8w, intgroup = "group")            # could use `ntop = 200` to specify the # of genes for analysis
## these look like inversions of each other...
## not showing create clustering...

