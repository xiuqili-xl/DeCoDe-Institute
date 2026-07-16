# Goal ----
# Explore Carpentries lesson, RNA-seq analysis with Bioconductor
# Note, this lesson is in beta
# Also, we'll start with the Episode, Exploratory analysis and quality control
# This is where DESeq2 starts
# https://carpentries-incubator.github.io/bioc-rnaseq/04-exploratory-qc.html


# Load libraries ----


library(tidyverse)
library(here)



# Import data ----
# The effect of upper-respiratory infection on transcriptomic changes in the CNS
# GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870
# Reconstruction of the count matrix from GEO: `xy-DESeq2-Carpentries-GEO-reconstruction.R`

# Actual data used: https://raw.githubusercontent.com/carpentries-incubator/bioc-rnaseq/main/episodes/data/GSE96870_counts_cerebellum.csv
# Downloaded locally for access and preservation

counts <- read_csv(here("data_other", "GSE96870_counts_cerebellum.csv")) %>%
  column_to_rownames(var = "gene")

metadata <- read_csv(here("data_other", "GSE96870_coldata_cerebellum.csv")) %>%
  column_to_rownames(var = "sample")

# gene info
gene_info <- read_delim(file = here("data_other", "GSE96870_rowranges.tsv"), 
                        delim = "\t", col_types = c("ENTREZID" = "c")) %>%
  column_to_rownames(var = "gene")



# Experimental design ----
glimpse(metadata)

## check if all samples in the count matrix represented in the metadata
all(colnames(counts) == rownames(metadata))

## look at a couple of variables 
unique(metadata$tissue)          # only cerebellum
unique(metadata$mouse) 
length(unique(metadata$mouse))   # 22, so unique id for animal


## experimental groups
metadata %>%
  count(tissue, age, sex, infection, time, name = "no_mice")
## the two main variables for subsequent analysis are sex and time



# Assemble Summarized Experiment ----
library(SummarizedExperiment)

## make show samples and genes are in the same order!
all.equal(colnames(counts), rownames(metadata))        # sample
all.equal(rownames(counts), rownames(gene_info))       # gene

## create the SummarizedExperiment object
summarized_exp <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowRanges = as(gene_info, "GRanges"),
  colData = metadata
) 


## access various slots in a SummarizedExperiment
## access the counts
dim(assay(summarized_exp))
head(assay(summarized_exp, "counts"))      # if there are more than one assay, we'll need to specify

## access sample annotation
colData(summarized_exp)

## access the gene annotation
head(rowData(summarized_exp))


## create a better `label` for samples in metadata (i.e., colData)
summarized_exp$time
summarized_exp$label <- paste(summarized_exp$sex, summarized_exp$time, summarized_exp$mouse, sep = "_")
colData(summarized_exp)

colnames(summarized_exp)
colnames(summarized_exp) <- summarized_exp$label


## reorder samples based on sex and time
summarized_exp$group <- paste(summarized_exp$sex, summarized_exp$time, sep = "_")
summarized_exp$group

summarized_exp$group <- factor(summarized_exp$group, levels = c("Female_Day0","Male_Day0", 
                                                                "Female_Day4","Male_Day4",
                                                                "Female_Day8","Male_Day8"))

summarized_exp <- summarized_exp[, order(summarized_exp$group)]
colData(summarized_exp)


## factor labels to keep them in order in the plot
summarized_exp$label <- factor(summarized_exp$label, levels = summarized_exp$label)


## export SummarizedExperiment as a RDS file
saveRDS(summarized_exp, here("data_other", "GEO96870_SummarizedExperiment.rds"))



# Gene Annotations ----
# not sure if this section is really needed...
# seems like a random deviation
library("AnnotationDbi")
library("org.Mm.eg.db")

mapIds(x = org.Mm.eg.db,             # annotation package
       keys = "497097",              # the IDs that we know
       column = "SYMBOL",            # the vallue we want 
       keytype = "ENTREZID")         # the type of key used




# RNA-Seq | Exploratory Analysis ----
library(DESeq2)

## Read in RDS object 
## This is what Carpenty lesson did... but it seems too complexed...
## summarized_exp <- readRDS(here("data_other", "GEO96870_SummarizedExperiment.rds"))


## Construct DESeq object ----
## using DESeqDataSetFromMatrix, might simplify the workflow

metadata_dds <- metadata %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")),
         time = factor(time, levels = c("Day0", "Day4", "Day8")))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata_dds,
                              rowRanges = as(gene_info, "GRanges"),
                              design = ~ sex + time)
dds       # 41,786 genes 


## Filter out lowly expressed genes ----
dds <- dds[rowSums(counts(dds)) > 5, ]
dds       # 27,430 genes

rowData(dds)
colData(dds)


## Library size difference ----
dds$animalID <- paste(dds$sex, dds$time, dds$mouse, sep = "_")
dds$group <- paste(dds$sex, dds$time, sep = "_")
dds$libSize <- colSums(counts(dds))
colData(dds)

ggplot(dat = colData(dds) %>% as.data.frame(),
       mapping = aes(x = animalID, y = libSize/1e6, fill = group)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Sample", y = "Total count in millions")

## so we do see variability in sample library size
## which we need to adjust for



## Estimate size factor ----
dds <- estimateSizeFactors(dds)
dds

colData(dds)           
## there is now a new col, sizeFactor -- to account for differences in seq depth
## if we later request counts(dds, normalized=TRUE), DESeq will divide the counts by size factor on the fly

head(counts(dds))      
head(counts)
## running estimateSizeFactors() does NOT change the raw count value

head(counts(dds, normalized = TRUE)) 
## now we get normalized counts


## plot sizeFactor against libSize
ggplot(dat = colData(dds) %>% as.data.frame(),
       mapping = aes(x = libSize, y = sizeFactor, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "Library size", y = "Size factor")







# Clearn environment before exiting ----
rm(list = ls())



