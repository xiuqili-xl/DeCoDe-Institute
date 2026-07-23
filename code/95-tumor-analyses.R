# Goal ----
# Follow tutorial from UT Austin https://nci-iteb.github.io/tumor_epidemiology_approaches/sessions/session_11/practical
# Should serve as a review of DESeq (starting with raw count data)
# Expand to gene set and pathway enrichment analysis


# Load libraries -----
library(tidyverse)
library(here)

library(DESeq2)


# Import data ----
htseq_output <- read_delim(here("data_raw", "TumorAnalyses", "Practical_11_htseq_counts.txt")) %>% 
  data.frame() %>%
  dplyr::rename(gene_name = UID) 

metadata <- read_delim(here("data_raw", "TumorAnalyses", "Practical_11_samples_clinical.txt"),
                       name_repair = "universal")

## inspect data
head(htseq_output)

head(metadata)
unique(metadata$Patient.ID)       # 10 patients, one tumor sample and one normal sample from each



# Normalization ----
## note, within provide script, there is script for TMP normalization

## we'll proceed with DESeq normalization
## "known as the median of ratios. This is a method where counts divided by sample-specific 
## size factors are determined by the median ratio of gene counts relative to geometric mean 
## per gene, taking into account sequencing depth and RNA composition."


## Clean up data and metadata ----
count_data <- htseq_output %>%
  column_to_rownames(var = "gene_name")

metadata_cleaned <- metadata %>%
  mutate(type = factor(Type, levels = c("Normal", "Tumor")),
         stage = case_when(str_detect(STAGE,'^I[AB]') ~ 1,
                           str_detect(STAGE,'^II[AB]') ~ 2,
                           str_detect(STAGE,'^III[AB]') ~ 3,
                           TRUE ~ 888),
         stage = factor(stage, levels = c("1", "2", "3"))) %>%
  select(UID, patient_id = Patient.ID, type, stage) %>%
  column_to_rownames(var = "UID")
  
## check that the rownames and colnames match
all(colnames(count_data) == rownames(metadata_cleaned))


## Create the DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata_cleaned,
                              design= ~ stage + type)
dds                     # dim: 59609 20 


## Filter out genes with low counts ----
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds                     # dim: 47835 20 


## Normalize! -----
dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name")

head(normalized_counts)        
## and we can export this, bc certain analysis will use normalize count
## but not DESeq, which uses the raw counts



# PCA ----
## Run a vst normalization ----
## vst = Variance stabilizing transformation
## will use due to its speed and effectiveness
## https://alexslemonade.github.io/refinebio-examples/03-rnaseq/00-intro-to-rnaseq.html#deseq2-transformation-methods

vsd <- vst(dds)

plotPCA(vsd, intgroup = c("type", "stage"))
plotPCA(vsd, intgroup = c("type"))
plotPCA(vsd, intgroup = c("stage"))

PCA_data <- plotPCA(vsd, intgroup = c("type", "stage"), returnData = TRUE)

ggplot(data = PCA_data, 
       mapping = aes(x = PC1, y = PC2, shape = type, color = stage)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_x_continuous(limits = c(-75, 60), expand = 0.01) +
  scale_y_continuous(limits = c(-75, 60), expand = 0.01) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "Tumor Analyses | PCA",
       x = "PC1 (63% variance)", y = "PC2 (6% variance)")

ggsave(path = here("graphs", "TumorAnalyses"), filename = "TumorAnalyses_PCA.png",
       width = 6, height = 5, dpi = 300, unit = "in")



# Differential Expression Analysis ----
## IMPORTANT! DESeq() require raw counts as input, not normalized/transformed count 

## Run DESeq() ----
dds <- DESeq(dds)
dds

assay(dds, "counts") %>% head()        # still the raw counts!
count_data %>% head()


## Pull results ----
## set threshold: FDR of 0.05 and logFC of at least 1
res <- results(dds, alpha = 0.05, lfcThreshold = 1)
res
summary(res)


## look at some of the genes that are most upregulated or downregulated
as.data.frame(res) %>%
  # filter for genes with adj pvalue < 0.05
  filter(padj < 0.05) %>%
  # filter for genes whose regulation was downregulated by 1/2
  filter(log2FoldChange < -1) %>%         # 1,425 genes at this point, matching summary(res)
  # arrange df so genes that are most suppressed show up on top
  arrange(log2FoldChange) %>%
  head(10)



# Visualization of DE results -----

## Volcano plot ----
## use lfcShrink to shrink the log2 fold change
res_shrink <- lfcShrink(dds,
                        contrast = c("type", "Tumor", "Normal"),
                        res = res, type = "normal")
res_shrink

library(EnhancedVolcano)
EnhancedVolcano(res_shrink,
                lab = rownames(res_shrink),
                pointSize = 1.5,
                labSize = 3,
                title = "Tumor Analyses | Volcano Plot",
                subtitle = NULL,
                x = "log2FoldChange",
                y = "pvalue")

ggsave(path = here("graphs", "TumorAnalyses"), filename = "TumorAnalyses_Volcano.png",
       width = 6.5, height = 7, dpi = 300, unit = "in")

## EnhancedVolcano is a neat package, but a little slow. probably prefer good old ggplot...



## DESeq2::plotCounts() ----
## plot normalized count for a single gene
plotCounts(dds, gene = "A2M", intgroup = "type")
normalized_counts %>% filter(gene_name == "A2M")

plotCounts(dds, gene = "TP53", intgroup = "type")
as.data.frame(res)["TP53", ]
## this is pretty neat, could be a quick viz after running estimateSizeFactors()
  

## Heatmap ----
## find top DE genes
de_genes <- as.data.frame(res) %>%
  # filter for genes with adj pvalue < 0.05, abs(logFC) > 1
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  .[1:1000, ] %>%
  rownames_to_column(var = "gene_name") %>%
  dplyr::pull(gene_name)

heatmap_data <- normalized_counts %>%
  filter(gene_name %in% de_genes) %>%
  pivot_longer(col = -gene_name, names_to = "sample", values_to = "exp") %>%
  left_join(metadata_cleaned %>% rownames_to_column(var = "sample"), by = "sample") %>%
  # calculate z-score for each gene
  group_by(gene_name) %>%
  mutate(exp = (exp - mean(exp)) / sd(exp)) %>%
  ungroup()


library(tidyHeatmap)
tidyHeatmap::heatmap(.data = heatmap_data,
                     .row = gene_name,
                     .column = sample,
                     .value = exp,
                     clustering_method_columns = "ward.D2",
                     clustering_method_rows = "ward.D2",
                     show_row_names = FALSE,
                     column_title = "Tumor Analyses | Heatmap") %>%
  annotation_tile (type) %>%
  annotation_tile (stage, palette = c("#1b9e77", "#d95f02", "#7570b3")) %>%
  save_pdf(here("graphs", "TumorAnalyses", "TumorAnalyses_Heatmap.pdf"),
           width = 7, height = 7)

ggsave(path = here("graphs", "TumorAnalyses"), filename = "TumorAnalyses_Heatmap.png",
       width = 7, height = 7, dpi = 300, unit = "in")

## I like tidyHeatmap!
## Probably don't need to graph the top 1,000 DE genes...



# Pathway analysis with GSEA ----
## https://www.genepattern.org/ no longer works

## look at the provided DESeq2_normalized_count data
norm_counts_for_GSEA <- read.delim(here("data_raw", "TumorAnalyses", "Practical_11_DESeq2_normalized_counts.gct"),
                                   skip = 2)

head(norm_counts_for_GSEA) %>% select(-UID)
head(normalized_counts)

round(norm_counts_for_GSEA[ , 3:22], 3) == round(normalized_counts[ , 2:21], 3)
all.equal(norm_counts_for_GSEA[ , 3:22], normalized_counts[ , 2:21])

