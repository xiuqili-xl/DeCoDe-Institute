# Goal -----
# Replicate DESeq2 analysis that was performed on Galaxy during the DeCoDe Institute


# Useful resources ----
# Tufts quick guide https://rtguides.it.tufts.edu/bio/tutorials/de-seq2-in-r-ood.html
# Schwartz blog https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial
# Full vignette: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# Dataset ----
# Transcriptome of Embryonic and Adult mouse cerebral cortex
# Dataset for teaching https://zenodo.org/records/20531535
# Note FM1 was not used in the Galaxy analysis; there was something wrong with it

# Load library ----
library(tidyverse)
library(here)
library(DESeq2)
library(pheatmap)
library(ggplotify)


# Import data ----
count_AM1 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_AM1.txt", 
                        sep = "\t", col.names = c("Feature", "AM1"))
count_AM2 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_AM2.txt", 
                        sep = "\t", col.names = c("Feature", "AM2"))
count_AM3 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_AM3.txt", 
                        sep = "\t", col.names = c("Feature", "AM3"))
count_FM3 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_FM3.txt", 
                        sep = "\t", col.names = c("Feature", "FM3"))
count_FM4 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_FM4.txt", 
                        sep = "\t", col.names = c("Feature", "FM4"))
count_FM5 <- read.delim(file = "https://zenodo.org/records/20531535/files/gene_count_FM5.txt", 
                        sep = "\t", col.names = c("Feature", "FM5"))
## note, all dataset has 32,711 rows

all(count_AM1$Feature == count_FM4$Feature)
## presumably, they all have the same feature


# Combine into Count Matrix ----
count_combined <- count_AM1 %>%
  inner_join(count_AM2, by = "Feature") %>%
  inner_join(count_AM3, by = "Feature") %>%
  inner_join(count_FM3, by = "Feature") %>%
  inner_join(count_FM4, by = "Feature") %>%
  inner_join(count_FM5, by = "Feature")

nrow(count_combined)         # inner_join() confirms that all datasets have the same number of features

## remove individual samples from environment
rm(count_AM1, count_AM2, count_AM3, count_FM3, count_FM4, count_FM5)

## construct count matrix for DESEq2
count_data <- count_combined %>%
  column_to_rownames(var = "Feature")

head(count_data)
class(count_data)


# Create metadata df ----
metadata_df <- data.frame(Sample = colnames(count_data)) %>%
  mutate(DevStage = case_when(str_detect(Sample, "AM") ~ "Adult",
                              str_detect(Sample, "FM") ~ "Embryonic"),
         DevStage = factor(DevStage, levels = c("Adult", "Embryonic"))) %>%
  column_to_rownames(var = "Sample")

metadata_df


## check the col of count_data is in the same order as rows of metadata_df
colnames(count_data) == rownames(metadata_df)



# Construct DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata_df,
                              design = ~ DevStage)
dds

## We can pre-filter to remove rows in which there are very few reads, thus reducing the required 
## memory, and increasing the speed. Prefiltering can also improve visualizations, as features 
## with no information for differential expression are not plotted.
## dds <- dds[rowSums(counts(dds)) > 10, ]

## However, the default in Galaxy does not perform pre-filtering
## So we'll first proceed without filtering



# Replicating the Galaxy analysis ----
## Get normalized counts ----
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

## check against galaxy output - appear to be the same
normalized_counts["Tfrc", ]

## can write to file for downstream analysis


## Differential Expression Analysis ----
dds <- DESeq(dds)

res <- results(dds, contrast = c("DevStage", "Adult", "Embryonic"),
               alpha = 0.05)       # note, unless specified, alpha defaults to 0.1
res
summary(res)           

res_df <- as.data.frame(res)

res_df["Tfrc", ]
## spot check matches galaxy output!



## MA-plot ----
plotMA(res)
plotMA(res, ylim = c(-13, 13), size = 1) 


## Plot dispersion ----
plotDispEsts(dds)


## PCA ----
vsd <- vst(dds, blind = TRUE)           # Variance-stabilizing transformation
plotPCA(vsd, intgroup = "DevStage") 

rld <- rlog(dds, blind = TRUE)          # Regularized log transformation
plotPCA(rld, intgroup = "DevStage")     # could use `ntop = 200` to specify the # of genes for analysis
## these two looks more similar than galaxy output

PCA_data <- plotPCA(rld, intgroup = "DevStage", returnData = TRUE) 
PCA_data

ggplot(data = PCA_data, 
       mapping = aes(x = PC1, y = PC2, color = DevStage, fill = DevStage)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) + 
  scale_x_continuous(limits = c(-40, 40), expand = 0.01) +
  scale_y_continuous(limits = c(-40, 40), expand = 0.01) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "Mouse Adult vs Embryonic Cortex | PCA",
       x = "PC1 (98% variance)", y = "PC2 (1% variance)")

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_PCA.png",
       width = 6, height = 4, dpi = 300, unit = "in")



## Sample-to-Sample Distance ----
## extact transformed matrics and transpose it
## use vst for very large datasets
sampleDists <- dist(t(assay(rld)))
sampleDists

## convert to a distance matrix
sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix

## create a sample-to-sample distance heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(30)) %>%
  as.ggplot() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Sample Distance Heatmap")
## visually identical to galaxy output

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Sample-Distance-Heatmap.png",
       width = 6, height = 6, dpi = 300, unit = "in", bg = "white")



## Volcano Plot ----
head(res_df, 10)

ggplot(data = res_df, 
       mapping = aes(x = log2FoldChange, y = -log(padj, 10))) +
  geom_point(size = 1, shape = 21) +
  theme_bw() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Volcano Plot")
## wonder if filtering out lowly expressed genes earlier would make a difference?

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Volcano-plot.png",
       width = 6, height = 6, dpi = 300, unit = "in", bg = "white")



## Heatmap ----
## work with log transformed data, so variance is approximately the same across different mean values
rld_df <- assay(rld) 
head(rld_df)

## find the top 100 varying genes
topVarGenes <- head(order(rowVars(rld_df), decreasing = TRUE), 200)
topVarGeneCounts <- rld_df[topVarGenes, ]

## plot using pheatmap
pheatmap(topVarGeneCounts,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         scale = "row",                          # scale by gene
         show_rownames = FALSE,
         fontsize = 6,
         #cutree_cols = 2,
         annotation_col = metadata_df,
         annotation_colors = list(
           DevStage = c("Adult" = "#009E73", "Embryonic" = "#CC79A7")
         )) %>%
  as.ggplot() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Heatmap",
       subtitle = "(top 200 varying gene) \n")

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Heatmap.png",
       width = 5, height = 5, dpi = 300, unit = "in", bg = "white")



# If we pre-filter the data ----
## Pre-filtering ----
dds_filtered <- dds[rowSums(counts(dds)) > 10, ]
## reduced the # of genes by 1/3


## Differential Expression Analysis ----
dds_filtered <- DESeq(dds_filtered)
res_filtered <- results(dds_filtered, contrast = c("DevStage", "Adult", "Embryonic"),
                        alpha = 0.05)          # raise threshold bc so many genes are significant

res_filtered
summary(res_filtered)           
## note, the results function default to adjusted p-value < 0.1
## we could specify cutoff using results(dds, alpha = 0.05)

res_filtered_df <- as.data.frame(res_filtered)


## MA-plot ----
plotMA(res_filtered)
plotMA(res_filtered, ylim = c(-13, 13), size = 1) 


## Plot dispersion ----
plotDispEsts(dds_filtered)


## PCA ----
rld_filtered <- rlog(dds_filtered, blind = TRUE)          # Regularized log transformation
plotPCA(rld_filtered, intgroup = "DevStage")

PCA_filtered_data <- plotPCA(rld_filtered, intgroup = "DevStage", returnData = TRUE) 
PCA_filtered_data

ggplot(data = PCA_filtered_data, 
       mapping = aes(x = PC1, y = PC2, color = DevStage, fill = DevStage)) +
  geom_point(shape = 21, size = 3, alpha = 0.7) + 
  scale_x_continuous(limits = c(-40, 40), expand = 0.01) +
  scale_y_continuous(limits = c(-40, 40), expand = 0.01) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "Mouse Adult vs Embryonic Cortex | PCA (filtered)",
       x = "PC1 (98% variance)", y = "PC2 (1% variance)")
# very sutble change

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_PCA_filtered.png",
       width = 6, height = 4, dpi = 300, unit = "in")



## Sample-to-Sample Distance ----
## extact transformed matrics and transpose it
## use vst for very large datasets
sampleDists_filtered <- dist(t(assay(rld_filtered)))
sampleDists_filtered

## convert to a distance matrix
sampleDistMatrix_filtered <- as.matrix(sampleDists_filtered)
sampleDistMatrix_filtered

## create a sample-to-sample distance heatmap
pheatmap(sampleDistMatrix_filtered,
         clustering_distance_rows = sampleDists_filtered,
         clustering_distance_cols = sampleDists_filtered,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(30)) %>%
  as.ggplot() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Sample Distance Heatmap (filtered)")
## visually identical to galaxy output

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Sample-Distance-Heatmap_filtered.png",
       width = 6, height = 6, dpi = 300, unit = "in", bg = "white")



## Volcano Plot ----
head(res_filtered_df, 10)

ggplot(data = res_filtered_df, 
       mapping = aes(x = log2FoldChange, y = -log(padj, 10))) +
  geom_point(size = 1, shape = 21) +
  theme_bw() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Volcano Plot (filtered)")
## wonder if filtering out lowly expressed genes earlier would make a difference?

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Volcano-plot_filtered.png",
       width = 6, height = 6, dpi = 300, unit = "in", bg = "white")



## Heatmap ----
## work with log transformed data, so variance is approximately the same across different mean values
rld_filtered_df <- assay(rld_filtered) 
head(rld_filtered_df)

## find the top 100 varying genes
topVarGenes_filtered <- head(order(rowVars(rld_filtered_df), decreasing = TRUE), 200)
topVarGeneCounts_filtered <- rld_filtered_df[topVarGenes_filtered, ]

## plot using pheatmap
pheatmap(topVarGeneCounts_filtered,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         scale = "row",                          # scale by gene
         show_rownames = FALSE,
         fontsize = 6,
         #cutree_cols = 2,
         annotation_col = metadata_df,
         annotation_colors = list(
           DevStage = c("Adult" = "#009E73", "Embryonic" = "#CC79A7")
         )) %>%
  as.ggplot() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Heatmap (filtered)",
       subtitle = "(top 200 varying gene) \n")

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Heatmap_filtered.png",
       width = 5, height = 5, dpi = 300, unit = "in", bg = "white")







# Clear environment at the end of the session ----
rm(list = ls())
