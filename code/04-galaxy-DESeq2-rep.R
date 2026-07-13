# Goal -----
# Replicate DESeq2 analysis that was performed on Galaxy during DeCoDe Institute

# Useful resources ----
# Full vignette: https://bioconductor.statistik.tu-dortmund.de/packages/3.5/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Tufts quick guide https://rtguides.it.tufts.edu/bio/tutorials/de-seq2-in-r-ood.html


# Dataset ----
# Transcriptome of Embryonic and Adult mouse cerebral cortex
# Dataset for teaching https://zenodo.org/records/20531535
# Note FM1 was not used in the Galaxy analysis; there was something wrong with it

# Load library ----
library(tidyverse)
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
## presumably they all have the same feature


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
count_matrix <- count_combined[ , 2:ncol(count_combined)]
rownames(count_matrix) <- count_combined$Feature
count_matrix <- as.matrix(count_matrix)

head(count_matrix)
class(count_matrix)
## to investigate for the future: do we really need a matrix or is a df fine


# Create metadata df ----
metadata_df <- data.frame(Sample = colnames(count_matrix)) %>%
  mutate(DevStage = case_when(str_detect(Sample, "AM") ~ "Adult",
                              str_detect(Sample, "FM") ~ "Embryonic"),
         DevStage = factor(DevStage, levels = c("Adult", "Embryonic")))

rownames(metadata_df) <- metadata_df$Sample
metadata_df <- metadata_df %>% select(-Sample)

metadata_df
glimpse(metadata_df)


## check the col of count_matrix is in the same order as rows of metadata_df
colnames(count_matrix) == rownames(metadata_df)



# Construct DESeqDataSet ----
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata_df,
                              design = ~ DevStage)
dds

## normally, we should prefilter to remove rows in which there are very few reads, thus
## reducing the required memory, and increasing the speed.
## prefiltering can also improve visualizations, as features with no information 
## for differential expression are not plotted.
## dds <- dds[rowSums(counts(dds)) > 10, ]

## However, the default in Galaxy does not perform prefiltering, so we won't here


# Differential Expression Analysis ----
dds <- DESeq(dds)

## Results ----
d_results <- results(dds, contrast = c("DevStage", "Adult", "Embryonic"))

summary(d_results)
d_results

d_results_df <- as.data.frame(d_results)
head(d_results_df)
## spot check matches galaxy output!


## Normalized counts ----
d_counts <- counts(dds, normalized = TRUE)

view(d_counts)
d_counts["Samd10", ]
## spot check matches galaxy output!



# MA-plot ----
plotMA(d_results)
plotMA(d_results, ylim = c(-13, 13), size = 1) 


# Plot dispersion ----
plotDispEsts(dds)


# PCA ----
d_vsd <- vst(dds, blind = TRUE)           # Variance-stabilizing transformation
plotPCA(d_vsd, intgroup = "DevStage") 

d_rld <- rlog(dds, blind = TRUE)          # Regularized log transformation
plotPCA(d_rld, intgroup = "DevStage") 
## these two looks more similar than galaxy output

PCA_data <- plotPCA(d_rld, intgroup = "DevStage", returnData = TRUE) 
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


# Sample-to-Sample Distance ----
## extact transformed matrics and transpose it
## use vst for very large datasets
d_rld_matrix <- assay(d_rld)
sampleDist <- dist(t(d_rld_matrix))
sampleDist

## convert to a distance matrix
sampleDistMatrix <- as.matrix(sampleDist)
sampleDistMatrix

# create a sample-to-sample distance heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(30)) %>%
  as.ggplot() +
  labs(title = "Mouse Adult vs Embryonic Cortex | Sample Distance Heatmap")
## visually identical to galaxy output

ggsave(path = here("graphs_other"), filename = "MouseDevCortex_Sample-Distance-Heatmap.png",
       width = 6, height = 6, dpi = 300, unit = "in", bg = "white")



# Volcano Plot ----





# Heatmap ----


# Clear environment at the end of the session ----
rm(list = ls())
