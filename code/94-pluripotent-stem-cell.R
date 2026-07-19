# Goal -----
# Analyze read counts from  GEO https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80550
# Which is the same .BAM dataset we were looking at https://doi.org/10.5281/zenodo.15359156 
# If all goes well, we can try to replicate some of the figures in https://pubmed.ncbi.nlm.nih.gov/27350605/ 


# Load libraries -----
library(tidyverse)
library(ggplot2)
library(here)

## specific libraries used for analysis are loaded as they are needed


# Import and wrangle data ----
## readcounts.txt downloaded from GEO "contains the absolute number of reads per gene for each sample"
read_counts <- read_delim(file = here("data_raw", "PluripotentStemCell", "GSE80550_readcounts.txt.gz"),
                          delim = "\t")
head(read_counts)
colnames(read_counts)
## so... this df contains both the read count and gene info


## gene info ----
gene_info <- read_counts %>% 
  select(-starts_with("3430_")) %>%
  dplyr::rename(seqnames = chromosome_name, start = start_position, end = end_position) %>%
  column_to_rownames(var = "ensembl_gene_id")

head(gene_info)


## read counts ----
read_counts <- read_counts %>%
  select(ensembl_gene_id, starts_with("3430_")) %>%
  column_to_rownames(var = "ensembl_gene_id")

head(read_counts)


## experimental info ----
sra_run_table <- read_csv(file = here("data_raw", "PluripotentStemCell", "GSE80550_SraRunTable.csv"),
                          name_repair = "universal") 

series_matrix <- read_delim(file = here("data_raw", "PluripotentStemCell", "GSE80550_series_matrix.txt.gz"),
                            delim = "\t", skip = 38) %>% 
  rename_with(~ "Variable", 1)

series_matrix_cleaned <- series_matrix %>%
  pivot_longer(cols = -Variable, names_to = "Sample", values_to = "value") %>%
  mutate(Variable = str_remove(Variable, "!Sample_")) %>%
  filter(Variable %in% c("geo_accession", "characteristics_ch1")) %>%
  pivot_wider(id_cols = Sample, names_from = Variable, values_from = value) %>%
  mutate(geo_accession = map_chr(geo_accession, 1),
         strain = map_chr(characteristics_ch1, 1),
         treatment = map_chr(characteristics_ch1, 2),
         days_treated = map_chr(characteristics_ch1, 3),
         cell_type = map_chr(characteristics_ch1, 4),
         knockdown = map_chr(characteristics_ch1, 5)) %>%
  select(-characteristics_ch1)


## metadata ----
meta_data <- data.frame(sample = colnames(read_counts)) %>%
  mutate(sample = str_replace(sample, "3430_", "Sample"),
         knockdown = case_when(str_detect(sample, "Scrambled") ~ "scrambled",
                               str_detect(sample, "Trim28") & str_detect(sample, "Setdb1") ~ "Trim28.Setdb1",
                               str_detect(sample, "Trim28") ~ "Trim28"),
         induced = case_when(str_detect(sample, "-DOX_") ~ "nodox",
                             TRUE ~ "dox"),
         group = paste(induced, knockdown, sep = "_"),
         group = factor(group, levels = c("dox_scrambled", "dox_Trim28", "dox_Trim28.Setdb1", "nodox_scrambled")), 
         sample = str_remove(sample, "_.*")) %>%
  column_to_rownames(var = "sample")

meta_data


## create read_counts_cleaned 
read_counts_cleaned <- read_counts
colnames(read_counts_cleaned)

colnames(read_counts_cleaned) <- rownames(meta_data)
head(read_counts_cleaned)


# Input data for DE analysis -----
library(DESeq2)

## Construct DESeq Object ----
dds <- DESeqDataSetFromMatrix(countData = read_counts_cleaned,
                              colData = meta_data,
                              rowRanges = as(gene_info, "GRanges"),
                              design = ~ group)
dds                  # 24,207 genes x 12 samples
dds$group            # using dox_scrambled as reference


## Pre-filtering ----
## let's first investigate the sum of genes across samples
rowSums(counts(dds))
rowSums(counts(dds)) %>% summary()

ggplot(data = NULL,
       mapping = aes(x = rowSums(counts(dds)))) +
  geom_histogram() +
  scale_y_continuous(expand = 0.01) +
  theme_bw() +
  labs(title = "Sum of reads for each gene across all samples")

ggplot(data = NULL,
       mapping = aes(x = rowSums(counts(dds)))) +
  geom_histogram() +
  scale_x_continuous(limits = c(-1, 100), breaks = seq(0, 100, by = 20), expand = 0.01) +
  scale_y_continuous(expand = 0.01) +
  theme_bw() +
  labs(title = "Sum of reads for each gene across all samples (zoomed in)")


## set a threshold 10 reads across all samples
dds <- dds[rowSums(counts(dds)) > 10, ]
dds                   # 18,176 genes



## Library size difference ----
colSums(counts(dds))                    # note colSums() is a base R function

dds$libSize <- colSums(counts(dds)) 
dds                                     # libSize has been added to colData

colData(dds) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>%
  mutate(order_no = str_remove(sample, "Sample"),
         order_no = as.numeric(order_no)) %>%
  ggplot(mapping = aes(x = reorder(sample, order_no), y = libSize / 1e6, fill = group)) +
  geom_col() +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(limits = c(0, 20), expand = 0.01, name = "Library size (million reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Total number of reads for each sample")

## confirming there are difference in library size


## Calculate size factor ----
dds <- estimateSizeFactors(dds)
dds         # check that sizeFactor has been added to colData

ggplot(data = colData(dds) %>% as.data.frame() %>% rownames_to_column("sample"),
       mapping = aes(x = libSize / 1e6, y = sizeFactor, color = group, fill = group)) +
  geom_point(size = 3, shape = 21, alpha = 0.5) +
  theme_bw() +
  labs(title = "Correlation between library size and estimateSizeFactors", 
       x = "Library size (million reads)", y = "Size factor")



## Explore normalized counts and explort ----
## this part is probably not needed for future analysis...
read_counts_norm <- counts(dds, normalized = TRUE) %>%
  as.data.frame()

head(read_counts_norm)

data.frame(libSize = colSums(read_counts_norm)) %>% 
  bind_cols(meta_data) %>%
  rownames_to_column(var = "sample") %>%
  mutate(order_no = str_remove(sample, "Sample"),
         order_no = as.numeric(order_no)) %>%
  ggplot(mapping = aes(x = reorder(sample, order_no), y = libSize / 1e6, fill = group)) +
  geom_col() +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(limits = c(0, 20), expand = 0.01, name = "Library size (million reads)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Total number of norm data for each sample")


read_counts_norm %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  write_csv(here("outputs_csv", "PluripotentStemCell", "GSE80550_norm_data.csv"))


# Exploratory analysis ----
## Use transformed data for exploratory analysis (like PCA, clustering, or heatmaps)
## but we need to use raw counts for DESeq()
library(ComplexHeatmap)
library(ggplotify)
library(RColorBrewer)


## Transform the data ----
## visualize sd and average read count 
library(vsn)
meanSdPlot(assay(dds), rank = FALSE)       # color of the dot here is a density heatmap

## let's apply the variance stabilizing transformation
## re vst() and rlog() --- use vst() as the default, esp with large number of samples
## rlog() may be desirable if sample size is small and library sizes are drastically different
## the way to check is meanSdPlot()

vsd <- vst(dds, blind = TRUE)
meanSdPlot(assay(vsd), ranks = FALSE)

rld <- rlog(dds, blind = TRUE)
meanSdPlot(assay(rld), ranks = FALSE)
## rld is slower


## Sample distances ----
dst <- dist(t(assay(vsd)))

Heatmap(
  as.matrix(dst), 
  col = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(30),
  name = "Distance",
  cluster_rows = hclust(dst),
  cluster_columns = hclust(dst),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = columnAnnotation(
    group = vsd$group,
    col = list(group = c(dox_scrambled = "#D55E00", dox_Trim28  = "#E69F00",
                         nodox_scrambled = "#000000", dox_Trim28.Setdb1 = "#F0E442"))
    )
) %>%
  as.ggplot() +
  labs(title = "Heatmap showing distances between all samples")

ggsave(path = here("graphs", "PluripotentStemCell"), filename = "GSE80550_Sample-Distance.png",
       width = 5.8, height = 4, dpi = 300, unit = "in", bg = "white")


## PCA ----
plotPCA(vsd, intgroup = "group")

pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_x_continuous(limits = c(-25, 45), expand = 0.01) +
  scale_y_continuous(limits = c(-25, 45), expand = 0.01) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "PCA of all samples", 
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"))

ggsave(path = here("graphs", "PluripotentStemCell"), filename = "GSE80550_PCA.png",
       width = 6, height = 4, dpi = 300, unit = "in")




# Differential Expression Analysis ----
## DESeq() requires raw, un-transformed count data (integers) as input. 
## It performs its own internal size factor normalization and uses the 
## negative binomial distribution to model the raw counts.


## Run DESeq() ----
dds <- DESeq(dds)
dds

dds$group

rowData(dds) %>% colnames()
## looks like there are three comparisons: dox_scrambled vs the other three groups

## the effect of dox OSKM induction 
res_dox <- results(dds, contrast = c("group", "dox_scrambled", "nodox_scrambled"), 
                   alpha = 0.05) 
res_dox
summary(res_dox)

write_csv(as.data.frame(res_dox) %>% rownames_to_column("ensembl_gene_id") %>% arrange(padj), 
          file = here("outputs_csv", "PluripotentStemCell", "GSE80550_DE_cox.csv"))

## what if we reverse the 2 groups in contrast
res_dox_rev <- results(dds, contrast = c("group", "nodox_scrambled", "dox_scrambled"), 
                       alpha = 0.05) 
summary(res_dox_rev)        ## genes up are now down, as expected


## the effect of Trim28 KD under dox OSKM induction 
res_trim28 <- results(dds, contrast = c("group", "dox_Trim28", "dox_scrambled"), 
                      alpha = 0.05) 
res_trim28
summary(res_trim28)

write_csv(as.data.frame(res_trim28) %>% rownames_to_column("ensembl_gene_id") %>% arrange(padj), 
          file = here("outputs_csv", "PluripotentStemCell", "GSE80550_DE_trim28.csv"))



## Look at dispersion ----
## note, bc we already ran DESeq(), we don't need to do run estimateDispersions()
## the latter is already included in DESeq()
plotDispEsts(dds)


## Log-fold shrinkage ----
## first look at MA-plots, this is essentially plotting the baseMean and the log2FoldChange colomns
## Genes with low expression tend to have bigger fold change... but this is bc of noise
plotMA(res_dox)
plotMA(res_trim28, ylim = c(-6, 6))        # note, Figure 3B in the manuscript, suggest shrinkage

## lfcShrink() looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution
resultsNames(dds)                          # name of coef to shrink, default type = apeglm requires coef
res_trim28_lfc <- lfcShrink(dds, 
                            coef = "group_dox_Trim28_vs_dox_scrambled", 
                            res = res_trim28)
res_trim28
res_trim28_lfc                             # change is happening in to log2FoldChange column

plotMA(res_trim28_lfc, ylim = c(-6, 6))



# Visualize results ----
## Heatmap of top DE genes ----
## the following lines of code are adapted from the Carpentry lesson 
## and uses the transformed counts vsd

## top DE genes
res_trim28_diff <- res_trim28 %>%
  as.data.frame() %>%
  filter(log2FoldChange > 2 | log2FoldChange < -2 ) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

nrow(res_trim28_diff)       # how many number of genes

res_trim28_diff_genes <- rownames(res_trim28_diff)

## subset vsd data
res_trim28_heatmap_data <- assay(vsd)[res_trim28_diff_genes, vsd$group %in% c("dox_scrambled", "dox_Trim28")]
head(res_trim28_heatmap_data)

## scale counts and clean up for visualization
res_trim28_heatmap_data_scaled <- t(scale(t(res_trim28_heatmap_data))) 

res_trim28_heatmap_data_cleaned <- res_trim28_heatmap_data_scaled %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl") %>%
  left_join(gene_info %>% rownames_to_column("ensembl") %>% select(ensembl, external_gene_id),
            by = "ensembl") %>%
  select(-ensembl) %>%
  column_to_rownames("external_gene_id")

res_trim28_heatmap_annot <- colData(vsd)[vsd$group %in% c("dox_scrambled", "dox_Trim28"), c("group")] 
res_trim28_heatmap_annot

## generate heatmap of DE genes
Heatmap(
  res_trim28_heatmap_data_cleaned,
  name = "z-score",
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(
    group = res_trim28_heatmap_annot,
    col = list(
      group = c(dox_scrambled = "#D55E00", dox_Trim28  = "#E69F00")
    )
  )
) %>%
  as.ggplot() +
  labs(title = "Heatmap showing z-score of differentially expressed genes",
       subtitle = "(abs log2FoldChange > 2 and adj p-value < 0.05)")

ggsave(path = here("graphs", "PluripotentStemCell"), 
       filename = "GSE80550_DE-Genes-Heatmap.png",
       width = 6, height = 12, dpi = 300, unit = "in", bg = "white")



## Volcanoe plots ----
library(plotly)
library(htmlwidgets)

## let's work with res_dox, bc this is what we have BAM file for at 
## https://doi.org/10.5281/zenodo.15359156 
sra_run_table %>%
  select(Run, Experiment, GEO_Accession = GEO_Accession..exp., knockdown, treatment)

## inspect res_dox
res_dox
summary(res_dox)
plotMA(res_dox, ylim = c(-10, 15))            # looks like cell is really ramping up expression...

## do some log fold shrinkage
resultsNames(dds)

res_dox_lfc <- lfcShrink(dds, type = "ashr",            # we'll try ashr with contrast 
                         contrast = c("group", "dox_scrambled", "nodox_scrambled"))
res_dox_lfc
plotMA(res_dox_lfc, ylim = c(-10, 15))   

dox_graph <- res_dox_lfc %>%
  as.data.frame() %>%
  rownames_to_column(var = "ensembl") %>%
  left_join(gene_info %>% rownames_to_column("ensembl") %>% select(ensembl, gene_id = external_gene_id),
            by = "ensembl") %>%
  filter(baseMean >= 100)           # cut down on the number of genes to graph to speed things up

dox_graph_subset1 <- dox_graph %>%
  filter(gene_id %in% c("Tfrc", "Trim28", "Pou5f1", "Sox2", "Klf4", "Mycbp")) %>%
  mutate(sig = (padj < 0.05))
# Oct4 is POU5F1 https://www.ncbi.nlm.nih.gov/datasets/gene/5460/ and mouse ortholog is Pou5f1

dox_graph_subset2 <- dox_graph %>%
  filter(gene_id %in% c("Actb", "Tuba1a", "Rpl13a", "Hprt")) %>%
  mutate(sig = (padj < 0.05))
# Tuba1a is tubulin alpha; Rpl13a is part of ribosomes; Hprt in a enzyme in purine salvage
  
dox_volcano <- ggplot(data = dox_graph, 
                      mapping = aes(x = log2FoldChange, y = -log10(padj),
                                    text = paste("Gene:", gene_id,
                                                 "\nlog2FC:", log2FoldChange,
                                                 "\npadj:", padj))
                      ) +
  geom_point(color = "grey60", fill = "grey80", alpha = 0.6, shape = 21, size = 3) +
  geom_point(data = dox_graph_subset2, 
             color = "royalblue3", fill = "royalblue3", shape = 21, size = 3) +
  geom_point(data = dox_graph_subset1, 
             color = "red3", fill = "red3", shape = 21, size = 3) +
  scale_x_continuous(limits = c(-12, 12), oob = scales::squish) +
  theme_bw()

dox_volcano_plotly <- ggplotly(dox_volcano,
                               tooltip = c("text"))
dox_volcano_plotly

saveWidget(dox_volcano_plotly, 
           file = here("graphs", "PluripotentStemCell", "GSE80550_Volcano_dox-treatment.html"),
           selfcontained = TRUE, 
           title = "Differentially expressed gene upon DOX OSKM induction")



# Clean environment before exiting ----
rm(list = ls())
