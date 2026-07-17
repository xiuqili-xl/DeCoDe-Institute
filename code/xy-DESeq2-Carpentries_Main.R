# Goal ----
# Explore Carpentries lesson, RNA-seq analysis with Bioconductor
# Note, this lesson is in beta
# Also, we'll start with the Episode, Exploratory analysis and quality control
# This is where DESeq2 starts
# https://carpentries-incubator.github.io/bioc-rnaseq/04-exploratory-qc.html


# Load libraries ----
library(tidyverse)
library(here)
library(ggplotify)



# Import data ----
# The effect of upper-respiratory infection on transcriptomic changes in the CNS
# GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870
# Reconstruction of the count matrix from GEO: `xy-DESeq2-Carpentries-GEO-reconstruction.R`

# Actual data used: https://raw.githubusercontent.com/carpentries-incubator/bioc-rnaseq/main/episodes/data/GSE96870_counts_cerebellum.csv
# Downloaded locally for access and preservation

count_data <- read_csv(here("data_other", "GSE96870_counts_cerebellum.csv")) %>%
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
all(colnames(count_data) == rownames(metadata))

## look at a couple of variables 
unique(metadata$tissue)          # only cerebellum
unique(metadata$mouse) 
length(unique(metadata$mouse))   # 22, so unique id for animal


## experimental groups
metadata %>%
  dplyr::count(tissue, age, sex, infection, time, name = "no_mice")
## the two main variables for subsequent analysis are sex and time



# Assemble Summarized Experiment ----
library(SummarizedExperiment)

## make show samples and genes are in the same order!
all.equal(colnames(count_data), rownames(metadata))        # sample
all.equal(rownames(count_data), rownames(gene_info))       # gene

## create the SummarizedExperiment object
summarized_exp <- SummarizedExperiment(
  assays = list(counts = as.matrix(count_data)),
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




# RNA-Seq Exploratory Analysis ----
library(DESeq2)

## Read in RDS object 
## summarized_exp <- readRDS(here("data_other", "GEO96870_SummarizedExperiment.rds"))
## This is what Carpenty lesson did... but it seems too complexed...


## Construct DESeq object ----
## using DESeqDataSetFromMatrix, might simplify the workflow
metadata <- metadata %>%
  mutate(sex = factor(sex, levels = c("Female", "Male")),
         time = factor(time, levels = c("Day0", "Day4", "Day8")))

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              rowRanges = as(gene_info, "GRanges"),
                              design = ~ sex + time)
dds       # 41,786 genes 


## Filter out lowly expressed genes ----
dds <- dds[rowSums(counts(dds)) > 5, ]
dds       # 27,430 genes

head(assay(dds))
## shortcut to grab count matrix
# can be used to grab normalized count using `normalized = TRUE`

head(counts(dds))       
# generic Bioconductor function, default to the first slot (i.e., the raw count)
# can be used to access other data types stored in the object `assay(dds, "vst)` 

rowData(dds)      # gene_info
colData(dds)      # metadata


## Library size difference ----
dds$animalID <- paste(dds$sex, dds$time, dds$mouse, sep = "_")
dds$group <- paste(dds$sex, dds$time, sep = "_")
dds$libSize <- colSums(counts(dds))

colData(dds)       # added three columns to colData within dds (i.e., metadata)

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
## there is now a new column, sizeFactor -- to account for differences in seq depth
## if we later request counts(dds, normalized=TRUE), DESeq will divide the counts by size factor on the fly

head(counts(dds))      
head(count_data)
## running estimateSizeFactors() does NOT change the raw count value

head(counts(dds, normalized = TRUE)) 
## now we get normalized counts


## plot sizeFactor against libSize
ggplot(dat = colData(dds) %>% as.data.frame(),
       mapping = aes(x = libSize, y = sizeFactor, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "Library size", y = "Size factor")



## Transform data ----
## For analysis, we want the variance of gene expression to be relatively independent of
## the average value. For RNA-seq read counts, this is not the case...
library(vsn)               # Variance Stabilization and Calibration

## visualize the relationship between mean and sd across rows
meanSdPlot(counts(dds), ranks = FALSE)
## note, color of the dot here is a density heatmap. it indicates the number of genes 
## falling into that specific area of the plot. The highest density of genes in RNA-seq data 
## typically lives at low-to-moderate expression levels (bottom left), while only a few 
## highly expressed genes fan out toward the upper right.


## perform variance stablizing transformation
vsd <- vst(dds, blind = TRUE)
vsd

meanSdPlot(assay(vsd), ranks = FALSE)



## Sample-to-sample distance | Heatmap and clustering ----
library(RColorBrewer)
library(ComplexHeatmap)         # a Bioconductor package

dst <- dist(t(assay(vsd)))
color_palette <- colorRampPalette(brewer.pal(9, "Blues"))(255)

ComplexHeatmap::Heatmap(
  as.matrix(dst), 
  col = rev(color_palette),
  row_names_gp = gpar(fontsize = 8),     # row label font size
  column_names_gp = gpar(fontsize = 8), # column label font size
  name = "Euclidean\ndistance",
  cluster_rows = hclust(dst),
  cluster_columns = hclust(dst),
  top_annotation = columnAnnotation(
    sex = vsd$sex,
    time = vsd$time,
    col = list(sex = c(Female = "red", Male = "blue"),
               time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple")))
) %>%
  as.ggplot() +
  labs(title = "GSE96870: Effect of Infection on Cerebellum | Sample-to-Sample Distance")

ggsave(path = here("graphs_other"), filename = "GSE96870_Sample-to-Sample-Distance.png",
       width = 6.2, height = 5, dpi = 300, unit = "in", bg = "white")




## PCA  ----
plotPCA(vsd, intgroup = c("sex", "time"))

pca_data <- plotPCA(vsd, intgroup = c("sex", "time"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(data = pca_data,
       mapping = aes(x = PC1, y = PC2, color = time, shape = sex)) +
  geom_point(size = 3, alpha = 0.5) +
  scale_x_continuous(limits = c(-14, 14), expand = 0.01) +
  scale_y_continuous(limits = c(-14, 14), expand = 0.01) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  labs(title = "GSE96870: Effect of Infection on Cerebellum | PCA", 
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"))

ggsave(path = here("graphs_other"), filename = "GSE96870_PCA.png",
       width = 6, height = 4, dpi = 300, unit = "in")



# Differential Gene Expression ----
## Dispersion ----
dds <- estimateDispersions(dds)
dds               # estimateDispersions() adds a new assay called `mu`

head(assay(dds))
head(assay(dds, "counts"))      # this is the same as head(assay(dds))
head(assay(dds, "mu"))

plotDispEsts(dds)


## Testing ----
dds <- nbinomWaldTest(dds)
dds             # nbinomWaldTest() adds 2 new assays called `H` and `cooks`



## Run DESeq() directly -----
dds2 <- DESeqDataSetFromMatrix(countData = count_data,
                               colData = metadata,
                               rowRanges = as(gene_info, "GRanges"),
                               design = ~ sex + time)
dds2       # 41,786 genes 

dds2 <- dds2[rowSums(counts(dds2)) > 5, ]      # filtering
dds2       # 27,430 genes

dds2 <- DESeq(dds2)
dds2
## note, sizeFactor was added to colData
## 4 assays: counts, mu, H, cooks
## aka... running DESeq() is the same as running estimateSizeFactors(), estimateDispersions(), and nbinomWaldTest()

all(dds$sizeFactor == dds2$sizeFactor)
all(assay(dds, "mu") == assay(dds2, "mu"))


## Explore results for specific contrast ----
results(dds)         
# bc of how we set up our factors, defaulting to comparing Day 8 vs Day 0
# since we set out design = ~ sex + time, we are accounting for sex already

res_time1 <- results(dds, contrast = c("time", "Day8", "Day0"), alpha = 0.05)     # Day 0 is the reference
summary(res_time1)
# note, by default results() uses `independentFiltering = TRUE`, 
# which is removing lowly expressed genes, hence the low counts value in summary 

res_time2 <- results(dds, contrast = c("time", "Day4", "Day0"), alpha = 0.05) 
summary(res_time2)

res_sex <- results(dds, contrast = c("sex", "Male", "Female"), alpha = 0.1) 
summary(res_sex)


## genes with the lowest adj pvalue
head(res_time1[order(res_time1$padj), ])
head(res_sex[order(res_sex$padj), ])

## by default, DESeq() uses Benjamini-Hochberg procedure to calculate padj and
## correct for multiple hypothesis testing


## MA-plot ----
plotMA(res_time1, ylim = c(-3, 3))
plotMA(res_sex, ylim = c(-3, 3))

## plot MA is essentially plotting the first 2 columns of res
ggplot(data = as.data.frame(res_time1) %>%
         mutate(sig = if_else(padj < 0.05, "sig", "not-sig")),
       mapping = aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_y_continuous(limit = c(-3, 3), expand = 0.01, oob = scales::squish) +
  scale_color_manual(values = c("sig" = "blue", "not-sig" = "grey50"), na.value = "grey50") +
  theme_bw()


## from the MA graph, we can see that genes with low mean counts to have big log fold change
## this is bc lowly expressed genes tend to be very noisy
## we can shrink the log fold changes of genes with low mean and high dispersion
## as they contain littel information
library(apeglm)          # Bioconductor package, needed for lfcShrink

resultsNames(dds)        # use to find coef below
res_time_lfc <- lfcShrink(dds, coef = "time_Day8_vs_Day0", res = res_time1)

res_time_lfc
plotMA(res_time_lfc)

ggplot(data = as.data.frame(res_time_lfc) %>%
         mutate(sig = if_else(padj < 0.05, "sig", "not-sig")),
       mapping = aes(x = baseMean, y = log2FoldChange, color = sig)) +
  geom_point(size = 0.5) +
  scale_x_log10() +
  scale_y_continuous(limit = c(-3, 3), expand = 0.01, oob = scales::squish) +
  scale_color_manual(values = c("sig" = "blue", "not-sig" = "grey50"), na.value = "grey50") +
  theme_bw()



## Heatmap of DE genes ----
## get top DE genes
top_genes <- as.data.frame(res_time1) %>%
  arrange(padj) %>%
  head(10) %>%
  rownames()

top_genes

heatmap_data <- assay(vsd)[top_genes, ]

## reorder columns to match that in Carpentry lesson -- only clustering on gene
sample_order <- colData(vsd) %>%
  as.data.frame() %>%
  arrange(time, sex) %>% 
  rownames()

heatmap_data <- heatmap_data[ , sample_order]

## scale counts for visualization
scale(t(heatmap_data))           # calculates the z-score
heatmap_data <- t(scale(t(heatmap_data)))


## generate heatmp
ComplexHeatmap::Heatmap(
  matrix = heatmap_data,
  top_annotation = HeatmapAnnotation(
    df = colData(vsd)[sample_order, c("time", "sex")],
    col = list(sex = c(Female = "red", Male = "blue"),
               time = c(Day0 = "yellow", Day4 = "forestgreen", Day8 = "purple"))
  ),
  name = "z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE
) %>%
  as.ggplot() +
  labs(title = "GSE96870: Effect of Infection on Cerebellum | Heatmap")

ggsave(path = here("graphs_other"), filename = "GSE96870_Heatmap.png",
       width = 8, height = 4, dpi = 300, unit = "in", bg = "white")



## export DESeq results ----
as.data.frame(res_time1) %>% 
  rownames_to_column(var = "gene") %>%
  write_csv(here("outputs_other", "GSE96870_DESeq_results.csv"))
## could combine with gene_info too


# Clean environment before exiting ----
rm(list = ls())



