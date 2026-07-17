# Goal -----
# Follow DESeq2 tutorial to understand workflow
# Tutorial https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/


# Dataset used for tutorial ----
# From EBI Expression Atlas: E-GEOD-50760
# Expression Atlas https://www.ebi.ac.uk/gxa/experiments/E-GEOD-50760/Downloads
# BioStudies https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-50760
# note, fastq file available --- could investigate


# Load libraries ----
library(DESeq2)
library(tidyverse)
library(here)


# Investigate dataset ----
# the tutorial provides a link to the raw counts, so does Expression Atlas
# start by comparing the two
rawCounts_tutorial <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")
rawCounts_original <- read.delim(here("data_raw", "E-GEOD-50760", "E-GEOD-50760-raw-counts.tsv")) %>%
  select(Gene.ID, Gene.Name, sort(names(.)))

glimpse(rawCounts_tutorial)
glimpse(rawCounts_original)
# sampleData_tutorial has ~7k more rows than sampleData_original

rawCounts_tutorial$Gene.ID[1:50] == rawCounts_original$Gene.ID[1:50]
rawCounts_tutorial$Gene.Name[1:50] == rawCounts_original$Gene.Name[1:50]
# the first 50 rows have identical Gene.ID and Gene.Name

rawCounts_tutorial$SRR975551[1:50]
rawCounts_original$SRR975551[1:50]
# but the values are very different for this one sample tested...

sum(is.na(rawCounts_tutorial$Gene.ID))
sum(is.na(rawCounts_original$Gene.ID))
# there are no NA for Gene.Name in either dataset

length(unique(rawCounts_tutorial$Gene.ID)) == nrow(rawCounts_tutorial)
length(unique(rawCounts_original$Gene.ID)) == nrow(rawCounts_original)
# Gene.ID in both datasets are unique

setdiff(rawCounts_tutorial$Gene.ID, rawCounts_original$Gene.ID)
setdiff(rawCounts_original$Gene.ID, rawCounts_tutorial$Gene.ID)
# so there are gene Gene.ID that's unique to either dataset
# not sure what to make of this...
# for now... proceed with using sampleData_tutorial to replicate tutorial

rm(rawCounts_original, rawCounts_tutorial)


# Input data ----
# read in the raw read counts
rawCounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")
head(rawCounts)

# read in the sample mappings (quick inspection matches that directly from the repository)
sampleData <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")
head(sampleData)



# Create a DESeqDataSet object ----
## convert count data into a *matrix* of appropriate form that DESeq2 can read
## (1) the matrix only stores numeric data, 
## (2) its rownames are genomic features (i.e., gene)
## (3) its colnames are sample names
geneID <- rawCounts$Gene.ID
geneID

rawCountsMatrix <- as.matrix(rawCounts[ , 3:ncol(rawCounts)])
colnames(rawCountsMatrix)
rownames(rawCountsMatrix) <- geneID
head(rawCountsMatrix)


## convert sample variable metadata to an appropriate form that DESeq2 can read
## note, only keep variables we care about and they should be factors
rownames(sampleData) <- sampleData$Run

sampleDataKey <- sampleData %>%
  select(tissueType = Sample.Characteristic.biopsy.site.,
         individualID = Sample.Characteristic.individual.) %>%
  mutate(individualID = factor(individualID), 
         tissueType = case_when(tissueType == "normal" ~ "normal-looking surrounding colonic epithelium",
                                tissueType == "primary tumor" ~ "primary colorectal cancer",
                                tissueType == "colorectal cancer metastatic in the liver" ~ "metastatic colorectal cancer to the liver"), 
         # set the control group to be first elemet in `levels`
         tissueType = factor(tissueType, levels = c("normal-looking surrounding colonic epithelium", 
                                                    "primary colorectal cancer", 
                                                    "metastatic colorectal cancer to the liver"))
         )

head(sampleDataKey)
glimpse(sampleDataKey)


## reorder rawCountsMatrix
## for DESeq2 to work properly, the column names of the count matrix must be 
## in the same order as the row names of metadata df 
rawCountsMatrix <- rawCountsMatrix[ , rownames(sampleDataKey)]
all(colnames(rawCountsMatrix) == rownames(sampleDataKey))


## finally, create the DESeqDataSet Object
deseq2Data <- DESeqDataSetFromMatrix(countData = rawCountsMatrix, 
                                     colData = sampleDataKey, 
                                     design= ~ individualID + tissueType)
deseq2Data



# Pre-filtering of data ----
## not necessary, but will reduce the size of DESeqDataSet object and speed up runtime
## here, we'll filter for genes that have more tan a sum total of 5 reads 

## let's first investigate the effect
dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])      # reduce no. rows from 65k to 35k

## apply the filter
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]



# Set up multi-core (optional)
library(BiocParallel)
register(MulticoreParam(4))
# to find out the number of cpu, use Terminal `sysctl -n hw.ncpu`
# when calling DESeq() and results(), add the argument `parallel = TRUE`



# Differential Expression Analysis ----
## run DESeq() - this might take a couple of minutes, parallel processing should hopeful speed things up
deseq2Data <- DESeq(deseq2Data, parallel = TRUE)


## extracting results with results()
## when we have more than 2 groups, we want to specify that comparisons to make
## `contrast = c(<variable>, <numerator aka comparator>, <denominator aka control>)`
deseq2Results <- results(deseq2Data, 
                         # for tissueType, compare primary cancer vs normal tissue
                         contrast = c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"),
                         parallel = TRUE)
summary(deseq2Results)


## MA-plot ----
## graphs log fold change against mean of normalized count
## in general, we expect expression of genes to remain consistent between conditions
## so MA plot should be simiar to the shape of a trumpet, with most point resideing around y = 0
plotMA(deseq2Results)


### we can also do this manually using ggplot
library(scales)       # for `oob = squish`
library(viridis)

deseq2ResDF <- as.data.frame(deseq2Results) %>%
  mutate(significant = if_else(padj < 0.1, "sig", "not-sig"))

head(deseq2ResDF)

ggplot(data = deseq2ResDF,
       mapping = aes(x = baseMean, y = log2FoldChange, 
                     color = significant, fill = significant)) +
  geom_point(shape = 21, size = 1, alpha = 0.3) +
  geom_hline(yintercept = 0, color = "tomato", linewidth = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-3, 3), oob = squish) +
  scale_color_manual(name = "q-value", values = c("sig" = "blue", "not-sig" = "grey"), 
                     na.value = "grey") +
  scale_fill_manual(name = "q-value", values = c("sig" = "blue", "not-sig" = "grey"), 
                     na.value = "grey") +
  theme_bw() +
  labs(title = "DESeq2 Tutorial | MA-plot",
       subtitle = "Colorectal Cancer Data (E-GEOD-50760)",
       x = "mean of normalized counts", y = "log fold change") 

ggsave(path = here("graphs", "E-GEOD-50760"), filename = "E-GEOD-50760_MA-plot.png",
       width = 6, height = 4, dpi = 300, units = "in")


## we could add a density map
ggplot(data = deseq2ResDF,
       mapping = aes(x = baseMean, y = log2FoldChange, color = padj)) +
  geom_point(size = 1) +
  geom_hline(yintercept = 0, color = "darkorchid4", linewidth = 1, linetype = "longdash") + 
  geom_density_2d(colour = "black", size = 1) +
  scale_x_log10() +
  scale_y_continuous(limits = c(-3, 3), oob = squish) +
  scale_colour_viridis(direction = -1, trans='sqrt') +
  scale_fill_manual(name = "q-value", values = c("sig" = "blue", "not-sig" = "grey"), 
                    na.value = "grey")  +
  theme_bw() +
  labs(title = "DESeq2 Tutorial | MA-plot",
       subtitle = "Colorectal Cancer Data (E-GEOD-50760)",
       x = "mean of normalized counts", y = "log fold change")

ggsave(path = here("graphs", "E-GEOD-50760"), filename = "E-GEOD-50760_MA-plot_density.png",
       width = 6, height = 4, dpi = 300, units = "in")



## Graph normalized count for a single geneID ----
## extract counts for the gene otop2
otop2Counts <- plotCounts(deseq2Data, 
                          gene = "ENSG00000183034", 
                          intgroup = c("tissueType", "individualID"), 
                          returnData = TRUE) %>%
  mutate(tissueType = fct_recode(tissueType, 
                                 "normal epithelium" = "normal-looking surrounding colonic epithelium",
                                 "primary tumor" = "primary colorectal cancer",
                                 "liver metastasis" = "metastatic colorectal cancer to the liver"))

otop2Counts

ggplot(data = otop2Counts, 
       mapping = aes(x = tissueType, y = count, 
                     colour = individualID, group = individualID)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  guides(colour = guide_legend(ncol = 2)) +
  labs(title = "DESeq2 Tutorial | Normalized Count of OTOP2",
       subtitle = "Colorectal Cancer Data (E-GEOD-50760)") 
  
ggsave(path = here("graphs", "E-GEOD-50760"), 
       filename = "E-GEOD-50760_Normalized Count of OTOP2.png",
       width = 6, height = 4, dpi = 300, units = "in")
  
  
## explore the stats and raw values for this gene
deseq2ResDF["ENSG00000183034", ]

rawCounts %>%
  filter(Gene.ID == "ENSG00000183034") %>%
  pivot_longer(cols = -c(Gene.ID, Gene.Name), names_to = "Run", values_to = "expression") %>%
  inner_join(select(sampleData, Run, tissueType = Sample.Characteristic.biopsy.site., individualID = Sample.Characteristic.individual.),
             by = "Run") %>%
  view()



## Plot heatmap ----
## plot heatmap of differentially expressed genes and perform unsupervised clustering

## transform count data using variance stabilizing transform
deseq2VST <- vst(deseq2Data)

## convert deseq2VST, which is a DESeqTransform object into a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

## keep only the significantly differentiated genes where the fold change was at least 3
deseq2ResDF %>%
  filter(padj < 0.05, log2FoldChange > 3)

deseq2VST_sig <- deseq2VST %>%
  filter(Gene %in% rownames(filter(deseq2ResDF, padj < 0.05, log2FoldChange > 3)))

## turn into long format for graphing
deseq2VST_sig_long <- deseq2VST_sig %>%
  pivot_longer(cols = -Gene, names_to = "Run", values_to = "value")

ggplot(data = deseq2VST_sig_long, 
       mapping = aes(x = Run, y = Gene, fill = value)) +
  geom_raster() +
  scale_fill_viridis(trans="sqrt") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


## Clustering ----
## compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(select(deseq2VST_sig, -Gene))
distanceSample <- dist(t(select(deseq2VST_sig, -Gene)))

## cluster based on the distance calculation
clusterGene <- hclust(distanceGene, method = "average")
clusterSample <- hclust(distanceSample, method = "average")

## construct a dendogram for samples
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))

sampleDendrogram <- ggplot(sampleDendrogramData) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_dendro()

sampleDendrogram


## recreate the heatmap, this time ordering the samples by clustering results
clusterSample$labels
clusterSample$order

clusterSample_reorderd <- clusterSample$labels[clusterSample$order]

heatmap1 <- deseq2VST_sig_long %>% 
  mutate(Run = factor(Run, levels = clusterSample_reorderd)) %>%
  ggplot(mapping = aes(x = Run, y = Gene, fill = value)) +
  geom_raster() +
  scale_fill_viridis(trans = "sqrt") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

heatmap1


## combine dendrogram with heatmap 
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap1, ncol = 1, heights = c(1, 5))


## there is probably an easier way to do this....
## Let's try using `pheatmap`
## following process from https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
library(pheatmap)

deseq2VST_sig2 <- select(deseq2VST_sig, -Gene)

pheatmap(deseq2VST_sig2,
         show_rownames = FALSE,
         fontsize = 6) 

heatmap_annotation_col <- sampleDataKey %>%
  mutate(tissueType = fct_recode(tissueType, 
                                 "normal epithelium" = "normal-looking surrounding colonic epithelium",
                                 "primary tumor" = "primary colorectal cancer",
                                 "liver metastasis" = "metastatic colorectal cancer to the liver"))

tumor_heatmap <- pheatmap(deseq2VST_sig2,
         color = viridis(50),
         clustering_method = "average",    # to match that method in tutorial
         show_rownames = FALSE,
         fontsize = 6,
         #cutree_cols = 2,
         annotation_col = heatmap_annotation_col,
         annotation_colors = list(
           tissueType = c("normal epithelium" = "lightgreen", "primary tumor" = "salmon", "liver metastasis" = "darkred")
         )) 

tumor_heatmap

png(filename = here("graphs", "E-GEOD-50760", "E-GEOD-50760_HeatmapVST-1.png"), 
    width = 8, height = 6, res = 300, units = "in")
tumor_heatmap
dev.off()
  

## explore additinal functionalities
## following instructions from https://bioinformatics.ccr.cancer.gov/docs/data-visualization-with-r/Lesson5_intro_to_ggplot/
library(ggplotify)

pheatmap(deseq2VST_sig2,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         scale = "row",            # scale by gene
         show_rownames = FALSE,
         fontsize = 6,
         #cutree_cols = 2,
         annotation_col = heatmap_annotation_col,
         annotation_colors = list(
           tissueType = c("normal epithelium" = "lightgreen", "primary tumor" = "salmon", "liver metastasis" = "darkred")
         )) %>%
  as.ggplot()

ggsave(path = here("graphs", "E-GEOD-50760"), filename = "E-GEOD-50760_HeatmapVST-2.png",
       width = 8, height = 6, dpi = 300, units = "in", bg = "white")


# Remaining questions ----
## Is there a way to export all normalized data for all genes after DESeq()?
## What dataset should be used to generate the heatmap?

## https://bioconductor.statistik.tu-dortmund.de/packages/3.5/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-quality-assessment-by-sample-clustering-and-visualization
## generates heatmap of different transformation of the data...

deseq2Data

deseq2Counts <- counts(deseq2Data, normalized = TRUE)
deseq2Counts
# normalize for sequencing depth and composition bias between samples
# does not normalize for gene length, meaning 
# the resulting values are comparable across samples for the same gene, 
# but not between different genes within the same sample
# Need to calculate TPM or FPKM for within-sample comparison


# Clear environment at the end of the session ----
rm(list = ls())
