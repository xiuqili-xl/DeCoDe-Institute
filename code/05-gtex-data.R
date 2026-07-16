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
train_data <- read_delim(here("data_raw", "gtex_tissue_train_200x100.tsv"))
test_data <- read_delim(here("data_raw", "gtex_tissue_test_labeled_100x100.tsv"))


# Explore data ----
glimpse(train_data)                # columns: SampleID, Tissue, and 100 genes

unique(train_data$SampleID)        # 200 unique Sample ID
## named GTEX-<DonorID>-<TissueID>

unique(train_data$Tissue)          # 10 different types of tissues
train_data %>% count(Tissue)       # 20 donors per tissue


# Galaxy activities ----
## Datamash ---- 
## mean of RBPJL by tissue
train_data %>%
  group_by(Tissue) %>%
  summarise(RBPJL_ave = mean(RBPJL)) %>%
  ungroup()

## mean and stdev of a couple of genes
train_data_summary_long <- train_data %>%
  pivot_longer(cols = -c("SampleID", "Tissue"), names_to = "Gene", values_to = "TPM") %>%
  group_by(Gene, Tissue) %>%
  summarise(TPM_ave = mean(TPM),
            TMP_sd = sd(TPM)) %>%
  ungroup()

### first three genes
train_data_summary_long %>%
  filter(Gene %in% colnames(train_data)[3:5]) %>%
  rename(ave = TPM_ave, sd = TMP_sd) %>%
  pivot_wider(id_cols = Tissue, names_from = Gene, 
              values_from = c("ave", "sd"), names_glue = "{Gene}_{.value}") %>%
  select(Tissue, starts_with("RBPJL"), starts_with("ENSG00000261012"), starts_with("FAM83C"))


## Heatmap ---- 
library(pheatmap)
library(ggplotify)

## create data matrix
train_data_graph <- train_data %>%
  select(-Tissue) %>%
  column_to_rownames(var = "SampleID") %>%
  t()

## create metadata
unique(train_data$Tissue)

train_data_graph_metadata <- train_data %>%
  select(SampleID, Tissue) %>%
  mutate(Tissue = case_when(Tissue == "Adipose_Subcutaneous" ~ "Adipose",
                            Tissue == "Brain_Cortex" ~ "Brain",
                            Tissue == "Heart_Left_Ventricle" ~ "Heart",
                            Tissue == "Muscle_Skeletal" ~ "Muscle",
                            Tissue == "Skin_Sun_Exposed_Lower_leg" ~ "Skin",
                            Tissue == "Whole_Blood" ~ "Blood",
                            TRUE ~ Tissue)) %>%
  column_to_rownames(var = "SampleID")


pheatmap(mat = train_data_graph,
         color=colorRampPalette(c("navy", "white", "red"))(50),
         fontsize_row = 3,
         show_colnames = FALSE,
         annotation_col = train_data_graph_metadata,
         annotation_colors = list(
           Tissue = c("Adipose" = "#000000", "Blood" = "#E69F00", "Brain" = "#56B4E9", 
                      "Heart" = "#009E73", "Liver" = "#F0E442", "Lung" = "#0072B2", 
                      "Muscle" = "#D55E00", "Pancreas" = "#CC79A7", "Skin" = "grey60",
                       "Spleen" = "white")
         )
         ) %>%
  as.ggplot() +
  labs(title = "GTEx Data Subset | Heatmap",
       subtitle = "(replicating analysis in Galaxy)\n")

ggsave(path = here("graphs_other"), filename = "GETxSubset_Heatmap.png",
       width = 10, height = 6, dpi = 300, unit = "in", bg = "white")

## note, further customization might require `ComplexHeatmap`
## https://stackoverflow.com/questions/63290621/can-i-change-the-column-names-to-the-top-in-pheatmap



## PCA ---- 
train_data_pca <- train_data %>%
  select(-Tissue) %>%
  column_to_rownames(var = "SampleID")

head(train_data_pca)

## perform PCA without scaling
train_data_pca_res <- prcomp(train_data_pca)         # not scaled...
train_data_pca_res
summary(train_data_pca_res)

## extract component scores and form df for plotting
train_data_pca_scores <- as.data.frame(train_data_pca_res$x) %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(train_data %>% select(SampleID, Tissue), by = "SampleID") %>%
  select(SampleID, Tissue, PC1, PC2, PC3) 


## create 3D graph using plotly
library(plotly)
library(htmlwidgets)

hover_labels <- paste(
  "<b>SampleID:</b>", train_data_pca_scores$SampleID, "<br>",
  "<b>Tissue:</b>", train_data_pca_scores$Tissue, "<br>",
  "<b>PC1:</b>", round(train_data_pca_scores$PC1, 2), "<br>",
  "<b>PC2:</b>", round(train_data_pca_scores$PC2, 2), "<br>",
  "<b>PC3:</b>", round(train_data_pca_scores$PC3, 2)
)


train_data_pca_plotly <- plot_ly(
  data = train_data_pca_scores,
  x = ~PC1, 
  y = ~PC2, 
  z = ~PC3, 
  color = ~Tissue,
  colors = c("#000000", "#E69F00", "#56B4E9","#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "grey60", "white"),
  type = "scatter3d", 
  mode = "markers",
  text = hover_labels,   # Inject your custom text vector
  hoverinfo = "text"     # Force plotly to ONLY display your custom text
) %>%
  layout(
    showlegend = FALSE, 
    scene = list(
      xaxis = list(title = 'PC1'),
      yaxis = list(title = 'PC2'),
      zaxis = list(title = 'PC3')
    )
  )

train_data_pca_plotly

saveWidget(train_data_pca_plotly, 
           file = here("graphs_other", "GETxSubset_PCA_plotly.html"),
           selfcontained = TRUE, 
           title = "GTEx Data Subset | PCA")

## this line of code often creates a temporary folder containing the JavaScript and CSS 
## assets during the build process. Since we have `selfcontained = TRUE`, once the process 
## finishes, we don't need to keep that extra folder
## The single HTML file should function perfectly on its own, and can be share or moved
## without losing functionality

## might want to look at https://plotly.com/r/pca-visualization/



