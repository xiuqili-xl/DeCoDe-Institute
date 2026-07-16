# Goal ----
# Investigate how the Carpentry counts data was constructed 
# See: https://carpentries-incubator.github.io/bioc-rnaseq/02-setup.html
# Aka try to recreate it based on GEO deposit 
# See: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870


# Load libraries ----
library(tidyverse)
library(here)


# Load Carpentry data ----
carpentry_counts <- read_csv(here("data_other", "GSE96870_counts_cerebellum.csv"))
carpentry_metadata <- read_csv(here("data_other", "GSE96870_coldata_cerebellum.csv"))



# Extract data from GEO ----
library(GEOquery)

# define list of GSM numbers
gsm_list <- carpentry_metadata$sample


# name of directory
gsm_dir <- here("data_other", "GSE96870_GEO-download")


# loop through each GSM and download its supplementary files
for (gsm in gsm_list) {
  tryCatch({
    ## print message to keep track
    message(paste("Downloading supplementary files for:", gsm))
    
    ## fetch suppl file from GEO, and download it to appropriate folder
    getGEOSuppFiles(GEO = gsm,
                    makeDirectory = FALSE,
                    baseDir = gsm_dir)
    
  }, error = function(e) {
    message(paste("Failed to download for", gsm, ":", e$message))
  })
}


# Combine GEO data into a single matrix  ----
gsm_names <- list.files(gsm_dir)

# set up a variable to store counts, for now just with the gene names
sample_combined_count <- carpentry_counts %>% select(gene)

for (i in 1:length(gsm_names)) {
  # store sample GSM number name for use later
  sample_name <- str_sub(gsm_names[i], start = 1, end = 10)
  
  # extract relevant columns with genes and counts
  sample_count <- read_delim(file = here(gsm_dir, gsm_names[i]), delim = "\t", skip = 1,
                             col_select = c(1, 7))
  # rename the columns
  colnames(sample_count) <- c("gene", sample_name) 
  
  # construct counts df
  sample_combined_count <- full_join(sample_combined_count, sample_count, by = "gene")
  
  # print message to keep track
  message(paste0("After adding ", sample_name, ", we have ", nrow(sample_combined_count), " rows"))
}



# Examine combined count data ----
head(sample_combined_count)
head(carpentry_counts)

sample_combined_count == carpentry_counts
all(sample_combined_count == carpentry_counts)
