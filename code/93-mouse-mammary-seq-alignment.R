# Goals -----
# The goals for this script is to test out alignment in an R environment
# Most alignment tools run in a linux environment, and they are very computationally intensive.
# Rsubread is the only aligner that can run in R

# This script follows along the following tutorial
# https://combine-australia.github.io/RNAseq-R/07-rnaseq-day2.html
# Which is the basis for the Galaxy tutorial 
# https://galaxyproject.github.io/training-material/topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.html

# If all goes according to plan, this script should produce a BAM file for us to look at with IGV..
# as it turns out: not really, bc IGV needs both a .bam and .bai file to visualize tracks
# also only a few reads are mapping to chr1 , so it's really hard to find one and interpret...


# Data sets ----
# See Figshare https://figshare.com/s/f5d63d8c265a05618137
# which is a subset of data from GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450

# Download dataset from Figshare and place them in `data_raw/MouseMammary`
# note, files in this folder are not tracked by git (due to it's size)
# in total, there are 12 fastq files with 1000 reads each, 4 index files for mm10 chr1,
# and target files with sample information



# Load libraries ----
library(tidyverse)
library(Rsubread)           # from Bioconductor


# Find fastq file path ----
fastq.files <- list.files(path = "data_raw/MouseMammary", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files


# Alignment ----

## Build the index ----
dir_path <- "data_raw/MouseMammary"

## buildindex(basename = "chr1_mm10", reference = "chr1.fa")

## buildindex(basename = file.path("/path/to/index_directory", "chr1_mm10"), reference= "chr1.fa")

## this line of code will take several minutes to run
## it will directly write files to the working directory
## aka the 4 files starting with chr1 already in `data_raw/MouseMammary` 
## maybe come back to this later


## Align reads to chr1 ----
output.files <- file.path(
  dir_path,
  "chr1_alignment",
  paste0(
    sub("\\.fastq\\.gz$", "", basename(fastq.files)),
    ".subread.BAM"
  )
)

align(index = paste0(dir_path, "/chr1_mm10"), 
      readfile1 = fastq.files,
      output_file = output.files)

## The code above aligns each of the 12 samples one after the other
## since there are only 1000 reads/sample, it should be fairly quick.
## However, to run full samples, it would take several hours per sample

## parameters one can change with align()
?Rsubread::align


## proportion of reads that mapped to the reference genome
bam.files <- list.files(path = file.path(dir_path, "chr1_alignment"), 
                        pattern = ".BAM$", full.names = TRUE)
bam.files

props <- propmapped(files=bam.files)
props


# To Be Continued ----
