# DeCoDe Institute

This repository contains my code, notes, outputs, and exploratory analysis from participating in the NIH Common Fund Data Ecosystem [DeCoDE Institute 2026](https://github.com/CFDETrainingCenter/decode-institute-2026). More information about the program is available from the [CFDE Training Center](https://www.orau.org/cfde-trainingcenter/training/decode-institute.html).


<br/>

## Repository Structure

This repo is **a personal working space** for learning and exploration during the DeCoDE Institute. Some scripts are follow-alongs from the training material, with modifications and notes added as I worked through them. Other scripts are independent explorations of datasets. Note, these scripts were never intended to be a not a polished or reproducible analysis package.


- `code/:` R scripts for follow-along, exploration, and independent analyses.
- `code_markdown/`: R Markdown training material. The main curriculum file is `06_R_Intro.Rmd`
- `data_raw/`: Raw metadata files used in exploration
- `data_processed/`: Processed data files copied from `outputs_csv` and saved for later analyses
- `outputs_csv/`: CSV outputs generated from data wrangling steps
- `graphs/`: Figures generated from MoTrPAC dataset exploration
- `graphs_other/`: Figures generated from additional learning exercises
- `galaxy/`: Galaxy-generated report files saved for reference

<br/>


## Main analysis threads

### R and tidyverse follow-alongs

These scripts follow the DeCoDE R training curriculum, including:

- exploring the MoTrPAC rat endurance training data package
- working with phenotype and sample metadata
- using tidyverse tools for filtering, grouping, joining, pivoting, and plotting
- generating exploratory figures with `ggplot2`

Relevant scripts include:

- `code/01-intro-R-tidyverse.R`
- `code/02-joining-pivot-plotting.R`
- `code/03-omics-bonus.R`


### DESeq2 and Galaxy replication

Some scripts are focused on understanding differential expression workflows in R and comparing them with Galaxy-based analyses.

Topics include:

- preparing count matrices and sample metadata
- running DESeq2 workflows
- generating MA plots, PCA plots, heatmaps, volcano plots, and sample-distance heatmaps

Relevant scripts include:

- `code/04-galaxy-DESeq2-rep.R`
- `code/xx-DESeq2-tutorial.R`



### MoTrPAC dataset exploration

Additional scripts exploring the MoTrPAC dataset, focusing on liver RNA-seq data at the moment.

Topics include:

- TBA

Relevant scripts include:

- `code/11-liver-analysis.R`



### Data sources

This work uses training datasets and public biomedical data resources introduced through the DeCoDE Institute, including:

- Allissa D. A comprehensive view of the transcriptome during development of the mouse cerebral cortex [Dataset]. Zenodo. 2026. [doi:10.5281/zenodo.20531535](https://doi.org/10.5281/zenodo.20531535)
- MoTrPAC Study Group. 2024. Temporal dynamics of the multi-omic response to endurance exercise training
Nature. Volume 629, pages 174–183 (2024). doi: 10.1038/s41586-023-06877-w
- MoTrPAC Study Group. Temporal dynamics of the multi-omic response to endurance exercise training. Nature. 2024 May;629(8010):174-183. doi: [10.1038/s41586-023-06877-w](https://doi.org/10.1038/s41586-023-06877-w).
  - Access through the `MotrpacRatTraining6moData` [package](https://motrpac.github.io/MotrpacRatTraining6moData/), [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA908279&o=acc_s%3Aa), [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242354)
  
<br/>




