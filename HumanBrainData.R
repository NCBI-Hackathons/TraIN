############################################
# Title: Aquisition of normalized RNAseq count data using recount2 - Immunological
# Author: Miranda Darby  (darby@hood.edu)
# Date: created on 5/9/19
# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# SRA entry for downloaded dataset: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP051688
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1
###########################################

### Install BiocManager. Only needs to be done once.
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

### Install recount2. Also only needs to be done once
# BiocManager::install("recount")

library(recount)

# to add brain specific metadata to Recount2 
brainStudies <- add_metadata()
## Subset out studies of interest
#### Find studies that include samples for Nucelus Accumbens, Putamen, Substantia Nigra, Hippocampus, Frontal Cortex, Caudate
idxT1 <- which(brainStudies$tissue_site_1=="Nucelus accumbens") | which(brainStudies$tissue_site_1=="Putamen") | which(brainStudies$tissue_site_1=="Substantia nigra") | which(brainStudies$tissue_site_1=="Hippocampus") | which(brainStudies$tissue_site_1=="Frontal Cortex")| which(brainStudies$tissue_site_1=="Caudate")

#### Find studies that include samples for dorsolateral prefrontal cortex, prefontal cortex
idxT2 <- which(brainStudies$tissue_site_3=="Dorsolateral prefrontal cortex") | which(brainStudies$tissue_site_3=="Prefontal cortex") 
#### Print info on workspace environment
sessionInfo()


