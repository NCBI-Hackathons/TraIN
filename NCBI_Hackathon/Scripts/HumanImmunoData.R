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

## Download the data if it is not there
if(!file.exists(file.path("SRP051688", "rse_gene.Rdata"))) {
  download_study("SRP051688", type = "rse-gene")
}

## 2019-05-08 10:11:16 downloading file rse_gene.Rdata to SRP051688

## Check that the file was downloaded
file.exists(file.path("SRP051688", "rse_gene.Rdata"))

## [1] TRUE

## Load the data
load(file.path("SRP051688", "rse_gene.Rdata"))

rse_gene
# class: RangedSummarizedExperiment 
# dim: 58037 56 
# metadata(0):
#   assays(1): counts
# rownames(58037): ENSG00000000003.14
# ENSG00000000005.5 ... ENSG00000283698.1
# ENSG00000283699.1
# rowData names(3): gene_id bp_length symbol
# colnames(56): SRR1740034 SRR1740035 ...
# SRR1740088 SRR1740089
# colData names(21): project sample ... title
# characteristics