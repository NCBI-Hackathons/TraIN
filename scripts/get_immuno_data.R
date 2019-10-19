# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# SRA entry for downloaded dataset: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP051688
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1

### Install BiocManager. Only needs to be done once.
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("recount")

library(recount)
library(tidyverse)
library(data.table)
## Download the data if it is not there, check if it has been downloaded
if(!file.exists(file.path("SRP051688", "rse_gene.Rdata"))) {
  download_study("SRP051688", type = "rse-gene")
}
file.exists(file.path("SRP051688", "rse_gene.Rdata"))

## Load the data
load(file.path("SRP051688", "rse_gene.Rdata"))

#get metadata
metadata <- colData(rse_gene) %>% as.data.frame() %>% rownames_to_column(var='accession')
counts <- assay(rse_gene)

fwrite(metadata,'immunoColData.csv',quote=FALSE)
fwrite(counts,'immunoGeneCounts.csv',quote=FALSE)