############################################
# Title: Aquisition of normalized RNAseq count data using recount2 - Immunological
# Author: Miranda Darby  (darby@hood.edu)
# Date: created on 5/9/19
# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# Recount-brain paper for future use: https://www.biorxiv.org/content/10.1101/618025v1
###########################################

### Install BiocManager. Only needs to be done once.
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

### Install recount2. Also only needs to be done once
# BiocManager::install("recount")

library(recount)

unique(brainStudies$tissue_site_1)
# [1] "Cerebral cortex"    NA                   "Mixed"              "Whole brain"       
# [5] "Cerebellum"         "Brainstem"          "Corpus callosum"    "Nucleus accumbens" 
# [9] "Dura mater"         "Putamen"            "Substantia nigra"   "Hippocampus"       
# [13] "Lumbar spinal cord" "Frontal Cortex"     "Caudate" 

# to add brain specific metadata to Recount2 
brainStudies <- add_metadata()
## Subset out studies of interest

#### Find studies that include samples for Nucelus Accumbens, Putamen, Substantia Nigra, Hippocampus, Frontal Cortex, Caudate
NucAcc <- which(brainStudies[, "tissue_site_1"] == "Nucleus accumbens")
Put <- which(brainStudies$tissue_site_1=="Putamen")
SubNig <- which(brainStudies$tissue_site_1=="Substantia nigra")
Hip <- which(brainStudies$tissue_site_1=="Hippocampus")
FC <- which(brainStudies$tissue_site_1=="Frontal Cortex")
Caud <- which(brainStudies$tissue_site_1=="Caudate")

#### Find studies that include samples for dorsolateral prefrontal cortex, prefontal cortex
unique(brainStudies$tissue_site_3)
# [1] "Superior temporal gyrus"       
# [2] "Dorsolateral prefrontal cortex"
# [3] NA                              
# [4] "Motor cortex"                  
# [5] "Superior frontal gyrus"        
# [6] "Anterior cingulate gyrus"      
# [7] "Somatosensory cortex"          
# [8] "Prefrontal cortex"             
# [9] "Anterior prefrontal cortex"    
DPC <-which(brainStudies$tissue_site_3=="Dorsolateral prefrontal cortex")
PC <- which(brainStudies$tissue_site_3=="Prefrontal cortex")

AllRegions <- c(NucAcc, Put, SubNig, Hip, FC, Caud, DPC, PC)
(AllRegions <- unique(AllRegions))
(BestRegions <- c(NucAcc, SubNig, Hip, FC))

#### Print info on workspace environment
sessionInfo()


