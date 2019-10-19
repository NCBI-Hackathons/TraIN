# Recount2 paper: https://www.biorxiv.org/content/10.1101/247346v2
# Recount2 published workflow: https://f1000research.com/articles/6-1558/v1
# Recount2 website: https://jhubiostatistics.shinyapps.io/recount/
# Recount-brain paper: https://www.biorxiv.org/content/10.1101/618025v1

## Install Packages
install.packages("BiocManager")
BiocManager::install('tidyverse')
BiocManager::install('data.table')
BiocManager::install("recount")
library(BiocManager)
library(recount)
library(data.table)
library(tidyverse)

### The data can be found on the GTEX page at https://jhubiostatistics.shinyapps.io/recount/
### Under the data, click brain, and download the V2 rse data.
load("/home/team9_hs19/data/raw/gtex.data.Rdata")

#get the counts and filter for rna quality.
#remove the cerebellar hemisphere, spinal cord and the cortex samples. 
#we are not interested in the spinal cord. The cortex and cerebellar hemisphere are sample duplicates
#https://gtexportal.org/home/faq

exclude_regions <- c("Brain - Cerebellar Hemisphere","Brain - Spinal cord (cervical c-1)","Brain - Cortex")

metadata <- colData(gtex_brain) %>% as.data.frame() %>% rownames_to_column(var='accession') %>%
  filter(!(smtsd %in% exclude_regions)) %>% filter(smrin >= 7.5)
counts <- assays(rse_gene)$counts
counts <- counts[,colnames(counts) %in% metadata$accession]

# write out the files independently..
amygdala_metadata <- metadata %>% filter(smtsd=="Brain - Amygdala")
amygdala_counts <- counts[,colnames(counts) %in% amygdala_metadata$accession]
write.csv(amygdala_counts, file="amygdala_counts.csv",quote=FALSE)
write.csv(amygdala_metadata, file="amygdala_sampleData.csv",quote=FALSE)

BA24_metadata <- metadata %>% filter(smtsd=="Brain - Anterior cingulate cortex (BA24)")
BA24_counts <- counts[,colnames(counts) %in% BA24_metadata$accession]
write.csv(BA24_counts, file="BA24_counts.csv",quote=FALSE)
write.csv(BA24_metadata, file="BA24_sampleData.csv",quote=FALSE)

caudate_metadata <- metadata %>% filter(smtsd == "Brain - Caudate (basal ganglia)")
caudate_counts <- counts[,colnames(counts) %in% caudate_metadata$accession]
write.csv(caudate_counts, file="caudate_counts.csv",quote=FALSE)
write.csv(caudate_metadata, file="caudate_sampleData.csv",quote=FALSE)

cerebellum_metadata <- metadata %>% filter(smtsd == "Brain - Cerebellum")
cerebellum_counts <- counts[,colnames(counts) %in% cerebellum_metadata$accession]
write.csv(cerebellum_counts, file="cerebellum_counts.csv", quote=FALSE)
write.csv(cerebellum_metadata, file="cerebellum_sampleData.csv",quote=FALSE)

BA9_metadata <- metadata %>% filter(smtsd == "Brain - Frontal Cortex (BA9)")
BA9_counts <- counts[,colnames(counts) %in% BA9_metadata$accession]
write.csv(BA9_counts, file="BA9_counts.csv",quote=FALSE)
write.csv(BA9_sampleData, file="BA9_sampleData.csv",quote=FALSE)

hippocampus_metadata <- metadata %>% filter(smtsd=="Brain - Hippocampus")
hippocampus_counts <- counts[,colnames(counts) %in% hippocampus_metadata$accession]
write.csv(hippocampus_counts, file="hippocampus_counts.csv",quote=FALSE)
write.csv(hippocampus_metadata, file="hippocampus_sampleData.csv",quote=FALSE)

hypothalamus_metadata <- metadata %>% filter(smtsd == "Brain - Hypothalamus")
hypothalamus_counts <- counts[,colnames(counts) %in% hypothalamus_metadata$accession]
write.csv(hypothalamus_counts, file="hypothalamus_counts.csv",quote=FALSE)
write.csv(hypothalamus_sampleData, file="hypothalamus_sampleData.csv",quote=FALSE)

nucleus_accumbens_metadata <- metadata %>% filter(smtsd == "Brain - Nucleus accumbens (basal ganglia)")
nucleus_accumbens_counts <-  counts[,colnames(counts) %in% nucleus_accumbens_metadata$accession]
write.csv(nucleus_accumbens_counts, file="nucleus_accumbens_counts.csv",quote=FALSE)
write.csv(nucleus_accumbens_sampleData, file="nucleus_accumbens_sampleData.csv",quote=FALSE)

putamen_metadata <- metadata %>% filter(smtsd == "Brain - Putamen (basal ganglia)")
putamen_counts <- counts[,colnames(counts) %in% putamen_metadata$accession]
write.csv(putamen_counts, file="putamen_counts.csv",quote=FALSE)
write.csv(putamen_sampleData, file="putamen_sampleData.csv",quote=FALSE)

substantia_nigra_metadata <- metadata %>% filter(smtsd == "Brain - Substantia nigra")
substantia_nigra_counts <- counts[,colnames(counts) %in% substantia_nigra_metadata$accession]
write.csv(substantia_nigra_counts, file="substantia_nigra_counts.csv",quote=FALSE)
write.csv(substantia_nigra_sampleData, file="substantia_nigra_sampleData.csv",quote=FALSE)
