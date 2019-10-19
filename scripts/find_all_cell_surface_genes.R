library(tidyverse)
download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", destfile = "gene2go.gz")

gene2go <- read.delim("~/Library/CloudStorage/iCloudDrive/Desktop/all_human_genes/gene2go.gz", header=TRUE)
human_protein_coding_genes <- read.delim("~/Library/CloudStorage/iCloudDrive/Desktop/all_human_genes/human_protein_coding_genes.tsv")

human_surface_genes <- gene2go %>%
  filter(X.tax_id == 9606) %>%
  ## find the genes with GO_ID 9986 and all the other children GO_ID
  filter(GO_ID == "GO:0009986" | GO_ID == "GO:0009897" |GO_ID == "GO:0031362" | GO_ID == "GO:0031232" | GO_ID == 
          "GO:0098591" | GO_ID == "GO:0031233" | GO_ID == "GO:0071575") %>%
  dplyr::select(GeneID, GO_ID, Evidence, GO_term, PubMed) %>%
  left_join(human_protein_coding_genes, by="GeneID")

write.table(human_surface_genes, file="/Users/yichengtsai/Desktop/all_human_genes/human_surface_genes.tsv", sep='\t', row.names=FALSE, quote=FALSE)
