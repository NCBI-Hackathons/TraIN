library(tidyverse)

download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", destfile = "gene2ensembl.gz")

gene2ensembl <- read.delim("~/Library/CloudStorage/iCloudDrive/Desktop/all_human_genes/gene2ensembl.gz", header=TRUE)

head(gene2ensembl %>% filter(X.tax_id == 9606))

human_protein_coding_gene <- gene2ensembl %>%
  filter(X.tax_id == 9606) %>%
  filter(protein_accession.version != "-") %>%
  select(GeneID, Ensembl_gene_identifier) %>%
  unique()

write.table(human_protein_coding_gene, file="/Users/yichengtsai/Desktop/all_human_genes/human_protein_coding_genes.tsv", sep='\t', row.names=FALSE, quote=FALSE)