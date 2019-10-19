library(tidyverse)

download.file("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", destfile = "gene2ensembl.gz")

gene2ensembl <- read.delim("gene2ensembl.gz", header=TRUE)

human_protein_coding_gene <- gene2ensembl %>%
  filter(X.tax_id == 9606) %>%
  filter(protein_accession.version != "-") %>%
  select(GeneID, Ensembl_gene_identifier) %>%
  unique()

write.table(human_protein_coding_gene, file="human_protein_coding_genes.tsv", sep='\t', row.names=FALSE, quote=FALSE)
