library(tidyverse)
library(data.table)
library(here)

URL <- "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.177/BIOGRID-ALL-3.5.177.tab2.zip"
download.file(URL,"BIOGRID-ALL-3.5.177.tab2.zip")

Biogrid <- unzip(here("BIOGRID-ALL-3.5.177.tab2.zip"))

Biogrid_table <- read.delim(Biogrid)

Biogrid_table <- Biogrid_table %>%
  filter(Organism.Interactor.A==9606,Organism.Interactor.B==9606)

view(Biogrid_table)
human_protein_coding_genes <- read.delim(here("human_protein_coding_genes.tsv"))
human_protein_coding_genes$GeneID <- as_factor(human_protein_coding_genes$GeneID)
human_prot_gene_A <- dplyr:::select(human_protein_coding_genes,Entrez.Gene.Interactor.A=GeneID,
                                    Ensembl_gene_identifier_A=Ensembl_gene_identifier)

human_prot_gene_B <- dplyr:::select(human_protein_coding_genes,Entrez.Gene.Interactor.B=GeneID,
                                Ensembl_gene_identifier_B=Ensembl_gene_identifier) 

Biogrid_table <- Biogrid_table %>%
  left_join(.,human_prot_gene_A)
Biogrid_table <- Biogrid_table %>%
  left_join(.,human_prot_gene_B)

write_csv(Biogrid_table,here("Biogrid_table_pp_interaction.csv"))
