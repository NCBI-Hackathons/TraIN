library(tidyverse)
library(data.table)
library(here)

URL <- "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.5.177/BIOGRID-ALL-3.5.177.tab2.zip"
download.file(URL,"BIOGRID-ALL-3.5.177.tab2.zip")

Biogrid <- unzip(here("BIOGRID-ALL-3.5.177.tab2.zip"))

Biogrid_table <- read.delim(Biogrid)

## Filter for human
Biogrid_table <- Biogrid_table %>%
  filter(Organism.Interactor.A==9606,Organism.Interactor.B==9606)

human_protein_coding_genes <- read.delim(here("human_protein_coding_genes.tsv"))
human_protein_coding_genes$GeneID <- as_factor(human_protein_coding_genes$GeneID)

## Rename the protein coding genes reference
human_prot_gene_A <- dplyr:::select(human_protein_coding_genes,Entrez.Gene.Interactor.A=GeneID,
                                    ENSEMBL_A=Ensembl_gene_identifier)

human_prot_gene_B <- dplyr:::select(human_protein_coding_genes,Entrez.Gene.Interactor.B=GeneID,
                                ENSEMBL_B=Ensembl_gene_identifier) 

# Join tables 
Biogrid_table <- Biogrid_table %>%
  left_join(.,human_prot_gene_A)
Biogrid_table <- Biogrid_table %>%
  left_join(.,human_prot_gene_B)

## Download surface genes
surface_genes <- read.delim(here("human_surface_genes.tsv"))

#Filter for surface genes
Biogrid_table <- Biogrid_table %>%
  filter( ENSEMBL_A %in% surface_genes$Ensembl_gene_identifier &
            ENSEMBL_B %in% surface_genes$Ensembl_gene_identifier)

# Write to file
write_csv(Biogrid_table,here("human_surface_pp_interaction.csv"))
