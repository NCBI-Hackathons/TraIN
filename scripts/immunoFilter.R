## Set your working directory as SOURCE FILE location (Session -> Set Working Directory -> Source File Location)
immunofile <- list.files(path="../immunodata/", pattern = "_TPM.csv", full.names = TRUE); immunofile
immunoreg <- gsub(list.files(path="../immunodata/", pattern = "_TPM.csv"), pattern = "_TPM.csv",replacement = ""); immunoreg

## FILTER BY PROTEIN CODING GENES
pcgenes <- data.frame(read.csv(file = "../reference/human_protein_coding_genes.tsv", header = TRUE, sep = "\t"))

for (i in 1:length(immunofile)){
  name <- immunoreg[i]
  file <- immunofile[i]
  genelist <- read.csv(file, header = TRUE, row.names = 1)
  genelist <- genelist[3:ncol(genelist)]
  #assign(name, genelist)
  
  # create table with mean values of gene expression (in cpm) for each region
  if (i==1){
    immunomean <- data.frame(rowMeans(genelist))
  }
  else {
    newdf <- data.frame(rowMeans(genelist))
    immunomean <- as.data.frame(cbind(immunomean, newdf))
  }
  
  if (i==length(immunofile)){
    colnames(immunomean) <- immunoreg
  }
}

# create filtered list of genes (protein coding only)
library(dplyr)
immunomean$ensembl_id <- gsub("\\..*","",rownames(immunomean))
immunomean <- distinct(.data = immunomean, ensembl_id,.keep_all = TRUE)
dim(immunomean)

rownames(immunomean) <- immunomean$ensembl_id
immunomean <- immunomean[,1:ncol(immunomean)-1]

(index <- which(rownames(immunomean) %in% pcgenes$Ensembl_gene_identifier))
immunomean <- immunomean[index,]

immunomean$max <- apply(immunomean[,], 1, max)

write.csv(immunomean, file = "../outputs/immuno_mean_counts.csv", quote = FALSE)

## CONVERT SURFACE MOLECULES FROM UNIPROT TO ENSEMBL
library(biomaRt)

surflist <- as.matrix(read.csv(file = "../reference/Hs_cell-surface-671.csv", header = TRUE)[,1])
uniprotid <- surflist[grep("UniProtKB:",surflist)]
uniprotid <- (gsub("UniProtKB:", "",uniprotid))

mart.hs <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneid <- getBM(attributes=c("uniprot_gn_id", "ensembl_gene_id_version", "ensembl_gene_id", "external_gene_name"), filter="uniprot_gn_id", values = uniprotid, mart=mart.hs)

## FILTER immuno DATA FOR SURFACE MOLECULES
immunosurf <- immunomean[rownames(immunomean) %in% geneid$ensembl_gene_id,]

write.csv(immunosurf, file = "../outputs/immuno_surf_counts.csv", quote = FALSE)

# LOAD PPI PAIRS
ppilist <- read.csv(file = "../reference/biogrid_table_pp_interaction.csv", header = TRUE)
ppicombA <- matrix(ppilist[, "ENSEMBL_A"])
ppicombB <- matrix(ppilist[, "ENSEMBL_B"])
ppicomb <- rbind(ppicombA,ppicombB); ppicomb

#FILTER DATA FOR PPI PAIRS AND GENERATE TWO MATRICES (EACH GENE IN PAIR)
immunoppiA <- immunosurf[0,]
immunoppiB <- immunosurf[0,]

immunosurf$ensembl_id <- rownames(immunosurf)
rownames(immunosurf) <- NULL
for (i in 1:nrow(ppilist)) {
  if((ppicombA[i] %in% immunosurf$ensembl_id)&(ppicombB[i] %in% immunosurf$ensembl_id)) {
    indexa <- which(immunosurf$ensembl_id %in% ppicombA[i])
    indexb <- which(immunosurf$ensembl_id %in% ppicombB[i])
    
    immunoppiA <- rbind(immunoppiA, immunosurf[indexa,])
    immunoppiB <- rbind(immunoppiB, immunosurf[indexb,])
  }
}

immunoppiA <- immunoppiA[ ,c(ncol(immunoppiA), 1:ncol(immunoppiA)-1)]
immunoppiB <- immunoppiB[ ,c(ncol(immunoppiB), 1:ncol(immunoppiB)-1)]

write.csv(immunoppiA, file = "../outputs/immuno_ppiA_counts.csv", quote = FALSE)
write.csv(immunoppiB, file = "../outputs/immuno_ppiB_counts.csv", quote = FALSE)

# CREATE MERGE MATRIX WITH BOTH GENE NAMES AND THE MINIMUM OF THE MAX REGION MEAN COUNTS
minppi <- data.frame("Amax"=immunoppiA$max, "Bmax"=immunoppiB$max); minppi

for (i in 1:nrow(immunoppiA)){
  if (immunoppiA$ensembl_id[i] %in% geneid$ensembl_gene_id) {
    immunoppiA$genename[i] <- geneid$external_gene_name[which(geneid$ensembl_gene_id %in% immunoppiA$ensembl_id[i])]
  }
  else {
    immunoppiA$genename[i] <- NA
  }
  if (immunoppiB$ensembl_id[i] %in% geneid$ensembl_gene_id) {
    immunoppiB$genename[i] <- geneid$external_gene_name[which(geneid$ensembl_gene_id %in% immunoppiB$ensembl_id[i])]
  }
  else {
    immunoppiB$genename[i] <- NA
  }
}

immunoppi <- data.frame("Ensembl_GeneID_A" = immunoppiA$ensembl_id,
                       "Ensembl_GeneID_B" = immunoppiB$ensembl_id,
                       "GeneA" = immunoppiA$genename,
                       "GeneB" = immunoppiB$genename,
                       "min" = apply(minppi, 1, max))

write.csv(immunoppi, file = "../outputs/immuno_ppi_mapping.csv", quote = FALSE, row.names = FALSE)