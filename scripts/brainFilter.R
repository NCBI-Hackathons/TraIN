## Set your working directory as SOURCE FILE location (Session -> Set Working Directory -> Source File Location)
brainfile <- list.files(path="../braindata/", pattern = "_counts.csv", full.names = TRUE); brainfile
brainreg <- gsub(list.files(path="../braindata/", pattern = "_counts.csv"), pattern = "_counts.csv",replacement = ""); brainreg

## FILTER BY PROTEIN CODING GENES
pcgenes <- data.frame(read.csv(file = "../reference/human_protein_coding_genes.tsv", header = TRUE, sep = "\t"))

for (i in 1:length(brainfile)){
  name <- brainreg[i]
  file <- brainfile[i]
  genelist <- read.csv(file, header = TRUE, row.names = 1)
  #assign(name, genelist)
  
  # create table with mean values of gene expression (in cpm) for each region
  if (i==1){
    brainmean <- data.frame(rowMeans(genelist))
  }
  else {
    newdf <- data.frame(rowMeans(genelist))
    brainmean <- as.data.frame(cbind(brainmean, newdf))
  }
  
  if (i==length(brainfile)){
    colnames(brainmean) <- brainreg
  }
}

# create filtered list of genes (protein coding only)
library(dplyr)
brainmean$ensembl_id <- gsub("\\..*","",rownames(brainmean))
brainmean <- distinct(.data = brainmean, ensembl_id,.keep_all = TRUE)
dim(brainmean)

rownames(brainmean) <- brainmean$ensembl_id
brainmean <- brainmean[,1:ncol(brainmean)-1]

(index <- which(rownames(brainmean) %in% pcgenes$Ensembl_gene_identifier))
brainmean <- brainmean[index,]

brainmean$max <- apply(brainmean[,], 1, max)

write.csv(brainmean, file = "../outputs/brain_mean_counts.csv", quote = FALSE)

## CONVERT SURFACE MOLECULES FROM UNIPROT TO ENSEMBL
library(biomaRt)

surflist <- as.matrix(read.csv(file = "../reference/Hs_cell-surface-671.csv", header = TRUE)[,1])
uniprotid <- surflist[grep("UniProtKB:",surflist)]
uniprotid <- (gsub("UniProtKB:", "",uniprotid))

mart.hs <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneid <- getBM(attributes=c("uniprot_gn_id", "ensembl_gene_id", "external_gene_name"), filter="uniprot_gn_id", values = uniprotid, mart=mart.hs)

# ## FILTER BRAIN DATA FOR SURFACE MOLECULES
# brainsurf <- brainmean[rownames(brainmean) %in% geneid$ensembl_gene_id,]
# 
# write.csv(brainsurf, file = "../outputs/brain_surf_counts.csv", quote = FALSE)
brainsurf <- brainmean
# LOAD PPI PAIRS
ppilist <- read.csv(file = "../reference/human_surface_pp_interaction.csv", header = TRUE)
ppicombA <- matrix(ppilist[, "ENSEMBL_A"])
ppicombB <- matrix(ppilist[, "ENSEMBL_B"])
ppicomb <- rbind(ppicombA,ppicombB); ppicomb

#FILTER DATA FOR PPI PAIRS AND GENERATE TWO MATRICES (EACH GENE IN PAIR)
brainppiA <- brainsurf[0,]
brainppiB <- brainsurf[0,]

brainsurf$ensembl_id <- rownames(brainsurf)
rownames(brainsurf) <- NULL
for (i in 1:nrow(ppilist)) {
  if((ppicombA[i] %in% brainsurf$ensembl_id)&(ppicombB[i] %in% brainsurf$ensembl_id)) {
    indexa <- which(brainsurf$ensembl_id %in% ppicombA[i])
    indexb <- which(brainsurf$ensembl_id %in% ppicombB[i])
    
    brainppiA <- rbind(brainppiA, brainsurf[indexa,])
    brainppiB <- rbind(brainppiB, brainsurf[indexb,])
  }
}

brainppiA <- brainppiA[ ,c(ncol(brainppiA), 1:ncol(brainppiA)-1)]
brainppiB <- brainppiB[ ,c(ncol(brainppiB), 1:ncol(brainppiB)-1)]

write.csv(brainppiA, file = "../outputs/brain_ppiA_counts.csv", quote = FALSE)
write.csv(brainppiB, file = "../outputs/brain_ppiB_counts.csv", quote = FALSE)

# CREATE MERGE MATRIX WITH BOTH GENE NAMES AND THE MINIMUM OF THE MAX REGION MEAN COUNTS
minppi <- data.frame("Amax"=brainppiA$max, "Bmax"=brainppiB$max); minppi

for (i in 1:nrow(brainppiA)){
  if (brainppiA$ensembl_id[i] %in% geneid$ensembl_gene_id) {
    brainppiA$genename[i] <- geneid$external_gene_name[which(geneid$ensembl_gene_id %in% brainppiA$ensembl_id[i])]
  }
  else {
    brainppiA$genename[i] <- NA
  }
  if (brainppiB$ensembl_id[i] %in% geneid$ensembl_gene_id) {
    brainppiB$genename[i] <- geneid$external_gene_name[which(geneid$ensembl_gene_id %in% brainppiB$ensembl_id[i])]
  }
  else {
    brainppiB$genename[i] <- NA
  }
}

brainppi <- data.frame("Ensembl_GeneID_A" = brainppiA$ensembl_id,
                       "Ensembl_GeneID_B" = brainppiB$ensembl_id,
                       "GeneA" = brainppiA$genename,
                       "GeneB" = brainppiB$genename,
                       "min" = apply(minppi, 1, max))

write.csv(brainppi, file = "../outputs/brain_ppi_mapping.csv", quote = FALSE, row.names = FALSE)