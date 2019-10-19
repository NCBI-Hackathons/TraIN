## Set your working directory as SOURCE FILE location (Session -> Set Working Directory -> Source File Location)
brainfile <- list.files(path="../braindata/", pattern = "_counts.csv", full.names = TRUE); brainfile
brainreg <- gsub(list.files(path="../braindata/", pattern = "_counts.csv"), pattern = "_counts.csv",replacement = ""); brainreg

## FILTER BY PROTEIN CODING GENES
pcgenes <- data.frame(read.csv(file = "../reference/human_protein_coding_genes.tsv", header = TRUE, sep = "\t"))

for (i in 1:length(brainfile)){
  name <- brainreg[i]
  file <- brainfile[i]
  genelist <- read.csv(file, header = TRUE, row.names = 1) 

  # create filtered list of genes (protein coding only)
  (index <- which(rownames(genelist) %in% pcgenes$Ensembl_gene_identifier))
  filterlist <- genelist[index,]
  assign(name, filterlist)

  # create table with mean values of gene expression (in cpm) for each region
  if (i==1){
    brainmean <- data.frame(rowMeans(filterlist))
  }
  else {
    newdf <- data.frame(rowMeans(filterlist))
    brainmean <- as.data.frame(cbind(brainmean, newdf))
  }
  
  if (i==length(brainfile)){
    colnames(brainmean) <- brainreg
  }
}

brainmean$max <- apply(brainmean[,], 1, max)

#write.csv(brainmean, file = "../outputs/brain_mean_counts.csv", quote = FALSE)

## CONVERT SURFACE MOLECULES FROM UNIPROT TO ENSEMBL
library(biomaRt)

surflist <- as.matrix(read.csv(file = "../reference/Hs_cell-surface-671.csv", header = TRUE)[,1])
uniprotid <- surflist[grep("UniProtKB:",surflist)]
uniprotid <- (gsub("UniProtKB:", "",uniprotid))

mart.hs <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneid <- getBM(attributes=c("uniprot_gn_id", "ensembl_gene_id", "external_gene_name"), filter="uniprot_gn_id", values = uniprotid, mart=mart.hs)

## FILTER BRAIN DATA FOR SURFACE MOLECULES
brainsurf <- brainmean[rownames(brainmean) %in% geneid$ensembl_gene_id,]

# LOAD PPI PAIRS
ppilist <- read.csv(file = "../reference/human-PPI-pairs.csv", header = TRUE)
ppicombA <- matrix(ppilist[, "ENSEMBL_A"])
ppicombB <- matrix(ppilist[, "ENSEMBL_B"])
ppicomb <- rbind(ppicombA,ppicombB); ppicomb

#FILTER DATA FOR GENES LIST
brainppi <- brainsurf[0,]
for (i in 1:nrow(ppilist)) {
  if((ppicombA[i] %in% rownames(brainsurf))&(ppicombB[i] %in% rownames(brainsurf))) {
    indexa <- which(rownames(brainsurf) %in% ppicombA[i])
    print(indexa)
    indexb <- which(rownames(brainsurf) %in% ppicombB[i])
    print(indexb)
    brainppi <- rbind(brainppi, brainsurf[indexa,])
    brainppi <- rbind(brainppi, brainsurf[indexb,])
  }
}