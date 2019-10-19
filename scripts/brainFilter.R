# set your working directory as SOURCE FILE location (Session -> Set Working Directory -> Source File Location)

brainfile <- list.files(path="../braindata/", pattern = "_counts.csv"); brainfile
brainreg <- gsub(brainfile,pattern = "_counts.csv",replacement = ""); brainreg
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
    brainmean <- data.matrix(rowMeans(filterlist))
  }
  else {
    newdf <- data.matrix(rowMeans(filterlist))
    brainmean <- as.data.frame(cbind(brainmean, newdf))
  }
  
  if (i==length(brainfile)){
    colnames(brainmean) <- brainreg
  }
}

brainmean$max <- apply(brainmean[,], 1, max)

write.csv(brainmean, file = "../output/brain_mean_counts", quote = FALSE)
