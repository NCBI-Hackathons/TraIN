col.data.path <- "/home/team9_hs19/immunoData/immunoColData.csv"
gene.counts.path <- "/home/team9_hs19/immunoData/immunoGeneCounts.csv"

# reading the annotations file
annotations <- read.csv(file=col.data.path, skip=1, header=FALSE, sep=",")
srr.to.celltype <- annotations[, c(1, 23)]

# indexing based on the rownames
rownames(srr.to.celltype) <- annotations[, 1]

# reading and indexing the gene counts
gene.counts <- read.csv(file=gene.counts.path, header=TRUE, sep=",", row.names=1)
cell.type.names <- make.names(srr.to.celltype[colnames(gene.counts), ]$V23)
colnames(gene.counts) <- cell.type.names

# aggregate the counts based on the cell type
aggregate <- t(rowsum(t(gene.counts), group=cell.type.names))
freqs <- table(cell.type.names)

# Extracting mean of the counts
mean <- data.frame(matrix(mapply('/', aggregate, freqs),  ncol=length(freqs)))
rownames(mean) <- rownames(aggregate)
colnames(mean) <- colnames(aggregate)
mean$colMax <- apply(mean, 2, function(x) max(x) )

write.csv(mean, file="immuno_mean_counts.csv",quote = FALSE)
