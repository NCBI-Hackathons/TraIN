
## import data
immcount <- read.csv("./immunoGeneCounts.csv", header=T) 
immcol_raw <- read.csv("./immunoColData.csv", header=F)
immcol <- read.csv("./immunoColData.csv", header=F)

# first row is the "run" column, therefore remove
immcol <- immcol[,-1] 

# the first row will be the header
colnames(immcol) <- as.character(unlist(immcol[1,])) 

# add in column name
immcol = immcol[-1, ]
names(immcol)[22]<-"cell_type" 
names(immcol)[23]<-"assay" 
names(immcol)[24]<-"donor" 

# check unique sample id for each row using run column
length(unique(immcol$run))
nrow(immcol)

immcol$cell_type <- substr(immcol$cell_type, 13, length(immcol$cell_type)-1)

## From immcol data: see how many sample id in each cell type
library(dplyr)

immcol %>% 
  group_by(cell_type) %>% 
  summarise(count = n())

### create list of sample id and cell type
id_cell_lookup <- immcol[,c("cell_type","run")]

## From immcount data: output sample id into a variable
sample_id <- data.frame(colnames(immcount)[-1])
names(sample_id)[1]<-"immcount_sample_id" 

## merge two data
merge_id_ctype <- merge(id_cell_lookup,sample_id, by.x = "run", by.y = "immcount_sample_id", all=T)

## transform data
immcount_t <- data.frame(t(immcount))
colnames(immcount_t) <- as.character(unlist(immcount_t[1,])) 
immcount_t <- immcount_t[-1,] 
immcount_t$run <- rownames(immcount_t)
rownames(immcount_t) <- c()
names(immcount_t)[58038]<-"run" 
 
### pull in cell type in immcount 
immcount_t$cell_type<-with(merge_id_ctype, cell_type[match(immcount_t$run, run)])
x=data.frame(immcount_t[1:10,58038:58039])  # check
x

# output data
#write.csv(immcount_t,"ImmunoGeneCounts_w_cell_type.csv", row.names = F)

# test
cell_b <- immcount_t[immcount_t$cell_type=="primary human B cells",]




