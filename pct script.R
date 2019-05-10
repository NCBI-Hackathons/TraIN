###########################################################################################################
### Input file: immunoGeneCounts.csv and immunoColData.csv from github @ 1:56 PM
### Output file: Expression percentage by cell type for each gene.csv
### Metric: 1) extract sample id and cell type from immunoColData.csv 
###                -> 7 cell types identified, each includes 8 sample id. All sample id are unique.
###         2) a. add cell type info based on sample id infomation in immunoGeneCounts.csv. 
###            b. immunoGeneCounts.csv has 56 sample id, all matches id in immunoColData.csv.
###            c. if count in immunoGeneCounts.csv>0 then recode as 1, otherwise 0.
###            d. calcualte expression percentage which using sum(# of 1)/9 for each cell type.
###               -> output file: for each row is cell type, each column is gene, cell is the percentage
###                 of expression.
### Creater: J.Gu
### Date: May 9, 2019
###########################################################################################################


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

# subset to cell_b
#cell_b <- immcount_t[immcount_t$cell_type=="primary human B cells",]
# not subset
cell_b <- immcount_t

# change factor to numeric for gene columns
cell_b_numeric <- data.frame(lapply(cell_b[,1:(ncol(cell_b)-2)], function(x) as.numeric(as.character(x))))

# if larger than 0 then assign to 1
cell_b_numeric_1 <- ifelse(cell_b_numeric==0, 0,1)

# combine numeric and check data types
cell_b_numeric_2 <- cbind(cell_b_numeric_1,cell_b[,(ncol(cell_b)-1):ncol(cell_b)])
class(cell_b_numeric_2[,(ncol(cell_b)-2)]); class(cell_b_numeric_2[,(ncol(cell_b)-1)]);class(cell_b_numeric_2[,ncol(cell_b)])

# create summary for total
cell_b_numeric_3 <- cell_b_numeric_2 %>% group_by(cell_type) %>% summarise_if(is.numeric, sum)

# since each cell has 8 samples, therefore divided by 8
cell_b_numeric_4 <- cbind(cell_b_numeric_3[,1], cell_b_numeric_3[,2:ncol(cell_b_numeric_3)]/8)

# output result 
write.csv(cell_b_numeric_4, "Expression percentage by cell type for each gene.csv", row.names=F)

# END

