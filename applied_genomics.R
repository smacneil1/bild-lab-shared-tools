# Demo for Applied Genomics Class

# Read in the file 


setwd("~/Documents/")


CCLE_Data=read.table("bild_signatures/CCLE_data/CCLE_Breast_RNAseq_TPMlog.txt",header = T,sep='\t', row.names=1)
View(CCLE)

#find the dimentions of the file
dim(CCLE_Data)

#just look at the top
head(CCLE_Data)

# look at the columns
head(CCLE_Data[,1])
colnames(CCLE_Data)

#look at the first row
CCLE_Data[1,]
rownames(CCLE_Data)

#look at both
dimnames(CCLE_Data)

#download the packages
source("http://www.Bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite( c("ShortRead", "edgeR") )

#seperate the data
FirstHalf=CCLE_Data[1:50,1:28]

View(FirstHalf)
SecondHalf=CCLE_Data[1:1,56]

heatmap.2(t(as.matrix(CCLE_Data)), col = bluered, trace='none',margins = c(9, 9), cexCol=1)


heatmap.2(t(as.matrix(CCLE_Data[1:20,])), col = bluered, trace='none',margins = c(9, 9), cexCol=1)

