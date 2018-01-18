library(genefilter)
setwd("/media/zhi/Drive3/GeneCoexpression/matlab_old");
data<-read.csv("./RNAdata.csv", header=T, stringsAsFactors=F)


# Step 1
RNA <- data.frame(data[1:dim(data)[1], 2:dim(data)[2]])

geneID <- data.frame(data[1:dim(data)[1], 1])
