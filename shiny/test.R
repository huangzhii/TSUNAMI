library(genefilter)
library(Biobase)
library(rPython)

# setwd("/media/zhi/Drive3/GeneCoexpression/matlab_old");
setwd("/Users/zhi/Desktop/GeneCoexpression/shiny"); #mac
# setwd("E:/GeneCoexpression/shiny"); #win
source("utils.R")

data<-read.csv("../matlab_old/RNAdata.csv", header=T, stringsAsFactors=F)

# Step 0
RNA <- as.matrix(data[1:dim(data)[1], 2:dim(data)[2]])
geneID <- data.frame(data[1:dim(data)[1], 1])
# Remove data with lowest 20% absolute exp value shared by all samples
RNA_filtered1 = RNA[apply(RNA,1,max) > quantile(RNA, 0.2)[[1]], ]
geneID_filtered1 = geneID[apply(RNA,1,max) > quantile(RNA, 0.2)[[1]], ]
# Remove data with lowest 10% variance across samples
index <- varFilter2(eset = RNA_filtered1, var.cutoff = 0.1)
RNA_filtered2 = RNA_filtered1[index, ]
geneID_filtered2 = geneID_filtered1[index]
expData <- RNA_filtered2
res <- highExpressionProbes(geneID_filtered2, geneID_filtered2, expData)
ind1 <- res$first
uniGene <- res$second
tmpExp <-expData[ind1,]
nSample <- ncol(tmpExp)
res <- sort.int(rowMeans(tmpExp), decreasing = TRUE, index.return=TRUE)
sortMean <- res$x
sortInd <- res$ix
topN <- min(20000, nrow(tmpExp))
finalExp <- tmpExp[sortInd[1:topN], ]
finalExp[is.nan(finalExp)] <- 0
finalSym <- uniGene[sortInd[1:topN]]
finalSymChar <- as.character(finalSym)

finalExp[is.na(finalExp)] <- 0
# Start the clock!
# ptm <- proc.time()

step1 = 1
gamma = 0.5
t = 1
lambda = 1
beta = 0.4
minClusterSize = 10
python.load("main.py")
# mergedCluster <- python.call("mainroutine", step1, as.vector(finalExp), nrow(finalExp), ncol(finalExp), gamma, t, lambda, beta, minClusterSize)

# cMatrix <- massivePCC_withoutNaN(finalExp)
# Stop the clock
# ptm <- proc.time() - ptm
# print(ptm)

# save(mergedCluster, file = "./mergedCluster.RData")
load(file = "./mergedCluster.RData")
geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))

text <- ""
for (i in 1:(length(mergedCluster))) {
  vector <- as.matrix(mergedCluster[[i]])
  vector <- vector + 1 # covert python indexing to R indexing
  geneChar <- finalSymChar[vector]
  geneCharVector[i] <- list(geneChar)
  text <- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
}

