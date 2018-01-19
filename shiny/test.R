library(genefilter)
library(Biobase)

# setwd("/media/zhi/Drive3/GeneCoexpression/matlab_old");
# setwd("/Users/zhi/Desktop/GeneCoexpression/shiny"); #mac
setwd("E:/GeneCoexpression/shiny"); #win
source("utils.R")
data<-read.csv("../matlab_old/RNAdata.csv", header=T, stringsAsFactors=F)

# Step 1
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
finalSym <- uniGene[sortInd[1:topN]]
# Start the clock!
ptm <- proc.time()

cMatrix <- massivePCC_withoutNaN(finalExp)
# Stop the clock
ptm <- proc.time() - ptm
print(ptm)

save(cMatrix, file = "../../cMatrix.RData")
