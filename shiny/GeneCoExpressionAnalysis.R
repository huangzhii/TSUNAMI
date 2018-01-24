# 01/17/2018 Zhi Huang
library(genefilter)
library(Biobase)
library(rPython)

lmQCM <- function(step1, data, gamma, lambda, t, beta, minClusterSize) {
  # print(sprintf("gamma: %.2f",gamma))
  # print(sprintf("lambda: %.2f",lambda))
  # print(sprintf("t: %.2f",t))
  # print(sprintf("beta: %.2f",beta))
  # print(sprintf("minClusterSize: %d",minClusterSize))
  
  # setwd("/media/zhi/Drive3/GeneCoexpression/matlab_old");
  # setwd("/Users/zhi/Desktop/GeneCoexpression/shiny"); #mac
  # setwd("E:/GeneCoexpression/shiny"); #win
  source("utils.R")
  
  # Step 0
  print("Preprocessing ...")
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
  finalExp[is.na(finalExp)] <- 0
  finalSymChar <- as.character(finalSym)
  print("Preprocessing Finished.")
  
  # Start the clock!
  # ptm <- proc.time()
  
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
    # geneCharVector[i] <- capture.output(cat(geneChar, sep=' '))
    # text <- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
  }
  data <- do.call(rbind, lapply(geneCharVector, 
                         function(x) paste(x,collapse=" ")))
  return(data)
}
