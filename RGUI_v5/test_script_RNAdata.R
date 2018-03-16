# 03/13/2017 Zhi Huang
library(shiny)
library(rsconnect)
library(plyr)
library(data.table)
library(genefilter)
library(Biobase)
library(rPython)
library(WGCNA)
library(GEOquery)
library(dplyr)
library(enrichR)
library(DT)
library(reticulate)
use_python("/Users/zhi/anaconda2/bin/python")
setwd('/Users/zhi/Desktop/GeneCoexpression/RGUI_v5/')

data <- NULL
GEO <- NULL
finalExp <- NULL
finalSym <- NULL
finalSymChar <- NULL
text <- NULL
geneCharVector_global <- NULL
eigengene_matrix <- NULL

#   +------------------------------------------------------------+
#   | 
#   | 
#   |                         Load Data
#   | 
#   |
#   +--------------------------------

data_in = as.matrix(readLines('../matlab_old/RNAdata.csv'), sep = '\n')
data_temp = strsplit(data_in, split=",")
max.length <- max(sapply(data_temp, length))
data_temp <- lapply(data_temp, function(v) { c(v, rep(NA, max.length-length(v)))})
data_temp <- data.frame(do.call(rbind, data_temp))
if(data_temp[dim(data_temp)[1],1] == "!series_matrix_table_end"){
  print("remove last row with \"!series_matrix_table_end\" ")
  data_temp = data_temp[-dim(data_temp)[1],]
}
data <- data_temp[-1,]
print("CSV / txt file Processed.")
print(sprintf("Number of Genes: %d",dim(data)[1]))
print(sprintf("Number of Samples: %d",(dim(data)[2]-1)))

#   +------------------------------------------------------------+
#   | 
#   | 
#   |                      Cleaning the Data
#   | 
#   |
#   +--------------------------------

source("utils.R")
RNA <- as.matrix(data[, 2:dim(data)[2]])
class(RNA) <- "numeric"
geneID <- data.frame(data[, 1])
print(dim(RNA))
print(dim(geneID))

# convert na to 0
RNA[is.na(RNA)] <- 0
# Remove data with lowest 20% absolute exp value shared by all samples
percentile <- 20/100.
# save(RNA, file="/Users/zhi/Desktop/RNA.Rdata")
if (percentile > 0){
  RNA_filtered1 = RNA[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
  geneID_filtered1 = geneID[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
} else {
  RNA_filtered1 = RNA
  geneID_filtered1 = geneID
}
print(dim(RNA_filtered1))
# Remove data with lowest 10% variance across samples
percentile <- 9.99/100.
if (percentile > 0){
  index <- varFilter2(eset = RNA_filtered1, var.cutoff = percentile)
  RNA_filtered2 = RNA_filtered1[index, ]
  geneID_filtered2 = geneID_filtered1[index]
} else {
  RNA_filtered2 = RNA_filtered1
  geneID_filtered2 = geneID_filtered1
}
print(dim(RNA_filtered2))

uniGene <- geneID_filtered2
tmpExp <- RNA_filtered2

nSample <- ncol(tmpExp)
res <- sort.int(rowMeans(tmpExp), decreasing = TRUE, index.return=TRUE)
sortMean <- res$x
sortInd <- res$ix
topN <- min(2000, nrow(tmpExp))
finalExp <- tmpExp[sortInd[1:topN], ]
finalSym <- uniGene[sortInd[1:topN]]
finalSymChar <- as.character(finalSym)

#   +------------------------------------------------------------+
#   | 
#   | 
#   |                         l m Q C M
#   | 
#   |
#   +--------------------------------

# cor_mat = cor(t(finalExp))


step1 = 1
gamma = 0.55
t = 1
lambda = 1
beta = 0.4
minClusterSize = 10
massiveCC = 'Pearson'

python.load("main.py")
mergedCluster <- python.call("mainroutine", step1, as.vector(finalExp), nrow(finalExp), ncol(finalExp), gamma, t, lambda, beta, minClusterSize, massiveCC)
geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
temp_eigengene <- matrix(0, nrow = length(mergedCluster), ncol = dim(finalExp)[2]) # Clusters * Samples

temptext <- ""
for (i in 1:(length(mergedCluster))) {
  vector <- as.matrix(mergedCluster[[i]])
  vector <- vector + 1 # covert python indexing to R indexing
  geneID <- vector
  print(i)
  print(vector)
  # ===== Calculate Eigengene Start
  X <- finalExp[geneID,]
  mu <- rowMeans(X)
  stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
  #normalize X:
  XNorm <- sweep(X,1,mu)
  XNorm <- apply(XNorm, 2, function(x) x/stddev)
  SVD <- svd(XNorm, LINPACK = FALSE)
  temp_eigengene[i,] <- t(SVD$v[,1])
  # ===== Calculate Eigengene Finished
  geneChar <- c(toString(i), finalSymChar[vector])
  geneCharVector[i] <- list(geneChar)
  temptext <- paste(temptext, capture.output(cat(geneChar, sep=',')), sep="\n")
}
temptext <- substring(temptext, 2) # remove first \n separater
geneCharVector_global <<- geneCharVector
text <<- temptext
eigengene_matrix <<- temp_eigengene

## Compute maximum length
max.length <- max(sapply(geneCharVector, length))
## Add NA values to list elements
geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind
geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))