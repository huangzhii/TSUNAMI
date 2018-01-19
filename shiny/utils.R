# 01/18/2018 Zhi Huang
library(genefilter)
library(nnet)
library(parallel)
varFilter2 <- function (eset, var.func = IQR, var.cutoff = 0.5, filterByQuantile = TRUE) {
  if (deparse(substitute(var.func)) == "IQR") {
    vars <- genefilter:::rowIQRs(eset)
  }
  else {
    vars <- apply(exprs(eset), 1, var.func)
  }
  if (filterByQuantile) {
    if (0 < var.cutoff && var.cutoff < 1) {
      quant = quantile(vars, probs = var.cutoff)
      selected = !is.na(vars) & vars > quant
    }
    else stop("Cutoff Quantile has to be between 0 and 1.")
  }
  else {
    selected <- !is.na(vars) & vars > var.cutoff
  }
  return(selected)
}

highExpressionProbes <- function (genes, probes, eset){
  meanData <- rowMeans(eset)
  res <- sort.int(genes, index.return=TRUE)
  sGenes <- res$x
  sInd <- res$ix
  sProbes <- probes[sInd]
  sMean <- meanData[sInd]
  unik <- !duplicated(sGenes)
  uniInd <- seq_along(sGenes)[unik]  ## indices
  uniGene <- sGenes[unik] ## the values
  uniInd <- append(uniInd, 0, after = 0)
  tmpInd <- vector(mode="numeric", length=length(uniGene))
  for (i in 1:(length(uniInd)-1)){
    a <- uniInd[i]+1
    b <- uniInd[i+1]
    # maxV <- max(sMean[a:b])
    tmpInd[i] <- uniInd[i] + which.is.max(sMean[a:b])
  }
  ind <- sInd[tmpInd]
  return(list(first=ind, second=uniGene))
}


massivePCC_withoutNaN <- function (eset){
  nRow <- nrow(eset)
  nCol <- ncol(eset)
  eset[is.na(eset)] <- 0
  if (nRow > 1) {
    sumX <- rowSums(eset)
    # meanX <- colMeans(eset)
    # A * B	Element-wise multiplication
    sumV <- sqrt( rowSums(eset * eset) - (sumX * sumX)/nCol )
    PCC_mat <- matrix(0, nrow = nRow, ncol = nRow)
    
    #starting parallel
    # no_cores <- detectCores() - 1
    # cl <- makeCluster(no_cores)
    
    for (i in 1:(nRow-1)){
      if (i %% round(nRow/100) == 0) {
        print(sprintf("Processing Massive PCC ... %.2f%%",(i/nRow*100)))
      }
      # %*% Matrix multiplication
      PCC_mat[i,(i+1):nRow] <- ( eset[(i+1):nRow,] %*% eset[i,] - sumX[i] * sumX[(i+1):nRow]/nCol ) / (sumV[i] * sumV[(i+1):nRow])
      # PCC_mat[i,(i+1):nRow] <- ( tcrossprod(eset[(i+1):nRow,],t(eset[i,])) - sumX[i] * sumX[(i+1):nRow]/nCol ) / (sumV[i] * sumV[(i+1):nRow])
      # array <- tcrossprod(eset[(i+1):nRow,],t(eset[i,]))
      # array <- array - sumX[i] * sumX[(i+1):nRow]/nCol
      # denom <- (sumV[i] * sumV[(i+1):nRow])
      # array <- array/denom
      # PCC_mat[i,(i+1):nRow] <- array
    }
    PCC_mat2 <- PCC_mat + t(PCC_mat)
    
  }
  return(PCC_mat)
}