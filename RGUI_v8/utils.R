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

getFileNameExtension <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split='/')[[1]]   
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
  # the extention must be the suffix of a non-empty name    
  ext
}