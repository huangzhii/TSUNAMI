localMaximumQCM <- function (cMatrix, gamma, t, lambda){
  C <- list()
  nRow <- nrow(cMatrix)
  maxV <- apply(cMatrix, 2, max)
  maxInd <- apply(cMatrix, 2, which.max) # several diferrences comparing with Matlab results
  # Step 1 - find the local maximal edges
  maxEdges <- matrix(0, nrow = 0, ncol = 2)
  maxW <- matrix(0, nrow = 0, ncol = 1)
  for (i in 1:nRow){
    if (maxV[i] == max(cMatrix[maxInd[i], ])) {
      maxEdges <- rbind(maxEdges, c(maxInd[i], i))
      maxW <- rbind(maxW, maxV[i])
    }
  }
  
  res <- sort.int(maxW, decreasing = TRUE, index.return=TRUE)
  sortMaxV <- res$x
  sortMaxInd <- res$ix
  sortMaxEdges <- maxEdges[sortMaxInd,]
  print(sprintf("Number of Maximum Edges: %d", length(sortMaxInd)))
  
  currentInit <- 1
  noNewInit <- 0
  
  nodesInCluster <- matrix(0, nrow = 0, ncol = 1)
  while ((currentInit <= length(sortMaxInd)) & (noNewInit == 0)) {
    if (sortMaxV[currentInit] < (gamma * sortMaxV[1]) ) {
      noNewInit <- 1
    }
    else {
      if ( (is.element(sortMaxEdges[currentInit, 1], nodesInCluster) == FALSE) & is.element(sortMaxEdges[currentInit, 2], nodesInCluster) == FALSE) {
        newCluster <- sortMaxEdges[currentInit, ]
        addingMode <- 1
        currentDensity <- sortMaxV[currentInit]
        nCp <- 2
        totalInd <- 1:nRow
        remainInd <- setdiff(totalInd, newCluster)
        # C = setdiff(A,B) for vectors A and B, returns the values in A that 
        # are not in B with no repetitions. C will be sorted.
        while (addingMode == 1) {
          neighborWeights <- colSums(cMatrix[newCluster, remainInd])
          maxNeighborWeight <- max(neighborWeights)
          maxNeighborInd <- which.is.max(neighborWeights)
          c_v = maxNeighborWeight/nCp;
          alphaN = 1 - 1/(2*lambda*(nCp+t));
          if (c_v >= alphaN * currentDensity) {
            newCluster <- c(newCluster, remainInd[maxNeighborInd])
            nCp <- nCp + 1
            currentDensity <- (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2)
            remainInd = setdiff(remainInd, remainInd[maxNeighborInd]);
          }
          else {
            addingMode <- 0
          }
        }
        nodesInCluster <- c(nodesInCluster, newCluster)
        C <- c(C, list(newCluster))
      }
    }
    currentInit <- currentInit + 1
  }
  return(C)
}