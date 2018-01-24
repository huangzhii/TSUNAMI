## Step 2 (around 10 min in Y710 windows)
# setwd("E:/GeneCoexpression/shiny"); #win
setwd("/Users/zhi/Desktop/GeneCoexpression/shiny"); #mac
load(file = "../../cMatrix.RData")
source("utils.R")
source("localMaximumQCM.R")


gamma <- 0.55
lambda <- 1
t <- 1
beta <- 0.4
minClusterSize <- 10

print(sprintf("gamma: %.2f",gamma))
print(sprintf("lambda: %.2f",lambda))
print(sprintf("t: %.2f",t))
print(sprintf("beta: %.2f",beta))
print(sprintf("minClusterSize: %d",minClusterSize))


C <- localMaximumQCM(abs(cMatrix), gamma, t, lambda)
C_save <- C
C <- C_save

# Step 3 - Merge the overlapped networks

sizeC <- matrix(0, nrow = 0, ncol = length(C))
for (i in 1:length(C)){
  sizeC[i] <- length(C[[i]])
}
res <- sort.int(sizeC, decreasing = TRUE, index.return=TRUE)
sortC <- res$x
sortInd <- res$ix

C <- C[sortInd] # Still C, but sorted based on number of elements in each cell

ind <- which(sortC >= minClusterSize)

mergedCluster <- C[ind]
mergeOccur <- 1
currentInd <- 0

print("start merge")
while (mergeOccur == 1) {
  mergeOccur <- 0
  while (currentInd < length(mergedCluster)){
    currentInd <- currentInd + 1
    if (currentInd < length(mergedCluster)){
      keepInd <- 1:currentInd
      for (j in (currentInd+1) : length(mergedCluster)) {
        interCluster <- intersect(mergedCluster[[currentInd]], mergedCluster[[j]]);
        if (length(interCluster) >= beta*min(length(mergedCluster[[j]]), length(mergedCluster[[currentInd]]))) {
          mergedCluster[currentInd] <- list(union(mergedCluster[[currentInd]], mergedCluster[[j]]))
          mergeOccur <- 1
        }
        else {
          keepInd <- c(keepInd, j)
        }
      }
      mergedCluster <- mergedCluster[keepInd]
      print(sprintf("The length of merged Cluster: %d", length(mergedCluster)))
    }
  }
  sizeMergedCluster <- matrix(0, nrow = 0, ncol = length(mergedCluster))
  for (i in 1 : length(mergedCluster)) {
    sizeMergedCluster[i] <- length(mergedCluster[[i]])
  }
  res <- sort.int(sizeMergedCluster, decreasing = TRUE, index.return=TRUE)
  sortSize <- res$x
  sortMergedInd <- res$ix
  mergedCluster <- mergedCluster[sortMergedInd]
  currentInd <- 0
}
