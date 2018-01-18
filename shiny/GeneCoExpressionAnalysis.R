# 01/17/2018 Zhi Huang
lmQCM_preprocess <- function(data) {
  # Step 1
  RNA <- data.frame(data[1:dim(data)[1], 2:dim(data)[2]])
  geneID <- data.frame(data[1:dim(data)[1], 1])
  
  
  return(666)
}



lmQCM <- function(cMatrix, gamma, lambda, t, beta, minClusterSize) {
  print(sprintf("gamma: %.2f",gamma))
  print(sprintf("lambda: %.2f",lambda))
  print(sprintf("t: %.2f",t))
  print(sprintf("beta: %.2f",beta))
  print(sprintf("minClusterSize: %d",minClusterSize))
  
  
  
  
  
  
  
  
  return(666)
}