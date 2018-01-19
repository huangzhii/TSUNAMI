## Step 2
load(file = "../../cMatrix.RData")
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