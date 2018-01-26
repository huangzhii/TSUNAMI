# Credit to wnchang@iu.edu 2017-10-11
# Edited by Zhi Huang 01/26/2018


#==========================first, log data===========================================
setwd("/Users/zhi/Desktop/GeneCoexpression/WGCNA"); #mac
data<-read.csv("../matlab_old/RNAdata.csv", header=T, stringsAsFactors=F)
RNA <- as.matrix(data[1:dim(data)[1], 2:dim(data)[2]])
geneID <- data[,1]
row.names(RNA) <- genename
datExpr <- t(RNA) # gene should be colnames, sample should be rownames
# datExpr <- log(datExpr + 1) # uncomment if don't need logarithm
dim(datExpr)


#=====================================================================================
#
#  Code chunk 1 : basic setting of WGCNA
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.1
workingDir = "/Users/zhi/Desktop/GeneCoexpression/WGCNA";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
allowWGCNAThreads() #enableWGCNAThreads()


#=====================================================================================
#
#  Code chunk 2 : choose the power
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(12, 9)
pdf(file = "COAD.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3 : cal net
#                 in this part, we need to pay attention to parameters which are
#                 'power', 'minModuliSize', and 'mergeCutHeight'
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,     #30,
                       reassignThreshold = 0, mergeCutHeight =  0.25,    # 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)

netcolors = net$colors
matrix<- data.frame(cbind(geneID, netcolors))

geneCharVector <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)

text <- ""
for (i in 1: length(unique(netcolors))-1){
  geneChar <- matrix[which(matrix$netcolors == i), 1]
  geneCharVector[i] <- list(geneChar)
  text <- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
}
#=====================================================================================
#
#  Code chunk 4 : plot clustring
#
#=====================================================================================


# open a graphics window
#sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main ="COAD, 
                    power= 6, minModuleSize= 30,mergeCutHeight= 0.25")

#=====================================================================================
#
#  Code chunk 5 : output and check
#
#=====================================================================================


table(net$colors)
table(mergedColors)

length(table(net$colors))
length(table(mergedColors))

colorNumber <- length(table(mergedColors))
par(mfrow = c(2,2))
#plot 
for (i in 0 : (colorNumber-1))
{
  #print(i)
  assign(paste("col_",i,sep=""), labels2colors(i))
  hist(cor(datExpr[,which(net$colors == i)]) ,xlim = c(-1,1),main = paste("Hist :color is",get(paste("col_",i,sep="")) ),col = "lightblue")
}
#save  
gene_name_list <- list()
for (i in 0: (colorNumber-1))
{
  ttt <- datExpr[,which(net$colors == i)]
  gene_name_list[[length(gene_name_list)+1]] <- colnames(ttt)	#store all gene name in order, 
  #larger group stored firstly, then small group
}										#length of this list is the number of group
#assign(paste(cancerStr,"_list",sep = ""),gene_name_list)
save( gene_name_list, file = "COAD_list.RData")

dev.off()


#big0c <- labels2colors(0)   #maybe noise ,like grey
#big1c <- labels2colors(1)
#big2c <- labels2colors(2)
#big3c <- labels2colors(3)

#par(mfrow = c(2,2))
#hist( cor(datExpr[,which(net$colors == 0)]) ,xlim = c(-1,1),main = paste("Hist :color is", big0c),col = "lightblue")
#hist( cor(datExpr[,which(net$colors == 1)]) ,xlim = c(-1,1),main = paste("Hist :color is", big1c),col = "lightblue")
#hist( cor(datExpr[,which(net$colors == 2)]) ,xlim = c(-1,1),main = paste("Hist :color is", big2c),col = "lightblue")
#hist( cor(datExpr[,which(net$colors == 3)]) ,xlim = c(-1,1),main = paste("Hist :color is", big3c),col = "lightblue")

###############################################################
#2017-11-06
#extracting top 5 clusters, compose a new matrix for down-stream analysis,
#e.g take LRR analysis for this new matrix, verify whether can get 5 low rank 

###############################################################


#'1' means turquoise, biggest cluster	
#'0' more like noise, we ignore it. The rest of group(different color) is 			
data_class1_trans <- datExpr[,which(net$colors == 1)]
data_class1 <- t(data_class1_trans)

data_class2_trans <- datExpr[,which(net$colors == 2)]
data_class2 <- t(data_class2_trans)

data_class3_trans <- datExpr[,which(net$colors == 3)]
data_class3 <- t(data_class3_trans)

data_class4_trans <- datExpr[,which(net$colors == 4)]
data_class4 <- t(data_class4_trans)

data_class5_trans <- datExpr[,which(net$colors == 5)]
data_class5 <- t(data_class5_trans)

##using rbind
data_5class <- rbind(data_class1, data_class2, data_class3, data_class4, data_class5)

#len1 <- nrow(data_class1)
#len2 <- nrow(data_class2)
#len3 <- nrow(data_class3)
#len4 <- nrow(data_class4)
#len5 <- nrow(data_class5)

#sumRow <- len1 + len2  + len3 + len4 + len5 
#sumCol <- ncol(data_class1)

#data_5class <- matrix(NA, nrow = sumRow, ncol = sumCol)

#data_5class[1:len1,] <- data_class1
#data_5class[(len1+1):(len1+len2),] <- data_class2
#data_5class[(len1+len2+1):(len1+len2+len3),] <- data_class3
#data_5class[(len1+len2+len3+1):(len1+len2+len3+len4),] <- data_class4
#data_5class[(len1+len2+len3+len4+1):(len1+len2+len3+len4+len5),] <- data_class5

#rownames(data_5class[1:nrow(data_class1),]) <- rownames(data_class1)

#save as RData
#save(data_5class, data_class1,data_class2,data_class3,data_class4,data_class5,file ="./COAD_5class.RData")
str_file <- 
  save(data_5class, file = "./COAD_5class.RData")

##OR save as csv
#write.csv(data_5class, file="./COAD_coexp_5class.csv")

#write.csv(data_class1, file="./COAD_coexp_class1.csv")
#write.csv(data_class2, file="./COAD_coexp_class2.csv")
#write.csv(data_class3, file="./COAD_coexp_class3.csv")
#write.csv(data_class4, file="./COAD_coexp_class4.csv")
#write.csv(data_class5, file="./COAD_coexp_class5.csv")




