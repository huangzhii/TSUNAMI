# Zhi Huang 04/09/2018

library(circlize)
library(openxlsx)
setwd("~/Desktop/TSUNAMI/shiny")
data <- read.xlsx("./MarkData.xlsx", sheet = 1, startRow = 1, colNames = T, rowNames = T)
# genes_str <- c("LOC102725121", "FAM138A", "RIMS2", "LINC01128", "MMP23A", "ULK4P1")
genes_str <- data[,1]




# import hg19 and hg38
load("./data/UCSC_hg19_refGene_20180330.Rdata") # varname: hg19
load("./data/UCSC_hg38_refGene_20180330.Rdata") # varname: hg38
hg19 <- data.frame(cbind(rownames(hg19), hg19, hg19[6]-hg19[5]))
hg38 <- data.frame(cbind(rownames(hg38), hg38, hg38[6]-hg38[5]))
colnames(hg38) = c("id","","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID","","","","length")
colnames(hg19) = c("id","","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","proteinID","alignID","","","","length")
hg19.ring <- hg19[!grepl("_", hg19$chrom),] # remove undefined chromosome
hg38.ring <- hg38[!grepl("_", hg38$chrom),]
hg19.ring <- hg19.ring[!grepl("chrM", hg19.ring$chrom),]
hg38.ring <- hg38.ring[!grepl("chrM", hg38.ring$chrom),]
hg19.matched <- hg19.ring[match(genes_str, hg19.ring$alignID, nomatch = 0), ]
hg38.matched <- hg38.ring[match(genes_str, hg38.ring$alignID, nomatch = 0), ]
hg19.ring.lengthsum <- aggregate(hg19.ring["length"],hg19.ring["chrom"],sum)
hg38.ring.lengthsum <- aggregate(hg38.ring["length"],hg38.ring["chrom"],sum)

factors_count = as.data.frame(hg38.ring.lengthsum)
factors = factor(factors_count[,1], levels = factors_count[,1])
xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
rownames(xlim) = factors_count[,1]
BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
BED.data$txStart <- as.numeric(sub('.*\\:', '', data$Location ))
BED.data$Type <- data$Type
BED.data$Alteration <- data$Alteration
BED.data$Drug <- data$Drug
BED.data$Location <- data$Location
BED.data$IDandDrug <- paste0(data$Gene, " (", data$Drug, ")")
for(i in 1:dim(data)[1]){
  if(is.na(data$Drug[i])){
    BED.data$IDandDrug[i] <- data$Gene[i]
  }
}






circlizeGenomics2 <- function(BED.data, factors, xlim, myTitle){
  par(mar = c(1, 1, 1, 1))
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0))
  # circos.initializeWithIdeogram(plotType = c("axis", "labels"))
  circos.initializeWithIdeogram(plotType = NULL)
  circos.genomicLabels(BED.data, labels.column = 10, side = "outside",
                       col = "black", line_col = "blue", cex = 0.5)
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
    circos.text(mean(xlim), mean(ylim), chr, cex = 0.4, col = "white", facing = "inside", niceFacing = TRUE)
  }, track.height = 0.08, bg.border = NA)

  title(myTitle)
  BED.data.link = BED.data[BED.data$Type == "FUSION",]
  for (i in 1:dim(BED.data.link)[1]){
    location = gsub(" ", "", BED.data.link$Location[i], fixed = TRUE)
    location = unlist(strsplit(location, "-"))
    chrom1 = unlist(strsplit(location[1], ":"))[1]
    pt1 = as.numeric(unlist(strsplit(location[1], ":")))[2]
    chrom2 = unlist(strsplit(location[2], ":"))[1]
    pt2 = as.numeric(unlist(strsplit(location[2], ":")))[2]
    circos.link(sector.index1=paste("chr", chrom1, sep=""), point1=pt1,
                sector.index2=paste("chr", chrom2, sep=""), point2=pt2,
                col = "pink", lwd = 5)
    # R color: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
  }
  circos.genomicTrack(BED.data, track.height = 0.05, bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
                        i = getI(...)
                        for(i in 1:dim(value)[1]){
                          current_value = value[i,]
                          current_region = region[i,]
                          if(current_value$Type == "CNA" & current_value$Alteration == "AMP"){
                            circos.genomicPoints(current_region, value = current_value, pch = 8, col="darkgreen", ...)
                          }
                          else if(current_value$Type == "CNA" & current_value$Alteration == "HOMDEL"){
                            circos.genomicPoints(current_region, value = current_value, pch = 8, col="red", ...)
                          }
                          else if(current_value$Type == "Inframe"){
                            circos.genomicPoints(current_region, value = current_value, pch = 16, col="darkgreen", ...)
                          }
                          else if(current_value$Type == "Missense"){
                            circos.genomicPoints(current_region, value = current_value, pch = 16, col="black", ...)
                          }
                          else if(current_value$Type == "Truncation"){
                            circos.genomicPoints(current_region, value = value[i,], pch = 16, col="red", ...)
                          }
                        }
                      })
  circos.genomicTrack(BED.data, track.height = 0.05, bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
                        i = getI(...)
                        for(i in 1:dim(value)[1]){
                          current_value = value[i,]
                          current_region = region[i,]
                          if(current_value$Type == "EXP mRNA" & current_value$Alteration == "UP"){
                            circos.genomicPoints(current_region, value = current_value, pch = 24, col="darkgreen", ...)
                          }
                          else if(current_value$Type == "EXP mRNA" & current_value$Alteration == "DOWN"){
                            circos.genomicPoints(current_region, value = current_value, pch = 25, col="red", ...)
                          }
                        }
                      })
  circos.genomicTrack(BED.data, track.height = 0.05, bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
                        i = getI(...)
                        for(i in 1:dim(value)[1]){
                          current_value = value[i,]
                          current_region = region[i,]
                          if(current_value$Type == "EXP Protein" & current_value$Alteration == "UP"){
                            circos.genomicPoints(current_region, value = current_value, pch = 15, col="darkgreen", ...)
                          }
                          else if(current_value$Type == "EXP Protein" & current_value$Alteration == "DOWN"){
                            circos.genomicPoints(current_region, value = current_value, pch = 15, col="red", ...)
                          }
                        }
                      })
  
  
  
  
    
  
  circos.clear()
}
circlizeGenomics2(BED.data, factors, xlim, "Cancer Genome")
