# 01/18/2018 Zhi Huang
library(genefilter)
library(nnet)
library(parallel)
library(circlize)

varFilter2 <- function (eset, var.cutoff = 0.5, filterByQuantile = TRUE) {
  # if (deparse(substitute(var.func)) == "IQR") {
  #   vars <- genefilter:::rowIQRs(eset)
  # }
  vars <- apply(eset, 1, var)
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

# smartModal is for showing processing windows. To close a smartModal, a "removeModal()" should followed.
smartModal <- function(error=c(T,F), title = "Title", content = "Content"){
  if(error){
    showModal(modalDialog(
      title = title, footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p(content),
          style = "margin: auto; text-align: center"
      )
    ))
  }
  else{
    showModal(modalDialog(
      title = title, footer = NULL,
      div(class = "busy",
          p(content),
          img(src="https://cdn.dribbble.com/users/503653/screenshots/3143656/fluid-loader.gif"),
          style = "margin: auto; text-align: center"
      )
    ))
  }
}

circlizeGenomics <- function(BED.data, factors, xlim, mySpecies, myTitle, circos_param_size, circos_param_genelink, circos_param_genesymbol){
  # save(BED.data, file ="~/Desktop/BEDdata.Rdata")
  par(mar = c(1, 1, 1, 1))
  # reference: http://zuguang.de/circlize_book/book/initialize-genomic-plot.html#initialize-cytoband
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0))
  # circos.initializeWithIdeogram(plotType = NULL)
  if (mySpecies=="hg19") {
    circos.initializeWithIdeogram(species = mySpecies)
  }
  else if (mySpecies=="hg38"){
    circos.initializeWithIdeogram(species = mySpecies, chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  }
  # text(0, 0, "Human Chromosomes", cex = 1)
  # Abbreviations of species. e.g. hg19 for human, mm10 for mouse.
  # If this value is specified, the function will download cytoBand.txt.gz
  # from UCSC website automatically. If there is no cytoband for user's species,
  # it will keep on trying to download chromInfo file.
  
  # we assume data is simply a data frame in BED format
  # (where the first column is the chromosome name, the
  # second and third column are start and end positions,
  # and the following columns are associated values)
  
  # circos.genomicRainfall(data.frame(hg38.ring[,c(4,6:7)]), pch = 16, cex = 0.4, col = "#FF000080")
  title(myTitle)
  circos.trackPlotRegion(factors = factors, ylim = c(0, 1), bg.border = NA,
                         bg.col = rep("grey", 24), track.height = 0.05,
                         panel.fun = function(x, y) {
                           sector.name = get.cell.meta.data("sector.index")
                           xlim = get.cell.meta.data("xlim")
                           # circos.text(mean(xlim), 2.5,
                           #             facing = "outside", niceFacing = T,
                           #             sector.name, cex = 0.8, adj = c(0.5, 0))
                         })
  circos.genomicTrack(BED.data, track.height = 0.01, bg.border = NA,
                      panel.fun = function(region, value, ...) {
                        cex = (value[[1]] - min(value[[1]]))/(max(value[[1]]) - min(value[[1]]))
                        i = getI(...)
                        # numeric.column is automatically passed to `circos.genomicPoints()`
                        circos.genomicPoints(region, value = 1, ...)
                      })
  if(circos_param_genelink){
    # circos.initializeWithIdeogram(plotType = NULL)
    for (i in 1:length(BED.data$chrom)){
      for (j in 1:length(BED.data$chrom)){
        if(j<i){
          circos.link(sector.index1=BED.data$chrom[i], point1=c(BED.data$txStart[i],BED.data$txEnd[i]),
                      sector.index2=BED.data$chrom[j], point2=c(BED.data$txStart[j],BED.data$txEnd[j]),
                      col = "cadetblue1")
          # R color: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
        }
      }
    }
  }
  if(circos_param_genesymbol){
    circos.genomicLabels(BED.data, labels.column = 5, side = "inside",
                       col = as.numeric(factor(BED.data[[1]])), line_col = as.numeric(factor(BED.data[[1]])))
  }
  
  
  circos.clear()
}