library(GEOquery)
setwd("~/Desktop/TSUNAMI/shiny/")
load("./data/GEO_20190303.Rdata")


annotation = NULL
bad.gse = NULL
i = 0
for (gse in GEO$Accession){
  i = i+1
  if (i > 10000){
    return()
  }
  print(gse)
  t <- try(gset <- getGEO(gse, GSEMatrix=TRUE, AnnotGPL=FALSE))
  if("try-error" %in% class(t)) {
    print("error occured")
    bad.gse = c(bad.gse, gse)
  } else {
    if (length(gset) == 0){
      print("This GSE accession doesn't contain any series matrix.")
      bad.gse = c(bad.gse, gse)
    }
    gset = gset[[1]] # only pick the first gset for testing
    annotation = c(annotation, gset@annotation)
  }
  rm(gset)
}

