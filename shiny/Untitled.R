load("./data/GEO_20180424.Rdata")
library(stringi)
GSElist = list()
for (i in 1:dim(GEO)[1]){
  if (stri_detect_fixed(GEO[i, 2], "microarray")) {
    
    myGSE = levels(droplevels(GEO[i,1]))
    t <- try(gset <- getGEO(myGSE, GSEMatrix=TRUE, AnnotGPL=FALSE)) #AnnotGPL default is FALSE
    if("try-error" %in% class(t)) {
      next()
    }
    if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
    if (length(idx) == 0){
      next()
    }
    gset <- gset[[idx]]
    edata <- exprs(gset) #This is the expression matrix
    if (dim(edata)[1] == 0 || is.null(dim(edata))){
      next()
    }
    print(myGSE)
    GSElist = c(GSElist, myGSE)
    
  }
  
}
