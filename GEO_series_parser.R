library(GEOquery)
library(stringr)
path = "/home/zhihuan/Desktop/TSUNAMI/NCBI_GEO_SERIES_20190303_DOWNLOADED/"
series_all = NULL
for (f in dir(path)){
  series = read.csv(paste0(path,f))
  series_all = rbind(series_all, series)
}


order =  sort.int(as.numeric(str_extract(series_all$Accession, "[[:digit:]]+")), index.return = T)$ix
series_all = series_all[order,]
GEO = series_all
save(GEO, file = paste0(path,"../GEO_20190303.Rdata"))



##########################################
gset_list = NULL
i = 0
for (accession in series_all$Accession){
  message("----------------")
  message(accession)
  message("----------------")
  i = i+1
  t <- try(gset <- getGEO(accession, GSEMatrix=TRUE, AnnotGPL=FALSE)) #AnnotGPL default is FALSE
  if("try-error" %in% class(t)) next
  gset_list[[i]] = gset
  
  if (length(gset) == 0){
    print("Object gset doesn't contain data.")
    next
  }
  
  # edata <- exprs(gset) #This is the expression matrix
  # if (dim(edata)[1] == 0 || is.null(dim(edata))){
  #   print("No expression data")
  #   next
  # }
  # # pdata <- pData(gset) # data.frame of phenotypic information.
  # fname = featureNames(gset) # e.g. 12345_at
  # data = cbind(fname, edata)
  # row.names(data) <- seq(1, length(fname))
  # 
  # annotation(gset)
  # 
  # series_all[series_all$Accession == accession, ]$dim1 = dim(data)[1]
  # series_all[series_all$Accession == accession, ]$dim2 = dim(data)[2]
}
