library(shiny)
library(rsconnect)
library(plyr)
library(data.table)
library(genefilter)
library(Biobase)
library(rPython)
library(WGCNA)
library(GEOquery)

options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB
options(shiny.sanitize.errors = FALSE)
# setwd("/Users/zhi/Desktop/GeneCoexpression/RGUI"); #mac #remove when deploy to shinyapps.io
source("./lmQCM/GeneCoExpressionAnalysis.R")
# ----------------------------------------------------
data <- NULL
GEO <- NULL
finalExp <- NULL
finalSym <- NULL
finalSymChar <- NULL
text <- NULL
eigengene_matrix <- NULL
fname <- NULL # e.g. 12345_at

function(input, output, session) {
  
  observeEvent(input$action1,{
      print('tab1')
      session$sendCustomMessage("myCallbackHandler", "tab1")
  })
  output$mytable1 <- DT::renderDataTable({
    showModal(modalDialog(
      title = "Loading NCBI GEO database", footer = modalButton("OK"), easyClose = TRUE,
      div(class = "busy",
          p("Loading NCBI GEO database. Last fetched version: 01/31/2018"),
          img(src="images/loading.gif"),
          style = "margin: auto"
      )
    ))
    # NCBI GEO
    
    load_data <- function(path) { 
      files <- dir(path, pattern = '\\.csv', full.names = TRUE)
      tables <- lapply(files, read.csv)
      do.call(rbind, tables)
    }
    
    GEO_DB_csvpath = "./NCBI_GEO_SERIES_20180131_DOWNLOADED"
    GEO <<- load_data(GEO_DB_csvpath)
    GEO[["Actions"]] <- paste0('
    <div>
    <button type="button" class="btn-analysis" id=analysis_',1:nrow(GEO),'>Analyze</button></div>
    ')
    print("GEO data loaded from ./NCBI_GEO_SERIES_20180131_DOWNLOADED")
    removeModal()
    DT::datatable(GEO, extensions = 'Responsive', escape=F, selection = 'none')
  })
  
observeEvent(input$lastClickId,{
  if (input$lastClickId%like%"analysis"){
  row_of_GEO <- as.numeric(gsub("analysis_","",input$lastClickId))
  myGSE <- as.character(GEO[row_of_GEO,1])
  print(myGSE)
  
  showModal(modalDialog(
    title = sprintf("Loading %s file from NCBI GEO Database", myGSE), footer = NULL,
    div(class = "busy",
        p("We are currently loading your selected file from NCBI GEO Database ..."),
        img(src="images/loading.gif"),
        style = "margin: auto"
    )
  ))
  t <- try(gset <- getGEO(myGSE, GSEMatrix=TRUE, AnnotGPL=FALSE)) #AnnotGPL default is FALSE
  if("try-error" %in% class(t)) {
    print("HTTP error 404")
    showModal(modalDialog(
      title = "Important message", footer = modalButton("OK"),
      sprintf("%s is not available in NCBI GEO Database, please try another!",myGSE), easyClose = TRUE
    ))
    return()
  }
  
  if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  edata <- exprs(gset) #This is the expression matrix
  # pdata <- pData(gset) # data.frame of phenotypic information.
  fname <<- featureNames(gset) # e.g. 12345_at
  data <<- cbind(fname, edata)
  row.names(data) <<- seq(1, length(fname))
  
  
  updateTextInput(session, "platform_text", value = paste(annotation(gset), input$controller))
  output$summary <- renderPrint({
    print(sprintf("Number of Genes: %d",dim(edata)[1]))
    print(sprintf("Number of Samples: %d",dim(edata)[2]))
    print(sprintf("Annotation Platform: %s",annotation(gset)))
  })
  
  print("GEO file is downloaded to server and processed.")
  removeModal()
  output$mytable4 <- DT::renderDataTable({
    # Expression Value
    DT::datatable(data[1:input$quicklook_row, 1:input$quicklook_col],
                  extensions = 'Responsive', escape=F, selection = 'none')
  })
  output$mytable5 <- DT::renderDataTable({
    # Expression Value
    DT::datatable(data[input$starting_row:input$quicklook_row, input$starting_col:input$quicklook_col],
                  , extensions = 'Responsive', escape=F, selection = 'none')
  })
  output$mytable6 <- renderTable({
    # Gene ID
    data[input$starting_gene_row:input$quicklook_row, 1]
  })
  print('tab2')
  session$sendCustomMessage("myCallbackHandler", "tab2")
  }
}
)
  
  output$readingcsv <- reactive({
    print("readingcsv = 1")
    return(is.null(input$csvfile))
  })
  
  
  observeEvent(input$action2,{
      if(is.null(input$csvfile)){
        print("no files!")
        showModal(modalDialog(
          title = "Important message", footer = modalButton("OK"),
          "No file uploaded! Please retry!", easyClose = TRUE
        ))
        return(NULL)
      }
      else {
        print("Reading file.")
        showModal(modalDialog(
          title = "Intepreting uploaded file in progress", footer = NULL,
          div(class = "busy",
              p("Intepreting ..."),
              img(src="images/loading.gif"),
              style = "margin: auto"
          )
        ))
        data <<- read.csv(input$csvfile$datapath,
                          header = input$header,
                          sep = input$sep,
                          quote = input$quote)
        print("CSV file Processed.")
        removeModal()
        
        output$summary <- renderPrint({
          print(sprintf("Number of Genes: %d",dim(data)[1]))
          print(sprintf("Number of Samples: %d",(dim(data)[2]-1)))
          print("Annotation Platform: Unknown")
        })
        
        output$mytable4 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[1:input$quicklook_row, 1:input$quicklook_col],
                        , extensions = 'Responsive', escape=F, selection = 'none')
        })
        output$mytable5 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[input$starting_row:input$quicklook_row, input$starting_col:input$quicklook_col],
                        , extensions = 'Responsive', escape=F, selection = 'none')
        })
        output$mytable6 <- renderTable({
          # Gene ID
          data[input$starting_gene_row:input$quicklook_row, 1]
        })
        print('tab2')
        session$sendCustomMessage("myCallbackHandler", "tab2")
      }
  })
  
  observeEvent(input$action_platform,{
    
    showModal(modalDialog(
      title = "Converting...", footer = NULL,
      div(class = "busy",
          p("We are currently converting probe ID to Gene ID..."),
          img(src="images/loading.gif"),
          style = "margin: auto"
      )
    ))
    platform_name <- gsub(" ", "", input$platform_text, fixed = TRUE)
    
    print(sprintf("Platform: %s",platform_name))
    try(gpl <- getGEO(platform_name))
    if("try-error" %in% class(t)) {
      removeModal()
      print("HTTP error 404")
      showModal(modalDialog(
        title = "Important message", footer = modalButton("OK"),
        sprintf("%s is not available in NCBI GEO Database, please try another!", platform_name), easyClose = TRUE
      ))
      return()
    }
    print("Platform Loaded.")
    
    #https://www.rdocumentation.org/packages/GEOquery/versions/2.38.4/topics/GEOData-class
    gpltable <- Table(gpl)
    fname2 <- fname
    for (i in 1:length(fname)){
      fname2[i] <- gpltable$`Gene Symbol`[which(gpltable$ID == fname[i])]
    }
    print(dim(data))
    print(length(fname2))
    data[,1] <<- fname2
    # row.names(data) <- seq(1, length(fname2))
    output$mytable4 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data[1:input$quicklook_row, 1:input$quicklook_col],
                    extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable5 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data[input$starting_row:input$quicklook_row, input$starting_col:input$quicklook_col],
                    extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable6 <- renderTable({
      # Gene ID
      data[input$starting_gene_row:input$quicklook_row, 1]
    })
    print("Platform Conversion finished")
    
    removeModal()
  })
  
  
  observeEvent(input$action3,{
      # source("/Users/zhi/Desktop/GeneCoexpression/RGUI/utils.R")
      source("utils.R")
      
      showModal(modalDialog(
        title = "Cleaning input data", footer = NULL,
        div(class = "busy",
            p("Cleaning ..."),
            img(src="images/loading.gif"),
            style = "margin: auto"
        )
      ))
      print(dim(data))
      # Step 0
      RNA <- as.matrix(data[input$starting_row:dim(data)[1], input$starting_col:dim(data)[2]])
      class(RNA) <- "numeric"
      geneID <- data.frame(data[input$starting_gene_row:dim(data)[1], 1])
      print(dim(RNA))
      print(dim(geneID))
      # Remove data with lowest 20% absolute exp value shared by all samples
      percentile <- input$absolute_expval/100.
      
      if (percentile > 0){
        RNA_filtered1 = RNA[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
        geneID_filtered1 = geneID[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
      }
      else {
        RNA_filtered1 = RNA
        geneID_filtered1 = geneID
      }
      # Remove data with lowest 10% variance across samples
      percentile <- input$variance_expval/100.
      if (percentile > 0){
        index <- varFilter2(eset = RNA_filtered1, var.cutoff = percentile)
        RNA_filtered2 = RNA_filtered1[index, ]
        geneID_filtered2 = geneID_filtered1[index]
      }
      else {
        RNA_filtered2 = RNA_filtered1
        geneID_filtered2 = geneID_filtered1
      }
      
      # expData <- RNA_filtered2
      # res <- highExpressionProbes(geneID_filtered2, geneID_filtered2, expData)
      # ind1 <- res$first
      # uniGene <- as.character(res$second)
      # tmpExp <- expData[ind1,]
      uniGene <- geneID_filtered2
      tmpExp <- RNA_filtered2
      
      if (input$checkbox_NA){
        tmpExp[is.na(tmpExp)] <- 0
      }
      if (input$checkbox_empty){
        print(sprintf("data dimension before remove gene with empty symbol: %d x %d",dim(tmpExp)[1],dim(tmpExp)[2]))
        uniGene_temp <- subset(uniGene, nchar(as.character(uniGene)) > 0)
        tmpExp <- subset(tmpExp, nchar(as.character(uniGene)) > 0)
        uniGene <- uniGene_temp
        print(sprintf("data dimension after remove gene with empty symbol: %d x %d",dim(tmpExp)[1],dim(tmpExp)[2]))
      }
      if (input$checkbox_duplicated){
        row2remove <- numeric()
        finalSymCharTable <- table(uniGene)
        for (i in 1:length(finalSymCharTable)){
          if (as.numeric(finalSymCharTable[i]) > 1){ # if exist duplicated Gene
            genename <- names(finalSymCharTable[i])
            idx_with_max_mean <- which.max(rowMeans(tmpExp[which(uniGene == genename),]))
            # print(idx_with_max_mean)
            row2remove <- c( row2remove, (which(uniGene == genename)[-idx_with_max_mean]) )
          }
        }
        if (length(row2remove) > 0){ # Otherwise numerical(0) will remove all data in tmpExp
          tmpExp <- tmpExp[-row2remove,]
          uniGene <- uniGene[-row2remove]
        }
        print(sprintf("data dimension after remove duplicated gene symbol: %d x %d",dim(tmpExp)[1],dim(tmpExp)[2]))
        
      }
      nSample <- ncol(tmpExp)
      res <- sort.int(rowMeans(tmpExp), decreasing = TRUE, index.return=TRUE)
      sortMean <- res$x
      sortInd <- res$ix
      topN <- min(input$max_gene_retain, nrow(tmpExp))
      finalExp <<- tmpExp[sortInd[1:topN], ]
      finalSym <<- uniGene[sortInd[1:topN]]
      finalSymChar <<- as.character(finalSym)
      removeModal()
      print('tab3')
      session$sendCustomMessage("myCallbackHandler", "tab3")
  })
  
  
  observeEvent(input$action4_lmQCM,{
      #lmQCM
      showModal(modalDialog(
        title = "Using lmQCM to calculate merged clusters", footer = NULL,
        div(class = "busy",
            p("Calculating ..."),
            img(src="images/loading.gif"),
            style = "margin: auto"
        )
      ))
      step1 = 1
      gamma = input$gamma
      t = input$t
      lambda = input$lambda
      beta = input$beta
      minClusterSize = input$minClusterSize
      massiveCC = input$massiveCC
      system("which python") # python version should be from Anaconda2
      # system("/Users/zhi/anaconda2/bin/python main.py")
      python.load("main_old.py")
      mergedCluster <- python.call("mainroutine", step1, as.vector(finalExp), nrow(finalExp), ncol(finalExp), gamma, t, lambda, beta, minClusterSize, input$massiveCC)
      geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
      temp_eigengene <- matrix(0, nrow = length(mergedCluster), ncol = dim(finalExp)[2]) # Clusters * Samples
      
      temptext <- ""
      for (i in 1:(length(mergedCluster))) {
        vector <- as.matrix(mergedCluster[[i]])
        vector <- vector + 1 # covert python indexing to R indexing
        geneID <- vector
        print(i)
        print(vector)
        # ===== Calculate Eigengene Start
        X <- finalExp[geneID,]
        mu <- rowMeans(X)
        stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
        #normalize X:
        XNorm <- sweep(X,1,mu)
        XNorm <- apply(XNorm, 2, function(x) x/stddev)
        SVD <- svd(XNorm, LINPACK = FALSE)
        temp_eigengene[i,] <- t(SVD$v[,1])
        # ===== Calculate Eigengene Finished
        
        geneChar <- c(toString(i), finalSymChar[vector])
        geneCharVector[i] <- list(geneChar)
        temptext <- paste(temptext, capture.output(cat(geneChar, sep=',')), sep="\n")
      }
      temptext <- substring(temptext, 2) # remove first \n separater
      text <<- temptext
      eigengene_matrix <<- temp_eigengene
      
      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))
      
      output$clusterResult <- renderTable({
        return(geneCharVector2)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)
      
      output$mytable7 <- renderTable({
        return(eigengene_matrix)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)
      
      removeModal()
      print('tab4')
      session$sendCustomMessage("myCallbackHandler", "tab4")
  })
  observeEvent(input$checkPower, {
    if (length(finalSym) > 0){
      #WGCNA
      showModal(modalDialog(
        title = "WGCNA step 1: Preview the power", footer = NULL,
        div(class = "busy",
            p("Calculating ..."),
            img(src="images/loading.gif"),
            style = "margin: auto"
        )
      ))
      
      row.names(finalExp) <- finalSym
      datExpr <- t(finalExp) # gene should be colnames, sample should be rownames
      # datExpr <- log(datExpr + 1) # uncomment if don't need logarithm
      print(dim(datExpr))
      # The following setting is important, do not omit.
      options(stringsAsFactors = FALSE);
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
      cex1 = 0.9;
      # Scale-free topology fit index as a function of the soft-thresholding power
      
      output$WGCNAPowerPlot1and2 <- renderUI({
        fluidRow(
          column(6, plotOutput('WGCNAPowerPlot1', height = "300px")),
          column(6, plotOutput('WGCNAPowerPlot2', height = "300px"))
        )
      })
      output$WGCNAPowerPlot1 <- renderPlot({
        plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence"))
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red")
        # this line corresponds to using an R^2 cut-off of h
        abline(h=0.9,col="blue")
      })
      output$WGCNAPowerPlot2 <- renderPlot({
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      })
      
      
      removeModal()
    }
  })
  
  observeEvent(input$action4_WGCNA,{
      #WGCNA
      showModal(modalDialog(
        title = "Using WGCNA to calculate merged clusters", footer = NULL,
        div(class = "busy",
            p("Calculating ..."),
            img(src="images/loading.gif"),
            style = "margin: auto"
        )
      ))
      
      row.names(finalExp) <- finalSym
      datExpr <- t(finalExp) # gene should be colnames, sample should be rownames
      # datExpr <- log(datExpr + 1) # uncomment if don't need logarithm
      # The following setting is important, do not omit.
      options(stringsAsFactors = FALSE);
      allowWGCNAThreads() #enableWGCNAThreads()
      
      #=====================================================================================
      #
      #  Code chunk 3 : cal net
      #                 in this part, we need to pay attention to parameters which are
      #                 'power', 'minModuliSize', and 'mergeCutHeight'
      #
      #=====================================================================================
      
      net = blockwiseModules(datExpr, power = input$power,
                             TOMType = "unsigned", minModuleSize = input$minModuleSize,     #30,
                             reassignThreshold = input$reassignThreshold, mergeCutHeight =  input$mergeCutHeight,    # 0.25,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             saveTOMs = FALSE,
                             saveTOMFileBase = "femaleMouseTOM", 
                             verbose = input$verbose)
      
      netcolors = net$colors
      matrixdata<- data.frame(cbind(finalSymChar, netcolors))
      geneCharVector <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)
      
      #=====================================================================================
      #
      #  Code chunk 4 : plot clustring
      #
      #=====================================================================================
      
      
      # open a graphics window
      #sizeGrWindow(12, 9)
      # Convert labels to colors for plotting
      mergedColors = labels2colors(net$colors)
      # print("Merged Colors:")
      # print(mergedColors)
      # print("net$blockGenes[[1]]:")
      # print(net$blockGenes[[1]])
      
      # Plot the dendrogram and the module colors underneath
      output$WGCNAresultUI <- renderUI({
        plotOutput('WGCNAresult', height = "500px")
      })
      output$WGCNAresult <- renderPlot({
          plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = sprintf("power= %d, minModuleSize= %d,mergeCutHeight= %.4f",
                                         input$power, input$minModuleSize, input$mergeCutHeight))
      })
      
      removeModal()
      geneCharVector <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)
      print("unique netcolors: ")
      print(unique(netcolors))
      # print(matrixdata)
      # print(finalSym)
      temptext <- ""
      temp_eigengene <- matrix(0, nrow = length(unique(netcolors))-1, ncol = dim(finalExp)[2]) # Clusters * Samples
      for (i in 1: (length(unique(netcolors))-1) ){
        geneID <- which(matrixdata$netcolors == i)
        # ===== Calculate Eigengene Start
        X <- finalExp[geneID,]
        mu <- rowMeans(X)
        stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
        #normalize X:
        XNorm <- sweep(X,1,mu)
        XNorm <- apply(XNorm, 2, function(x) x/stddev)
        SVD <- svd(XNorm, LINPACK = FALSE)
        temp_eigengene[i,] <- t(SVD$v[,1])
        # ===== Calculate Eigengene Finished
        geneChar <- c(toString(i), matrixdata[geneID, 1])
        geneCharVector[i] <- list(geneChar)
        temptext <- paste(temptext, capture.output(cat(geneChar, sep=',')), sep="\n")
      }
      temptext <- substring(temptext, 2) # remove first \n separater
      text <<- temptext
      eigengene_matrix <<- temp_eigengene
      
      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))
      
      output$clusterResult <- renderTable({
        return(geneCharVector2)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)
      
      output$mytable7 <- renderTable({
        return(eigengene_matrix)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)
      
      print('tab4')
      session$sendCustomMessage("myCallbackHandler", "tab4")
  })
  
  output$downloadData1 <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      name = "mergedCluster"
      # name = paste("lmQCMresult","gamma",gamma(),"lambda",lambda(),"t",t(),"beta",beta(),"minClusterSize",minClusterSize(), sep = "_", collapse = NULL)
      paste(name, input$filetype1, sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype1, "csv" = 0, "txt" = 1)
      if (sep == 0){
        text_download = gsub(",", ",", noquote(text))
        write.table(text_download, file, eol = "\r\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
      if (sep == 1){
        text_download = gsub(",", "\t", noquote(text))
        write.table(text_download, file, eol = "\r\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
    }
  )
  output$downloadData2 <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      name = "EigengeneMatrix"
      # name = paste("lmQCMresult","gamma",gamma(),"lambda",lambda(),"t",t(),"beta",beta(),"minClusterSize",minClusterSize(), sep = "_", collapse = NULL)
      paste(name, input$filetype2, sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      separator <- switch(input$filetype2, "csv" = ',', "txt" = '\t')
      write.table(eigengene_matrix, file = file, append = FALSE, quote = TRUE, sep = separator,
                  eol = "\r\n", na = "NA", dec = ".", row.names = F,
                  col.names = F, qmethod = c("escape", "double"),
                  fileEncoding = "")
    }
  )
}
