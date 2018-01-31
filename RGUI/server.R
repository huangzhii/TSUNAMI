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
setwd("/Users/zhi/Desktop/GeneCoexpression/RGUI"); #mac
source("../lmQCM/GeneCoExpressionAnalysis.R")
# ----------------------------------------------------
data <- NULL
finalExp <- NULL
finalSym <- NULL
finalSymChar <- NULL
text <- ""

function(input, output, session) {
  observe({
    if(input$action1 > 0){
      
      print('tab1')
      session$sendCustomMessage("myCallbackHandler", "tab1")
    }
  })
  output$readingcsv <- reactive({
    print("readingcsv = 1")
    return(is.null(input$csvfile))
  })
  observe({
    if(input$action2 > 0) {
      
      if(is.null(input$csvfile)){
        print("no files!")
        showModal(modalDialog(
          title = "Important message", footer = modalButton("OK"),
          "No csv file uploaded! Please retry!", easyClose = TRUE
        ))
        return(NULL)
      }
      else {
        print("Reading CSV file.")
        showModal(modalDialog(
          title = "Intepreting CSV in progress", footer = NULL,
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
        output$mytable4 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[1:input$quicklook_row, 1:input$quicklook_col])
        })
        output$mytable5 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[input$starting_row:input$quicklook_row, input$starting_col:input$quicklook_col])
        })
        output$mytable6 <- renderTable({
          # Gene ID
          data[input$starting_gene_row:input$quicklook_row, 1]
        })
        print('tab2')
        session$sendCustomMessage("myCallbackHandler", "tab2")
      }
    }
  })
  observe({
    if(input$action3 > 0){
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
      geneID <- data.frame(data[input$starting_gene_row:dim(data)[1], 1])
      print(dim(RNA))
      print(dim(geneID))
      # Remove data with lowest 20% absolute exp value shared by all samples
      percentile <- input$absolute_expval/100.
      RNA_filtered1 = RNA[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
      geneID_filtered1 = geneID[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
      
      # Remove data with lowest 10% variance across samples
      percentile <- input$variance_expval/100.
      index <- varFilter2(eset = RNA_filtered1, var.cutoff = percentile)
      RNA_filtered2 = RNA_filtered1[index, ]
      geneID_filtered2 = geneID_filtered1[index]
      
      expData <- RNA_filtered2
      res <- highExpressionProbes(geneID_filtered2, geneID_filtered2, expData)
      ind1 <- res$first
      uniGene <- res$second
      tmpExp <-expData[ind1,]
      nSample <- ncol(tmpExp)
      res <- sort.int(rowMeans(tmpExp), decreasing = TRUE, index.return=TRUE)
      sortMean <- res$x
      sortInd <- res$ix
      topN <- min(input$max_gene_retain, nrow(tmpExp))
      finalExp <<- tmpExp[sortInd[1:topN], ]
      
      print(nrow(tmpExp))
      finalSym <<- uniGene[sortInd[1:topN]]
      finalSymChar <<- as.character(finalSym)
      if (input$NAconverter){
        finalExp[is.na(finalExp)] <- 0
      }
      removeModal()
      print('tab3')
      session$sendCustomMessage("myCallbackHandler", "tab3")
    }
  })
  
  observe({
    if(input$action4_lmQCM > 0){
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
      python.load("main.py")
      mergedCluster <- python.call("mainroutine", step1, as.vector(finalExp), nrow(finalExp), ncol(finalExp), gamma, t, lambda, beta, minClusterSize)
      geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
      
      for (i in 1:(length(mergedCluster))) {
        vector <- as.matrix(mergedCluster[[i]])
        vector <- vector + 1 # covert python indexing to R indexing
        geneChar <- finalSymChar[vector]
        geneCharVector[i] <- list(geneChar)
        text <<- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
      }
      geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
      
      text <- ""
      for (i in 1:(length(mergedCluster))) {
        vector <- as.matrix(mergedCluster[[i]])
        vector <- vector + 1 # covert python indexing to R indexing
        geneChar <- finalSymChar[vector]
        geneCharVector[i] <- list(geneChar)
        text <- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
      }
      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))
      
      output$clusterResult <- renderTable({
        return(geneCharVector2)
      },rownames = TRUE, colnames = FALSE, na = "", bordered = TRUE)
      
      removeModal()
      print('tab4')
      session$sendCustomMessage("myCallbackHandler", "tab4")
    }
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
      # Plot the results:
      sizeGrWindow(12, 9)
      pdf(file = "COAD.pdf")
      par(mfrow = c(1,2));
      cex1 = 0.9;
      # Scale-free topology fit index as a function of the soft-thresholding power
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
  observe({
    if(input$action4_WGCNA > 0){
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
      matrix<- data.frame(cbind(finalSym, netcolors))
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
      # Plot the dendrogram and the module colors underneath
      
      output$WGCNAresult <- renderPlot({
          plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05,
                          main = sprintf("COAD, 
                          power= %d, minModuleSize= 30,mergeCutHeight= 0.25", input$power))
      })
      
      removeModal()
      geneCharVector <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)
      
      text <- ""
      for (i in 1: length(unique(netcolors))-1){
        geneChar <- matrix[which(matrix$netcolors == i), 1]
        geneCharVector[i] <- list(geneChar)
        text <- paste(text, capture.output(cat(geneChar, sep=' ')), sep="\n")
      }
      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))
      
      output$clusterResult <- renderTable({
        return(geneCharVector2)
      },rownames = TRUE, colnames = FALSE, na = "", bordered = TRUE)
      
      print('tab4')
      session$sendCustomMessage("myCallbackHandler", "tab4")
    }
  })
  
}
