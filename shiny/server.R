library(shiny)
library(rsconnect)
library(plyr)
library(data.table)
library(genefilter)
library(Biobase)
library(lmQCM)
library(WGCNA)
library(GEOquery)
library(dplyr)
library(enrichR)
library(DT)
library(plotly)
library(openxlsx)
library(survival)
library(naturalsort)
# circos plot
library(circlize)
# library(topGO) # Somehow conflict with WGCNA and GEOquery. Deprecated.

options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB
options(shiny.sanitize.errors = FALSE)
source("utils.R")
# setwd("/Users/zhi/Desktop/GeneCoexpression/RGUI"); #mac #remove when deploy to shinyapps.io
# source("./lmQCM/GeneCoExpressionAnalysis.R")
# ----------------------------------------------------
data <- NULL
data_final <- NULL
GEO <- NULL
GSE_name_title <-NULL
finalExp <- NULL
finalSym <- NULL
finalSymChar <- NULL
text <- NULL
geneCharVector_global <- NULL
eigengene_matrix <- NULL
fname <- NULL # e.g. 12345_at
# listEnrichrDbs()
enrichr_dbs <- c("GO_Biological_Process_2017b",
                 "GO_Molecular_Function_2017b",
                 "GO_Cellular_Component_2017b",
                 "Jensen_DISEASES",
                 "Reactome_2016",
                 "KEGG_2016",
                 "Transcription_Factor_PPIs",
                 "Genome_Browser_PWMs",
                 "TRANSFAC_and_JASPAR_PWMs",
                 "ENCODE_TF_ChIP-seq_2015",
                 "Chromosome_Location",
                 "miRTarBase_2017",
                 "TargetScan_microRNA_2017",
                 "ChEA_2016")
enriched <- NULL # all enrichr results
final_genes_str <- NULL

function(input, output, session) {
  observeEvent(input$action1,{
      print('tab1')
      session$sendCustomMessage("myCallbackHandler", "tab1")
  })
  
  output$NCBI_GEO_Release_Date <- renderPlotly({
    
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    progress$set(message = "Loading NCBI GEO database. Last fetched version: 01/31/2018", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    # NCBI GEO
    load("./data/GEO_20180131.Rdata")
    GEO[["Actions"]] <- paste0('<div><button type="button" class="btn-analysis" id=dataset_analysis_',1:nrow(GEO),'>Analyze</button></div>')
    colnames(GEO)[which(names(GEO) == "Sample.Count")] <- "Samples"
    GEO <- GEO %>% select("Samples", everything())
    GEO <- GEO %>% select("Actions", everything())
    GEO <<- GEO
    # Basic process
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(GEO, extensions = 'Responsive', escape=F, selection = 'none', rownames = F)
    })
    
    # Other fancy processes
    # number of samples
    output$NCBI_GEO_Sample_Histogram <- renderPlotly({
      GEO_sample_number <- GEO$Samples
      plot_ly(x = GEO_sample_number, type = "histogram") %>% 
        layout(title = "Samples Histogram (log scale)",font=list(size = 10), 
               xaxis = list(title = "Number of Samples"),
               yaxis = list(title = "Occurence", type = "log")) %>% 
        layout(plot_bgcolor='transparent') %>% 
        layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)
    })
    
    #release date
    GEO_release_date <- as.Date(GEO$Release.Date,format = "%B %d, %Y")
    GEO_release_date_unique = cumsum(table(GEO_release_date)) #cummulative sum
    
    font <- list(family = "Courier New, monospace", size = 12, color = "black")
    x_axis <- list(title = "", titlefont = font, tickangle = 45, zeroline = TRUE)
    y_axis <- list(title = "Number of GSE data", titlefont = font)
    plot_ly(x = names(GEO_release_date_unique), y = as.numeric(GEO_release_date_unique),
            type = 'scatter', mode = 'lines') %>%
      layout(title = "Number of GSE Data Growing",font=list(size = 10), xaxis = x_axis, yaxis = y_axis) %>% 
      layout(plot_bgcolor='transparent') %>% 
      layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)
  })

observeEvent(input$dataset_lastClickId,{
  if (input$dataset_lastClickId%like%"dataset_analysis"){
  row_of_GEO <- as.numeric(gsub("dataset_analysis_","",input$dataset_lastClickId))
  myGSE <- GEO$Accession[row_of_GEO]
  GSE_name_title <<- GEO$Title[row_of_GEO]
  print(myGSE)
  smartModal(error=F, title = sprintf("Loading %s file from NCBI GEO Database", myGSE),
                         content = "We are currently loading your selected file from NCBI GEO Database ...")
  t <- try(gset <- getGEO(myGSE, GSEMatrix=TRUE, AnnotGPL=FALSE)) #AnnotGPL default is FALSE
  if("try-error" %in% class(t)) {
    removeModal()
    print("HTTP error 404")
    smartModal(error=T, title = "HTTP error 404", content = sprintf("%s is not available in NCBI GEO Database, please try another! (Hint: Maybe bad Internet connection)",myGSE))
    return()
  }

  if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  edata <- exprs(gset) #This is the expression matrix
  if (dim(edata)[1] == 0){
    removeModal()
    print("No expression data")
    smartModal(error=T, title = "Important message", content = sprintf("%s contains no expression data, please try another!", myGSE))
    return()
  }
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
    DT::datatable(data[1:ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1:ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                  extensions = 'Responsive', escape=F, selection = 'none')
  })
  output$mytable5 <- DT::renderDataTable({
    # Expression Value
    DT::datatable(data[ifelse(is.na(input$starting_row),1,input$starting_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), ifelse(is.na(input$starting_col),2,input$starting_col):ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                  extensions = 'Responsive', escape=F, selection = 'none')
  })
  output$mytable6 <- renderTable({
    # Gene ID
    data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1]
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
        smartModal(error=T, title = "Important message", content = "No file uploaded! Please retry!")
        return(NULL)
      }
      else {
        print("Reading file.")
        smartModal(error=F, title = "Intepreting uploaded file in progress", content = "Intepreting ...")
        fileExtension <- getFileNameExtension(input$csvfile$datapath)
        if(fileExtension == "csv"){
          data <<- read.csv(input$csvfile$datapath,
                            header = input$header,
                            sep = input$sep,
                            quote = input$quote)
          print("csv file Processed.")
        }
        else if(fileExtension == "txt"){
          data_temp = as.matrix(readLines(input$csvfile$datapath), sep = '\n')
          data_temp = strsplit(data_temp, split=input$sep)
          max.length <- max(sapply(data_temp, length))
          data_temp <- lapply(data_temp, function(v) { c(v, rep(NA, max.length-length(v)))})
          data_temp <- data.frame(do.call(rbind, data_temp))
          if(data_temp[dim(data_temp)[1],1] == "!series_matrix_table_end"){
            print("remove last row with \"!series_matrix_table_end\" ")
            data_temp = data_temp[-dim(data_temp)[1],]
          }
          data <<- data_temp
          print("txt file Processed.")
        } else if(fileExtension == "xlsx" || fileExtension == "xls"){
          data_temp <- read.xlsx(input$csvfile$datapath, sheet = 1, startRow = 1, colNames = TRUE)
          if(data_temp[dim(data_temp)[1],1] == "!series_matrix_table_end"){
            print("remove last row with \"!series_matrix_table_end\" ")
            data_temp = data_temp[-dim(data_temp)[1],]
          }
          data <<- data_temp
          print("xlsx / xls file Processed.")
        }
        removeModal()

        output$summary <- renderPrint({
          print(sprintf("Number of Genes: %d",dim(data)[1]))
          print(sprintf("Number of Samples: %d",(dim(data)[2]-1)))
          print("Annotation Platform: Unknown")
        })
        if ((dim(data)[2]-1) == 0){
          print("Number of samples 0!")
          smartModal(error=T, title = "Important message", content = sprintf("Target file contains no sample. This problem could because user pick a not matched separator (default: Comma), please try another seperator (e.g. Tab or Space)."))
          return()
        }
        output$mytable4 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[1:ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1:ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                        , extensions = 'Responsive', escape=F, selection = 'none')
        })
        output$mytable5 <- DT::renderDataTable({
          # Expression Value
          DT::datatable(data[ifelse(is.na(input$starting_row),1,input$starting_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), ifelse(is.na(input$starting_col),2,input$starting_col):ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                        , extensions = 'Responsive', escape=F, selection = 'none')
        })
        output$mytable6 <- renderTable({
          # Gene ID
          data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1]
        })
        print('tab2')
        session$sendCustomMessage("myCallbackHandler", "tab2")
      }
  })


  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                       Platform Convert
  #   |
  #   |
  #   +--------------------------------

  observeEvent(input$action_platform,{
    if(is.null(data)){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
    smartModal(error=F, title = "Converting...", content = "We are currently converting probe ID to Gene Symbol...")
    platform_name <- gsub(" ", "", input$platform_text, fixed = TRUE)

    print(sprintf("Platform: %s",platform_name))
    t <- try(gpl <- getGEO(platform_name))
    if("try-error" %in% class(t)) {
      removeModal()
      print("HTTP error 404")
      smartModal(error=c(T,F), title = "HTTP error 404", content = sprintf("Platform %s is not available in NCBI GEO Database, please try another!", platform_name))
      return()
    }
    print("Platform Loaded.")

    #https://www.rdocumentation.org/packages/GEOquery/versions/2.38.4/topics/GEOData-class
    gpltable <- Table(gpl)
    if (is.null(fname)){
      # data is not from GEO
      print("data is self-uploaded, so no fname defined.")
      fname <- data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):dim(data)[1], 1]# Gene ID
      fname <- gsub("\"","",fname) # convert "\"1553418_a_at\"" to "1553418_a_at"
      # save(fname,file="/Users/zhi/Desktop/fname.Rdata")
    }
    fname2 <- fname
    if (!is.null(gpltable$`Gene Symbol`)){
      print("load GPL table with name \"Gene Symbol\"")
      fname2 <- gpltable$`Gene Symbol`[match(fname,gpltable$ID)]
    }
    
    else if (!is.null(gpltable$`GENE_SYMBOL`)){
      print("load GPL table with name \"GENE_SYMBOL\"")
      fname2 <- gpltable$`GENE_SYMBOL`[match(fname,gpltable$ID)]
    }
    else {
        removeModal()
        print("HTTP error 404")
        smartModal(error=T, title = "Important message", content = sprintf("Error occured while loading %s. This issue could be different name defined on Gene Symbol (GENE_SYMBOL or others) in the platform.", platform_name))
        return()
    }

    print(dim(data))
    print(length(fname2))
    data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):dim(data)[1],1] <<- fname2
    # row.names(data) <- seq(1, length(fname2))
    output$mytable4 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data[1:ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1:ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                    extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable5 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data[ifelse(is.na(input$starting_row),1,input$starting_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), ifelse(is.na(input$starting_col),2,input$starting_col):ifelse(is.na(input$quicklook_col),10,input$quicklook_col)],
                    extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable6 <- renderTable({
      # Gene ID
      data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):ifelse(is.na(input$quicklook_row),100,input$quicklook_row), 1]
    })
    print("Platform Conversion finished")

    removeModal()
  })

  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                      Cleaning the Data
  #   |
  #   |
  #   +--------------------------------

  observeEvent(input$action3,{
      if(is.null(data)){
        smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
        return()
      }
      smartModal(error=F, title = "Preprocessing input data", content = "Preprocessing ...")
      print(dim(data))
      # Step 0
      RNA <- as.matrix(data[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data)[1], ifelse(is.na(input$starting_col),2,input$starting_col):dim(data)[2]])
      class(RNA) <- "numeric"
      geneID <- data.frame(data[ifelse(is.na(input$starting_gene_row),1,input$starting_gene_row):dim(data)[1], 1])
      print(dim(RNA))
      print(dim(geneID))
      
      # convert na to 0
      if (input$checkbox_NA){RNA[is.na(RNA)] <- 0}

      # Remove data with lowest 20% absolute exp value shared by all samples
      percentile <- ifelse(is.na(input$absolute_expval),0,input$absolute_expval)/100.
      print(sprintf("percentile 1: %f",percentile))
      # save(RNA, file="/Users/zhi/Desktop/RNA.Rdata")
      if (percentile > 0){
        RNA_filtered1 = RNA[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
        geneID_filtered1 = geneID[apply(RNA,1,max) > quantile(RNA, percentile)[[1]], ]
      }
      else {
        RNA_filtered1 = RNA
        geneID_filtered1 = geneID
      }
      
      print("after remove lowest k% abs exp value:")
      print(dim(RNA_filtered1))
      # Remove data with lowest 10% variance across samples
      percentile <- ifelse(is.na(input$variance_expval),0,input$variance_expval)/100.
      print(sprintf("percentile 2: %f",percentile))
      if (percentile > 0){
        index <- varFilter2(eset = RNA_filtered1, var.cutoff = percentile)
        RNA_filtered2 = RNA_filtered1[index, ]
        geneID_filtered2 = geneID_filtered1[index]
      }
      else {
        RNA_filtered2 = RNA_filtered1
        geneID_filtered2 = geneID_filtered1
      }
      
      print("after remove lowest l% var exp value:")
      print(dim(RNA_filtered2))

      # expData <- RNA_filtered2
      # res <- highExpressionProbes(geneID_filtered2, geneID_filtered2, expData)
      # ind1 <- res$first
      # uniGene <- as.character(res$second)
      # tmpExp <- expData[ind1,]
      uniGene <- geneID_filtered2
      tmpExp <- RNA_filtered2

      if (input$checkbox_logarithm){
        tmpExp[tmpExp <= 0] <- 0.000001
        tmpExp <- log(tmpExp)
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
      topN <- min(ifelse(is.na(input$max_gene_retain),Inf,input$max_gene_retain), nrow(tmpExp))
      #Remove gene symbol after vertical line: from ABC|123 to ABC:
      uniGene <- gsub("\\|.*$","", uniGene)

      finalExp <- tmpExp[sortInd[1:topN], ]
      finalSym <- uniGene[sortInd[1:topN]]
      finalSymChar <- as.character(finalSym)
      
      output$summary_advanced <- renderPrint({
        print(sprintf("Number of Genes: %d",dim(finalExp)[1]))
        print(sprintf("Number of Samples: %d",dim(finalExp)[2]))
      })
      
      # advanced
      if (input$sorting_adv_checkbox){
        finalPValue <- matrix(0, ncol = 0, nrow = length(finalSym))
        OS_IND <- as.numeric(data[input$row_osefs_ind, input$sort_col_start:dim(data)[2]])
        OS <- as.numeric(data[input$row_osefs, input$sort_col_start:dim(data)[2]])
        # print(OS_IND)
        # print(OS)
        # pvalue
        for(i in 1:dim(finalExp)[1]){
          rr = finalExp[i,]
          rr_sorted_list = sort(rr, decreasing = FALSE, index.return=T)
          rr_sorted = rr_sorted_list$x
          rr_sorted_idx = rr_sorted_list$ix
          medianB <- rep(0, length(rr))
          medianB[ rr_sorted_idx > length(rr)/2 ] = 1
          ss <- survdiff(Surv(OS, OS_IND) ~ medianB)
          finalPValue[i] <- 1-pchisq(ss$chisq, 1)
        }
        # print("final P value:")
        # print(finalPValue)
        finalPValue <- as.numeric(finalPValue)
        save(finalPValue,file="~/Desktop/finalPValue.Rdata")
        if (input$select_pval_adv_checkbox){
          finalExp <- finalExp[finalPValue<=input$advance_selection_pvalue,]
          finalSym <- finalSym[finalPValue<=input$advance_selection_pvalue]
          finalSymChar <- finalSymChar[finalPValue<=input$advance_selection_pvalue]
          finalPValue <- finalPValue[finalPValue<=input$advance_selection_pvalue]
        }
        
        data_final <<- data.frame(cbind(finalSym,finalPValue,finalExp))
        colnames(data_final)[1:2] = c("Gene_Symbol", sprintf("P-value of %s",input$choose_OS_EFS))
      }
      else{
        data_final <<- data.frame(cbind(finalSym,finalExp))
        colnames(data_final)[1] = "Gene_Symbol"
      }
      #finally no matter if just basic or advanced:
      finalExp <<- finalExp
      finalSym <<- finalSym
      finalSymChar <<- finalSymChar
      output$mytable_finaldata <- DT::renderDataTable({
        DT::datatable(data_final,
                      extensions = 'Responsive', escape=F, selection = 'none')
      })
      
      removeModal()
      print('tab3')
      session$sendCustomMessage("download_finaldata_ready","-")
      session$sendCustomMessage("myCallbackHandler", "tab3")
  })
  

  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                         l m Q C M
  #   |
  #   |
  #   +--------------------------------

  observeEvent(input$action4_lmQCM,{
    if(is.null(finalExp)){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
      #lmQCM
      smartModal(error=F, title = "Using lmQCM to calculate merged clusters", content = "Calculating. This could be several seconds to several minutes depends on number of genes. Please be patient.")
      mergedCluster <- lmQCM(finalExp, input$gamma, input$t, input$lambda, input$beta, input$minClusterSize, input$massiveCC)
      geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
      temp_eigengene <- matrix(0, nrow = length(mergedCluster), ncol = dim(finalExp)[2]) # Clusters * Samples

      temptext <- ""
      for (i in 1:(length(mergedCluster))) {
        vector <- as.matrix(mergedCluster[[i]])
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
      geneCharVector_global <<- geneCharVector
      text <<- temptext
      eigengene_matrix <<- temp_eigengene

      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))

      output$clusterResult <- DT::renderDataTable({
        geneCharVector2[["Actions"]] <- paste0('<div style="text-align: center"><button type="button" class="btn-analysis" id=go_analysis_',1:nrow(geneCharVector2),'>GO</button></div>')
        geneCharVector2[["Plots"]] <- paste0('<div style="text-align: center"><button type="button" class="btn-plot" id=circos_',1:nrow(geneCharVector2),'>Circos</button></div>')
        geneCharVector2 <- geneCharVector2 %>%
          select(c("Actions","Plots"), everything())
        colnames(geneCharVector2)[3:4] <- c("Cluster ID", "Genes")
        # print(head(geneCharVector2))
        DT::datatable(geneCharVector2, selection="none", escape=FALSE,
                      options = list(paging = F, searching = F, dom='t',ordering=T),
                      rownames = F#, colnames = NULL
        )
      })

      output$mytable7 <- renderTable({
        return(eigengene_matrix)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)

      removeModal()
      print('tab4')
      session$sendCustomMessage("download_cluster_ready","-")
      session$sendCustomMessage("myCallbackHandler", "tab4")
  })


  #=================================================
  #===================           ===================
  #==================  W G C N A  ==================
  #===================           ===================
  #=================================================

  observeEvent(input$checkPower, {
    if(is.null(finalExp)){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
    if (length(finalSym) > 0){
      #WGCNA
      smartModal(error=F, title = "Preview the power", content = "Running ...")

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
      if(is.null(finalExp)){
        smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
        return()
      }
      #WGCNA
      smartModal <- function(error=F, title = "Using WGCNA to calculate merged clusters", content = "Calculating. This could take a while depend on number of genes. Please be patient.")

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
      geneCharVector_withoutColor <- matrix(0, nrow = 0, ncol = length(unique(netcolors))-1)
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
        geneChar <- c(toString(i), labels2colors(i), matrixdata[geneID, 1])
        geneCharVector[i] <- list(geneChar)
        geneCharVector_withoutColor[i] <- list(c(toString(i), matrixdata[geneID, 1]))
        temptext <- paste(temptext, capture.output(cat(geneChar, sep=',')), sep="\n")
      }
      temptext <- substring(temptext, 2) # remove first \n separater
      geneCharVector_global <<- geneCharVector_withoutColor # remove the color label
      text <<- temptext
      eigengene_matrix <<- temp_eigengene

      ## Compute maximum length
      max.length <- max(sapply(geneCharVector, length))
      ## Add NA values to list elements
      geneCharVector2 <- lapply(geneCharVector, function(v) { c(v, rep(NA, max.length-length(v)))})
      ## Rbind
      # save(geneCharVector2,file = "/Users/zhi/Desktop/geneCharVector2.Rdata")
      geneCharVector2 <- data.frame(do.call(rbind, geneCharVector2))

      # output$clusterResult <- renderTable({
      #   return(geneCharVector2)
      # },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)

      output$clusterResult <- DT::renderDataTable({
        geneCharVector2[["Actions"]] <- paste0('<div style="text-align: center"><button type="button" class="btn-analysis" id=go_analysis_',1:nrow(geneCharVector2),'>GO</button></div>')
        geneCharVector2[["Plots"]] <- paste0('<div style="text-align: center"><button type="button" class="btn-plot" id=circos_',1:nrow(geneCharVector2),'>Circos</button></div>')
        geneCharVector2 <- geneCharVector2 %>%
          select(c("Actions","Plots"), everything())
        colnames(geneCharVector2)[3:5] <- c("Cluster ID", "Color Name", "Genes")
        # print(head(geneCharVector2))
        DT::datatable(geneCharVector2, selection="none", escape=FALSE,
                      options = list(paging = F, searching = F, dom='t',ordering=F),
                      rownames = F#, colnames = NULL
        )
      })

      output$mytable7 <- renderTable({
        return(eigengene_matrix)
      },rownames = FALSE, colnames = FALSE, na = "", bordered = TRUE)

      print('tab4')
      session$sendCustomMessage("download_cluster_ready","-")
      session$sendCustomMessage("myCallbackHandler", "tab4")
  })

  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                  GO Enrichment Analysis
  #   |
  #   |
  #   +--------------------------------
  observeEvent(input$action_finaldata4enrichr,{
    if(is.null(finalSym)){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
    genes_str = finalSym
    genes_str <- unlist(strsplit(genes_str, " /// "))
    # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
    print("genes_str for enrich analysis: ")
    print(genes_str)
    final_genes_str <<- genes_str
    updateTextAreaInput(session, "textareainput_GOEA",
                        label = paste(sprintf("Number of Genes: %d", length(final_genes_str)), input$controller),
                        value = paste(paste(final_genes_str, collapse = '\n'), input$controller))
    
    enriched <<- enrichr(final_genes_str, enrichr_dbs)
    
    Map(function(id) {
      dbres = enriched[[enrichr_dbs[id]]]
      dbres = dbres[ , -which(names(dbres) %in% c("Old.P.value","	Old.Adjusted.P.value"))]
      output[[paste("mytable_Enrichr",id,sep="_")]] <- DT::renderDataTable({DT::datatable(dbres, selection="none", escape=FALSE, pageLength = 100,
                                                                                          options = list(paging = F, searching = T, dom='t',ordering=T), extensions = 'Responsive',
                                                                                          rownames = T) #%>% formatRound(colnames(dbres)[3:dim(dbres)[2]], digits=8)
      })
    }, 1:length(enrichr_dbs))
    removeModal()
    print('tab5')
    session$sendCustomMessage("download_go_ready","-")
    session$sendCustomMessage("myCallbackHandler", "tab5")
  })
    
  observeEvent(input$button_lastClickId,{
    print("lastClickedId received.")
    if (input$button_lastClickId%like%"go_analysis"){
      cluster <- as.numeric(gsub("go_analysis_","",input$button_lastClickId))

      smartModal(error=F, title = "Performing GO Enrichment Analysis...", content = "Loading GO analysis results. This could take a while depends on number of genes. Please be patient.")
      print("row_of_final_cluster:")
      print(cluster)
      # print(geneCharVector_global[[cluster]])
      genes_str <- geneCharVector_global[[cluster]]
      genes_str <- unlist(strsplit(genes_str, " /// "))
      # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
      print("genes_str for enrich analysis: ")
      print(genes_str[-1])
      final_genes_str <<- genes_str[-1]
      updateTextAreaInput(session, "textareainput_GOEA",
                          label = paste(sprintf("Number of Genes: %d", length(final_genes_str)), input$controller),
                          value = paste(paste(final_genes_str, collapse = '\n'), input$controller))

      enriched <<- enrichr(final_genes_str, enrichr_dbs)

      Map(function(id) {
        dbres = enriched[[enrichr_dbs[id]]]
        dbres = dbres[ , -which(names(dbres) %in% c("Old.P.value","	Old.Adjusted.P.value"))]
        output[[paste("mytable_Enrichr",id,sep="_")]] <- DT::renderDataTable({DT::datatable(dbres, selection="none", escape=FALSE,
                                                                                            options = list(paging = T, pageLength = 100, searching = T, dom='t',ordering=T), extensions = 'Responsive',
                                                                                            rownames = T) #%>% formatRound(colnames(dbres)[3:dim(dbres)[2]], digits=8)
        })
      }, 1:length(enrichr_dbs))
      removeModal()
      print('tab5')
      session$sendCustomMessage("download_go_ready","-")
      session$sendCustomMessage("myCallbackHandler", "tab5")
    }
    #   +------------------------------------------------------------+
    #   |
    #   |
    #   |                        Circos Plot
    #   |
    #   |
    #   +--------------------------------
    
    
    if (input$button_lastClickId%like%"circos"){
      smartModal(error=F, title = "Processing", content = "We are working on your customized circos plot ...")
      cluster <- as.numeric(gsub("circos_","",input$button_lastClickId))
      genes_str <- geneCharVector_global[[cluster]]
      genes_str <- unlist(strsplit(genes_str, " /// "))
      # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
      print("genes_str for circos plot: ")
      print(genes_str[-1])
      
      # import hg19 and hg38
      load("./data/UCSC_hg19_refGene_20180330.Rdata") # varname: hg19
      load("./data/UCSC_hg38_refGene_20180330.Rdata") # varname: hg38
      # genes_str <- c("LOC102725121", "FAM138A", "RIMS2", "LINC01128", "MMP23A", "ULK4P1")
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
      
      output$circos_plot_component_hg38 <- renderPlot({
        factors_count = as.data.frame(hg38.ring.lengthsum)
        factors = factor(factors_count[,1], levels = factors_count[,1])
        xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
        rownames(xlim) = factors_count[,1]
        BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
        circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genomics (GRCh38/hg38)")
      })
      
      output$circos_plot_component_hg19 <- renderPlot({
        factors_count = as.data.frame(hg19.ring.lengthsum)
        factors = factor(factors_count[,1], levels = factors_count[,1])
        xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
        rownames(xlim) = factors_count[,1]
        BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
        circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genomics (GRCh37/hg19)")
      })
      removeModal()
      print('tab4_functional_plots')
      session$sendCustomMessage("myCallbackHandler", "tab4_functional_plots")
    }
    
  })
  
  observeEvent(input$circos_button_update_1,{
    smartModal(error=F, title = "Processing", content = "We are working on your customized circos plot ...")
    cluster <- as.numeric(gsub("circos_","",input$button_lastClickId))
    genes_str <- strsplit(input$textareainput_circos, "\n")[[1]]
    genes_str <- unlist(strsplit(genes_str, " /// "))
    print("genes_str for circos plot: ")
    print(genes_str)
    
    # import hg19 and hg38
    load("./data/UCSC_hg19_refGene_20180330.Rdata") # varname: hg19
    load("./data/UCSC_hg38_refGene_20180330.Rdata") # varname: hg38
    # genes_str <- c("LOC102725121", "FAM138A", "RIMS2", "LINC01128", "MMP23A", "ULK4P1")
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
    
    output$circos_plot_component_hg38 <- renderPlot({
      factors_count = as.data.frame(hg38.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genomics (GRCh38/hg38)")
    })
    
    output$circos_plot_component_hg19 <- renderPlot({
      factors_count = as.data.frame(hg19.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genomics (GRCh37/hg19)")
    })
    removeModal()
  })
  
  
  
  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                     Download Handler
  #   |
  #   |
  #   +--------------------------------

  output$downloadData1 <- downloadHandler(
    filename = function() {
      name = "mergedCluster"
      # name = paste("lmQCMresult","gamma",gamma(),"lambda",lambda(),"t",t(),"beta",beta(),"minClusterSize",minClusterSize(), sep = "_", collapse = NULL)
      paste(name, input$filetype1, sep = ".")
    },
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
    filename = function() {
      name = "EigengeneMatrix"
      paste(name, input$filetype2, sep = ".")
    },
    content = function(file) {
      separator <- switch(input$filetype2, "csv" = ',', "txt" = '\t')
      write.table(eigengene_matrix, file = file, append = FALSE, quote = TRUE, sep = separator,
                  eol = "\r\n", na = "NA", dec = ".", row.names = F,
                  col.names = F, qmethod = c("escape", "double"),
                  fileEncoding = "")
    }
  )

  output$download_finaldata <- downloadHandler(
    filename = function() {
      name = "finaldata.csv"
    },
    content = function(file) {
      write.table(data_final, file = file, append = FALSE, quote = TRUE, sep = ',',
                  eol = "\r\n", na = "NA", dec = ".", row.names = F,
                  col.names = T, qmethod = c("escape", "double"),
                  fileEncoding = "")
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename = 'GO_results.zip',
    content = function(fname) {
      separator <- switch(input$filetype3, "csv" = ',', "txt" = '\t')
      if(separator == ','){
        fs <- paste0(enrichr_dbs, '.csv')
      }
      else{
        fs <- paste0(enrichr_dbs, '.txt')
      }
      fs <- c(fs, 'genes_list.txt')
      for(i in 1:length(enrichr_dbs)){
        write.table(enriched[[enrichr_dbs[i]]], file = fs[i], sep = separator)
        print(fs[i])
      }
      if(length(final_genes_str) > 0){
        write(final_genes_str, file = 'genes_list.txt', sep = "\n")
      }
      else{
        write("You cannot believe I disabled this function. Haha.", file = 'genes_list.txt')
      }
      zip(zipfile=fname, files=fs)
      if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}

    },
    contentType = "application/zip"
  )
}