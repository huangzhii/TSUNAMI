options(repos = BiocManager::repositories())
getOption("repos")
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
library(shinyWidgets)
library(progress)
# circos plot
library(circlize)
# library(topGO) # Somehow conflict with WGCNA and GEOquery. Deprecated.

options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB
options(shiny.sanitize.errors = FALSE)
options(stringAsFactors = FALSE)
source("utils.R")
# setwd("/Users/zhi/Desktop/GeneCoexpression/RGUI"); #mac #remove when deploy to shinyapps.io
# source("./lmQCM/GeneCoExpressionAnalysis.R")
# ----------------------------------------------------
setClass("QCMObject", representation(clusters.id = "list", clusters.names = "list",
                                     eigengene.matrix = "data.frame"))
# listEnrichrDbs()
path2TCGAdata <<- "http://web.ics.purdue.edu/~huang898/TSUNAMI_data/20190303_TCGA_gdac_mRNAseq/"
enrichr_dbs <<- c("GO_Biological_Process_2018",
                 "GO_Molecular_Function_2018",
                 "GO_Cellular_Component_2018",
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
TCGA.gdac.list <<- read.csv("./data/TCGA_gdac_metadata.csv", header = T)

server <- function(input, output, session) {
  
  data <- reactiveVal(0)
  data_final <- reactiveVal(0)
  GEO <- reactiveVal(0)
  GSE_name_title <- reactiveVal(0)
  mygset <- reactiveVal(0)
  finalExp <- reactiveVal(0)
  finalSym <- reactiveVal(0)
  sampleID <- reactiveVal(0)
  finalSymChar <- reactiveVal(0)
  text.final <- reactiveVal(0)
  geneCharVector_global <- reactiveVal(0)
  eigengene_matrix <- reactiveVal(0)
  fname <- reactiveVal(0) # e.g. 12345_at
  enriched <- reactiveVal(0) # all enrichr results
  final_genes_str <- reactiveVal(0)
  
  observeEvent(input$action1,{
      print('tab1')
      session$sendCustomMessage("myCallbackHandler", "tab1")
  })
  
  output$NCBI_GEO_Release_Date <- renderPlotly({
    
    # Create a Progress object
    progress <- shiny::Progress$new(session)
    progress$set(message = "Loading NCBI GEO database. Last fetched version: 03/03/2019", value = 100)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    # NCBI GEO
    load("./data/GEO_20190303.Rdata")
    GEO.temp <- GEO
    GEO.temp[["Actions"]] <- paste0('<div><button type="button" class="btn-analysis" id=GEO_dataset_analysis_',1:nrow(GEO.temp),'>Analyze</button></div>')
    colnames(GEO.temp)[which(names(GEO.temp) == "Sample.Count")] <- "Samples"
    GEO.temp <- GEO.temp %>% select("Samples", everything())
    GEO.temp <- GEO.temp %>% select("Actions", everything())
    GEO(GEO.temp)
    # Basic process
    TCGA.gdac.list[["Actions"]] <- paste0('<div><button type="button" class="btn-analysis" id=TCGA_dataset_analysis_',1:nrow(TCGA.gdac.list),'>Analyze</button></div>')
    TCGA.gdac.list <- TCGA.gdac.list %>% select("Actions", everything())
    output$mytable0 <- DT::renderDataTable({
      DT::datatable(TCGA.gdac.list, extensions = 'Responsive', escape=F, selection = 'none',
                    options = list(paging = F, searching = T, dom='t',ordering=T), rownames = F)
    })
    output$mytable1 <- DT::renderDataTable({
      DT::datatable(GEO(), extensions = 'Responsive', escape=F, selection = 'none', rownames = F)
    })
    
    # Other fancy processes
    # number of samples
    output$NCBI_GEO_Sample_Histogram <- renderPlotly({
      GEO_sample_number <- GEO()$Samples
      plot_ly(x = GEO_sample_number, type = "histogram") %>% 
        layout(title = "Samples Histogram (log scale)",font=list(size = 10), 
               xaxis = list(title = "Number of Samples"),
               yaxis = list(title = "Occurence", type = "log")) %>% 
        layout(plot_bgcolor='transparent') %>% 
        layout(paper_bgcolor='transparent') %>% config(displayModeBar = F)
    })
    
    #release date
    GEO_release_date <- as.Date(GEO()$Release.Date,format = "%B %d, %Y")
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


observeEvent(input$dataset_lastClickId_mytable0,{
  if (input$dataset_lastClickId_mytable0%like%"TCGA_dataset_analysis"){
    row_of_TCGA <- as.numeric(gsub("TCGA_dataset_analysis_","",input$dataset_lastClickId_mytable0))
    print(TCGA.gdac.list[row_of_TCGA,])
    cancer.name = TCGA.gdac.list[row_of_TCGA,2]
    smartModal(error=F, title = sprintf("Loading TCGA %s Data", cancer.name),
               content = sprintf("Loading TCGA %s Data ...", cancer.name))
    load(url(paste0(path2TCGAdata, cancer.name, ".Rdata")))
    TCGA.mRNAseq = data.frame(cbind(rownames(TCGA.mRNAseq), TCGA.mRNAseq))
    colnames(TCGA.mRNAseq)[1] = "Gene"
    data(TCGA.mRNAseq)
    updateTextInput(session, "platform_text", value = "")
    output$summary <- renderPrint({
      print(sprintf("Number of Genes: %d",dim(data())[1]))
      print(sprintf("Number of Samples: %d",dim(data())[2]))
    })
    
    print("GEO file is downloaded to server and processed.")
    output$mytable4 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data(),extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable5 <- DT::renderDataTable({
      # Expression Value
      verified_data = data()[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data())[1],
                             c(1, ifelse(is.na(input$starting_col),2,input$starting_col):dim(data())[2])]
      colnames(verified_data)[1] <- "Gene"
      DT::datatable(verified_data,extensions = 'Responsive', escape=F, selection = 'none')
    })
    
    print('tab2')
    removeModal()
    session$sendCustomMessage("myCallbackHandler", "tab2")
  }
})


observeEvent(input$dataset_lastClickId_mytable1,{
  if (input$dataset_lastClickId_mytable1%like%"GEO_dataset_analysis"){
    row_of_GEO <- as.numeric(gsub("GEO_dataset_analysis_","",input$dataset_lastClickId_mytable1))
    myGSE <- GEO()$Accession[row_of_GEO]
    GSE_name_title(GEO()$Title[row_of_GEO])
    message(GSE_name_title())
    print(myGSE)
    smartModal(error=F, title = sprintf("Loading %s file from NCBI GEO Database", myGSE),
               content = "We are currently loading your selected file from NCBI GEO Database ...")
    t <- try(gset <- getGEO(myGSE, GSEMatrix=TRUE, AnnotGPL=FALSE)) #AnnotGPL default is FALSE
    if("try-error" %in% class(t)) {
      removeModal()
      print("HTTP error 404")
      smartModal(error=T, title = "HTTP error 404", content = sprintf("%s is not available in NCBI GEO Database, please try other available GSE data (e.g., GSE17537, GSE73119). (Hint: Maybe bad Internet connection)",myGSE))
      return()
    }
    
    # if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
    if (length(gset) == 0){
      removeModal()
      print("This GSE accession doesn't contain any series matrix.")
      smartModal(error=T, title = "This GSE accession doesn't contain any series matrix.", content = sprintf("%s doesn't contain any series matrix. Please try other available GSE data (e.g., GSE17537, GSE73119).",myGSE))
      return()
    }
    mygset(gset[[1]])
    # select index
    if (length(gset) > 1){
      gset.names = names(gset)
      print(gset.names)
      gset.names = paste0(rep(1:length(gset.names)), ". ", gset.names)
      gset.names = paste0(gset.names, collapse = "\n")
      removeModal()
      inputSweetAlert(session, inputId = "whichgset", title = "This data contains multiple serie matrices.\nInput the ID (e.g. 1) to select desired one.",
                      text = gset.names,
                      type = "info", btn_labels = "Ok")
      observeEvent(input$whichgset, {
        mygset(gset[[as.numeric(input$whichgset)]])
      })
    }
  }
})


observeEvent(data(),{
  if (typeof(data()) == "double"){
    return()
  }
  samples = colnames(data())[input$starting_col:dim(data())[2]]
  output$data_sample_subgroup_ui <- renderUI({
    checkboxGroupInput("data_sample_subgroup", "Choose samples:",
                       choiceNames = samples,
                       choiceValues = samples,
                       selected = samples
    )
  })
})


observeEvent(mygset(),{
  if (typeof(mygset()) == "double"){
    return()
  }
  smartModal(error=F, title = sprintf("Selecting table from NCBI GEO Database"),
             content = "We are currently loading your selected file from NCBI GEO Database ...")
  edata <- exprs(mygset()) #This is the expression matrix
  if (dim(edata)[1] == 0 || is.null(dim(edata))){
    removeModal()
    print("No expression data")
    smartModal(error=T, title = "Important message", content = sprintf("%s doesn't contain any expression data, please try other available GSE data (e.g., GSE17537, GSE73119).", myGSE))
    return()
  }
  # pdata <- pData(mygset()) # data.frame of phenotypic information.
  fname(featureNames(mygset())) # e.g. 12345_at
  data(cbind(fname(), edata))
  data.temp <- data()
  row.names(data.temp) <- seq(1, length(fname()))
  data(data.temp)
  
  updateTextInput(session, "platform_text", value = paste(annotation(mygset()), input$controller))
  output$summary <- renderPrint({
    print(sprintf("Number of Genes: %d",dim(edata)[1]))
    print(sprintf("Number of Samples: %d",dim(edata)[2]))
    print(sprintf("Annotation Platform: %s",annotation(mygset())))
  })
  
  print("GEO file is downloaded to server and processed.")
  output$mytable4 <- DT::renderDataTable({
    # Expression Value
    DT::datatable(data(),extensions = 'Responsive', escape=F, selection = 'none')
  })
  output$mytable5 <- DT::renderDataTable({
    # Expression Value
    verified_data = data()[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data())[1],
                           c(1, ifelse(is.na(input$starting_col),2,input$starting_col):dim(data())[2])]
    colnames(verified_data)[1] <- "Gene"
    DT::datatable(verified_data,extensions = 'Responsive', escape=F, selection = 'none')
  })
  print('tab2')
  removeModal()
  session$sendCustomMessage("myCallbackHandler", "tab2")
})

output$readingcsv <- reactive({
  print("readingcsv = 1")
  return(is.null(input$csvfile))
})


observeEvent(input$action2,{
  options(stringsAsFactors = FALSE)
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
        data.temp <- read.csv(input$csvfile$datapath,
                          header = input$header,
                          sep = input$sep,
                          quote = input$quote)
        data(data.temp)
        print("csv file Processed.")
      }
      else if(fileExtension == "txt"){
        data_temp = as.matrix(readLines(input$csvfile$datapath), sep = '\n')
        data_temp = strsplit(data_temp, split=input$sep)
        max.length <- max(sapply(data_temp, length))
        data_temp <- lapply(data_temp, function(v) { c(v, rep(NA, max.length-length(v)))})
        data_temp <- data.frame(do.call(rbind, data_temp))
        colnames(data_temp) = data_temp[1,]
        data_temp = data_temp[2:dim(data_temp)[1],]
        if (is.null(dim(data_temp))){
          removeModal()
          sendSweetAlert(session, title = "Error", text = "Input file extention TXT found, but cannot construct matrix. Please confirm your data format and separator.", type = "error",
                         btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
          return()
        }
        if(data_temp[dim(data_temp)[1],1] == "!series_matrix_table_end"){
          print("remove last row with \"!series_matrix_table_end\" ")
          data_temp = data_temp[-dim(data_temp)[1],]
        }
        # data_temp <- print.data.frame(data.frame(data_temp), quote=FALSE)
        data(data_temp)
        print("txt file Processed.")
      } else if(fileExtension == "xlsx"){
        data_temp <- read.xlsx(input$csvfile$datapath, sheet = 1, startRow = 1, colNames = TRUE)
        if(data_temp[dim(data_temp)[1],1] == "!series_matrix_table_end"){
          print("remove last row with \"!series_matrix_table_end\" ")
          data_temp = data_temp[-dim(data_temp)[1],]
        }
        # data_temp <- print.data.frame(data.frame(data_temp), quote=FALSE)
        data(data_temp)
        print("xlsx file Processed.")
      } else if(fileExtension == "xls"){
        sendSweetAlert(session, title = "Error", text = "We discontinued to support XLS format. Please resubmit with another file format.", type = "error",
                       btn_labels = "Ok", html = FALSE, closeOnClickOutside = TRUE)
      }
      removeModal()
      
      # # first column as gene
      # data.temp = data()#[,2:dim(data())[2]]
      # rownames(data.temp) = data()[,1]
      # data(data.temp)

      output$summary <- renderPrint({
        print(sprintf("Number of Genes: %d",dim(data())[1]))
        print(sprintf("Number of Samples: %d",(dim(data())[2]-1)))
        print("Annotation Platform: Unknown")
      }, quoted = FALSE)
      if ((dim(data())[2]-1) == 0){
        print("Number of samples 0!")
        smartModal(error=T, title = "Important message", content = sprintf("Target file contains no sample. This problem could because user pick a not matched separator (default: Comma), please try another seperator (e.g. Tab or Space)."))
        return()
      }
      output$mytable4 <- DT::renderDataTable({
        # Expression Value
        DT::datatable(data(),extensions = 'Responsive', escape=F, selection = 'none')
      })
      output$mytable5 <- DT::renderDataTable({
        # Verified data
        verified_data = data()[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data())[1],
                             c(1, ifelse(is.na(input$starting_col),2,input$starting_col):dim(data())[2])]
        colnames(verified_data)[1] <- "Gene"
        DT::datatable(verified_data,extensions = 'Responsive', escape=F, selection = 'none')
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
    if(is.null(data())){
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
    if (is.null(fname())){
      # data is not from GEO
      print("data is self-uploaded, so no fname defined.")
      fname(data()[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data())[1], 1])# Gene ID
      fname.temp <- fname()
      fname.temp <- noquote(fname.temp) # convert "\"1553418_a_at\"" to "1553418_a_at" (safer)
      fname.temp <- gsub("\"","",fname.temp) # convert "\"1553418_a_at\"" to "1553418_a_at"
      fname(fname.temp)
      # save(fname,file="/Users/zhi/Desktop/fname.Rdata")
    }
    fname2 <- fname()
    if (!is.null(gpltable$`Gene Symbol`)){
      print("load GPL table with name \"Gene Symbol\"")
      fname2 <- gpltable$`Gene Symbol`[match(fname(), gpltable$ID)]
    }
    else if (!is.null(gpltable$`GENE_SYMBOL`)){
      print("load GPL table with name \"GENE_SYMBOL\"")
      fname2 <- gpltable$`GENE_SYMBOL`[match(fname(), gpltable$ID)]
    }
    else if (!is.null(gpltable$`CLONE_ID`)){
      print("load GPL table with name \"CLONE_ID\"")
      fname2 <- gpltable$`CLONE_ID`[match(fname(), gpltable$ID)]
    }
    else {
        removeModal()
        print("HTTP error 404")
        smartModal(error=T, title = "Important message", content = sprintf("Error occured while loading %s. This issue could be different name defined on Gene Symbol (GENE_SYMBOL or others) in the platform.", platform_name))
        return()
    }

    print(dim(data()))
    print(length(fname2))
    data.temp <- data()
    data.temp[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data.temp)[1],1] <- fname2
    rownames(data.temp)[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data.temp)[1]] <- fname2
    data(data.temp)
    # row.names(data) <- seq(1, length(fname2))
    output$mytable4 <- DT::renderDataTable({
      # Expression Value
      DT::datatable(data(),extensions = 'Responsive', escape=F, selection = 'none')
    })
    output$mytable5 <- DT::renderDataTable({
      # Expression Value
      
      verified_data = data()[ifelse(is.na(input$starting_row),1,input$starting_row):dim(data())[1],
                           c(1, ifelse(is.na(input$starting_col),2,input$starting_col):dim(data())[2])]
      colnames(verified_data)[1] <- "Gene"
      DT::datatable(verified_data,extensions = 'Responsive', escape=F, selection = 'none')
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

  observe({
    if(typeof(data()) != "double"){
      samples = colnames(data())[input$starting_col:dim(data())[2]]
      if(input$data_sample_subgroup_selectall == 0) return(NULL) 
      else if (input$data_sample_subgroup_selectall%%2 == 1)
      {
        updateCheckboxGroupInput(session,"data_sample_subgroup",choices=samples)
      }
      else
      {
        updateCheckboxGroupInput(session,"data_sample_subgroup",choices=samples, selected=samples)
      }
    }
  })

  observeEvent(input$action3,{
      if(is.null(data())){
        smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
        return()
      }
      smartModal(error=F, title = "Preprocessing input data", content = "Preprocessing ...")
      
      RNA = data()
      if(!is.null(input$data_sample_subgroup)){
        RNA = RNA[, colnames(RNA) %in% input$data_sample_subgroup]
      }
      print(dim(RNA))
      
      withProgress(message = 'Preprocessing input data', value = 0, {
      # Step 0
      # Increment the progress bar, and update the detail text.
      incProgress(1/5, detail = "Parsing Input Data")
        
      if (!is.null(input$starting_col) && input$starting_col >= 2){
        geneID <- data.frame(RNA[,input$starting_col-1])
      } else{ # use rownames as gene ID
        geneID <- data.frame(rownames(RNA)[ifelse(is.na(input$starting_row),1,input$starting_row):dim(RNA)[1]])
      }
      if(!is.null(input$data_sample_subgroup)){
        RNA <- as.matrix(RNA[ifelse(is.na(input$starting_row),1,input$starting_row):dim(RNA)[1],])
      } else{
        RNA <- as.matrix(RNA[ifelse(is.na(input$starting_row),1,input$starting_row):dim(RNA)[1], ifelse(is.na(input$starting_col),2,input$starting_col):dim(RNA)[2]])
      }
      class(RNA) <- "numeric"
      print(dim(RNA))
      print(dim(geneID))
      
      
      # convert na to 0
      if (input$checkbox_NA){RNA[is.na(RNA)] <- 0}

      # Remove data with lowest 20% mean exp value shared by all samples
      percentile <- ifelse(is.na(input$mean_expval),0,input$mean_expval)/100.
      percentile.mean <- percentile
      print(sprintf("percentile 1: %f",percentile))
      # save(RNA, file="~/Desktop/RNA.Rdata")
      # save(geneID, file="~/Desktop/geneID.Rdata")
      if (percentile > 0){
        RNAmean = apply(RNA,1,mean)
        RNA_filtered1 = RNA[RNAmean > quantile(RNAmean, percentile)[[1]], ]
        geneID_filtered1 = geneID[RNAmean > quantile(RNAmean, percentile)[[1]], ]
      }
      else {
        RNA_filtered1 = RNA
        geneID_filtered1 = as.matrix(geneID)
      }
      
      print("after remove lowest k% mean exp value:")
      print(dim(RNA_filtered1))
      incProgress(1/5, detail = "Remove lowest k% mean exp value")
      
      # Remove data with lowest 10% variance across samples
      percentile <- ifelse(is.na(input$variance_expval),0,input$variance_expval)/100.
      percentile.var <- percentile
      print(sprintf("percentile 2: %f",percentile))
      if (percentile > 0){
        if (dim(RNA_filtered1)[2] > 3){
          index <- varFilter2(eset = RNA_filtered1, var.cutoff = percentile)
          RNA_filtered2 = RNA_filtered1[index, ]
          geneID_filtered2 = geneID_filtered1[index]
        }
        else{
          smartModal(error=F, title = "Preprocessing input data", content = "Preprocessing ... \n Cannot calculate order statistic on object with less than 3 columns, will not remove data based on variance.")
          RNA_filtered2 = RNA_filtered1
          geneID_filtered2 = geneID_filtered1
        }
      }
      else {
        RNA_filtered2 = RNA_filtered1
        geneID_filtered2 = geneID_filtered1
      }
      
      print("after remove lowest l% var exp value:")
      print(dim(RNA_filtered2))
      incProgress(1/5, detail = "Remove lowest l% var exp value")

      # expData <- RNA_filtered2
      # res <- highExpressionProbes(geneID_filtered2, geneID_filtered2, expData)
      # ind1 <- res$first
      # uniGene <- as.character(res$second)
      # tmpExp <- expData[ind1,]
      uniGene <- geneID_filtered2
      tmpExp <- RNA_filtered2

      if (input$checkbox_logarithm){
        # tmpExp[tmpExp <= 0] <- 0.000001
        tmpExp <- log2(tmpExp+1)
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
      
      incProgress(1/5, detail = "Sort genes based on mean.")
      
      nSample <- ncol(tmpExp)
      res <- sort.int(rowMeans(tmpExp), decreasing = TRUE, index.return=TRUE)
      sortMean <- res$x
      sortInd <- res$ix
      # topN <- min(ifelse(is.na(input$max_gene_retain),Inf,input$max_gene_retain), nrow(tmpExp))
      #Remove gene symbol after vertical line: from ABC|123 to ABC:
      uniGene <- gsub("\\|.*$","", uniGene)

      finalExp(tmpExp[sortInd, ])
      finalSym(uniGene[sortInd])
      finalSymChar(as.character(finalSym()))
      
      # save(finalExp, file = "~/Desktop/finalExp.Rdata")
      # save(finalSym, file = "~/Desktop/finalSym.Rdata")
      # save(finalSymChar, file = "~/Desktop/finalSymChar.Rdata")
      
      incProgress(1/5, detail = "\nPost-processing...")
      
      data_final(data.frame(cbind(finalSym(), finalExp())))
      data_final.temp <- data_final()
      colnames(data_final.temp)[1] = "Gene_Symbol"
      data_final(data_final.temp)
      #finally no matter if just basic or advanced:
      sampleID(colnames(finalExp()))
      output$mytable_finaldata <- DT::renderDataTable({
        DT::datatable(data_final(),
                      extensions = 'Responsive', escape=F, selection = 'none')
      })
      removeModal()
      }) # progress bar finished.
      
      # showModal(modalDialog(
      #   title = "Data After Preprocessing", footer = modalButton("OK"), easyClose = TRUE,
      #   div(class = "busy",
      #       p(sprintf("%d genes ---- Original.", dim(RNA)[1])),
      #       p(sprintf("%d genes remained after remove lowest %.2f%% means.", dim(RNA_filtered1)[1], percentile.mean*100)),
      #       p(sprintf("%d genes remained after remove lowest %.2f%% variances.", dim(RNA_filtered2)[1], percentile.var*100)),
      #       style = "margin: auto; text-align: left"
      #   )
      # ))
      
      sendSweetAlert(session, title = "Data After Preprocessing",
                     text = sprintf("%d genes ---- Original.\n%d genes remained after remove lowest %.2f%% means.\n%d genes remained after remove lowest %.2f%% variances.",
                                    dim(RNA)[1], dim(RNA_filtered1)[1], percentile.mean*100, dim(RNA_filtered2)[1], percentile.var*100),
                     type = "success",
                     btn_labels = "OK", html = FALSE, closeOnClickOutside = TRUE)

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
    if(is.null(finalExp())){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
      #lmQCM
      smartModal(error=F, title = "lmQCM [1/2]", content = "Calculating massive correlation coefficient matrix...")

      data_in = finalExp()
      gamma = input$gamma
      t = input$t
      lambda = input$lambda
      beta = input$beta
      minClusterSize = input$minClusterSize
      CCmethod = tolower(input$massiveCC)
      normalization = input$lmQCM_weight_normalization
      
#----------------------------------------------------------------------------
# lmQCM start
#----------------------------------------------------------------------------
      withProgress(message = 'lmQCM [1/2]: ', value = 0, {
        incProgress(1/16, detail = "Calculating massive correlation coefficient matrix...")
        cMatrix <- cor(t(data_in), method = CCmethod)
        diag(cMatrix) <- 0
        incProgress(1/2, detail = "Find the local maximal edges")
        if(normalization){
          # Normalization
          cMatrix <- abs(cMatrix)
          D <- rowSums(cMatrix)
          D.half <- 1/sqrt(D)
          
          cMatrix <- apply(cMatrix, 2, function(x) x*D.half )
          cMatrix <- t(apply(cMatrix, 1, function(x) x*D.half ))
        }
  
        C <- list()
        nRow <- nrow(cMatrix)
        maxV <- apply(cMatrix, 2, max)
        maxInd <- apply(cMatrix, 2, which.max) # several diferrences comparing with Matlab results
        
        # Step 1 - find the local maximal edges
        lm.ind <- which(maxV == sapply(maxInd, function(x) max(cMatrix[x,])))
        maxEdges <- cbind(maxInd[lm.ind], lm.ind)
        maxW <- maxV[lm.ind]
        
        res <- sort.int(maxW, decreasing = TRUE, index.return=TRUE)
        sortMaxV <- res$x
        sortMaxInd <- res$ix
        sortMaxEdges <- maxEdges[sortMaxInd,]
        message(sprintf("Number of Maximum Edges: %d", length(sortMaxInd)))
        
        incProgress(7/16, detail = sprintf("Number of Maximum Edges: %d", length(sortMaxInd)))
        
        currentInit <- 1
        noNewInit <- 0
        
        nodesInCluster <- matrix(0, nrow = 0, ncol = 1)
        
      })
      removeModal()
      
      
      pb <- progress_bar$new(format = " Calculating [:bar] :percent eta: :eta",
                             total = length(sortMaxInd), clear = F, width=60)
      iter = 0
      
      
      progressSweetAlert(
        session = session, id = "myprogress",
        title = "lmQCM [2/2] Merging ...",
        display_pct = TRUE, value = 0
      )
      
      while ((currentInit <= length(sortMaxInd)) & (noNewInit == 0)) {
        pb$tick()
        iter = iter + 1
        updateProgressBar(
          session = session,
          id = "myprogress",
          value = iter/length(sortMaxInd)*100
        )
        if (sortMaxV[currentInit] < (gamma * sortMaxV[1]) ) {
          noNewInit <- 1
        }
        else {
          if ( (is.element(sortMaxEdges[currentInit, 1], nodesInCluster) == FALSE) & is.element(sortMaxEdges[currentInit, 2], nodesInCluster) == FALSE) {
            newCluster <- sortMaxEdges[currentInit, ]
            addingMode <- 1
            currentDensity <- sortMaxV[currentInit]
            nCp <- 2
            totalInd <- 1:nRow
            remainInd <- setdiff(totalInd, newCluster)
            # C = setdiff(A,B) for vectors A and B, returns the values in A that
            # are not in B with no repetitions. C will be sorted.
            while (addingMode == 1) {
              neighborWeights <- colSums(cMatrix[newCluster, remainInd])
              maxNeighborWeight <- max(neighborWeights)
              maxNeighborInd <- which.max(neighborWeights)
              c_v = maxNeighborWeight/nCp;
              alphaN = 1 - 1/(2*lambda*(nCp+t));
              if (c_v >= alphaN * currentDensity) {
                newCluster <- c(newCluster, remainInd[maxNeighborInd])
                nCp <- nCp + 1
                currentDensity <- (currentDensity*((nCp-1)*(nCp-2)/2)+maxNeighborWeight)/(nCp*(nCp-1)/2)
                remainInd <- setdiff(remainInd, remainInd[maxNeighborInd]);
              }
              else {
                addingMode <- 0
              }
            }
            nodesInCluster <- c(nodesInCluster, newCluster)
            C <- c(C, list(newCluster))
          }
        }
        currentInit <- currentInit + 1
      }
      
      if(length(C) == 0) {
        removeModal()
        sendSweetAlert(
          session = session,
          title ="Clusters size = 0 before merging. Please try other set of parameters. Program stopped.",
          type = "error"
        )
        return()
      }
      
      closeSweetAlert(session = session)
      sendSweetAlert(
        session = session,
        title ="lmQCM Calculation completed !",
        type = "success"
      )
      
      t <- try(clusters <- merging_lmQCM(C, beta, minClusterSize))
      if("try-error" %in% class(t)) {
        removeModal()
        sendSweetAlert(
          session = session,
          title ="Too few genes to perform lmQCM clustering and merging.",
          type = "error"
        )
        clusters <- C
      }
      
      # map rownames to clusters
      clusters.names = list()
      for (i in 1:length(clusters)){
        mc = clusters[[i]]
        clusters.names[[i]] = rownames(data_in)[mc]
      }
      # calculate eigengene
      eigengene.matrix <- matrix(0, nrow = length(clusters), ncol = dim(data_in)[2]) # Clusters * Samples
      
      for (i in 1:(length(clusters.names))) {
        geneID <- as.matrix(clusters.names[[i]])
        X <- data_in[geneID,]
        mu <- rowMeans(X)
        stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
        XNorm <- sweep(X,1,mu) # normalize X
        XNorm <- apply(XNorm, 2, function(x) x/stddev)
        SVD <- svd(XNorm, LINPACK = FALSE)
        eigengene.matrix[i,] <- t(SVD$v[,1])
      }
      eigengene.matrix = data.frame(eigengene.matrix)
      colnames(eigengene.matrix) = colnames(data_in)
      
      mergedCluster <- methods::new("QCMObject", clusters.id = clusters, clusters.names = clusters.names,
                                eigengene.matrix = eigengene.matrix)
      
      message("Done.")
      
#----------------------------------------------------------------------------
# lmQCM finished
#----------------------------------------------------------------------------
        
        
        
      # if("try-error" %in% class(t)) {
      #   removeModal()
      #   smartModal(error=T, title = "Error in lmQCM",
      #              content = sprintf("Too few genes to do lmQCM merging."))
      #   return()
      # }
      mergedCluster <- mergedCluster@clusters.id
      geneCharVector <- matrix(0, nrow = 0, ncol = length(mergedCluster))
      temp_eigengene <- matrix(0, nrow = length(mergedCluster), ncol = dim(finalExp())[2]) # Clusters * Samples

      temptext <- ""
      for (i in 1:(length(mergedCluster))) {
        vector <- as.matrix(mergedCluster[[i]])
        geneID <- vector
        print(i)
        print(vector)
        # ===== Calculate Eigengene Start
        X <- finalExp()[geneID,]
        mu <- rowMeans(X)
        stddev <- rowSds(as.matrix(X), na.rm=TRUE) # standard deviation with 1/(n-1)
        #normalize X:
        XNorm <- sweep(X,1,mu)
        XNorm <- apply(XNorm, 2, function(x) x/stddev)
        SVD <- svd(XNorm, LINPACK = FALSE)
        temp_eigengene[i,] <- t(SVD$v[,1])
        # ===== Calculate Eigengene Finished
        geneChar <- c(toString(i), finalSymChar()[vector])
        geneCharVector[i] <- list(geneChar)
        temptext <- paste(temptext, capture.output(cat(geneChar, sep=',')), sep="\n")
      }
      temptext <- substring(temptext, 2) # remove firstfinal_genes_str \n separater
      geneCharVector_global(geneCharVector)
      text.final(temptext)
      
      colnames(temp_eigengene) <- sampleID()
      eigengene_matrix(temp_eigengene)
      output$eigengene_matrix_select_row_ui <- renderUI({
        selectInput(inputId="eigengene_matrix_select_row", label="Please select row:",
                    choices = 1:dim(eigengene_matrix())[1], selected = 1, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      })

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
        return(eigengene_matrix())
      },rownames = TRUE, colnames = TRUE, na = "", bordered = TRUE, digits = 4)

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
    if(is.null(finalExp())){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
    if (length(finalSym()) > 0){
      #WGCNA
      smartModal(error=F, title = "Preview the power", content = "Running ...")
      finalExp.temp <- finalExp()
      row.names(finalExp.temp) <- finalSym()
      finalExp(finalExp.temp)
      datExpr <- t(finalExp()) # gene should be colnames, sample should be rownames
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
             xlab="Soft threshold (power)",ylab="Scale-free topology model fit, signed R²",type="n",
             main = paste("Scale independence"))
        text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red")
        # this line corresponds to using an R^2 cut-off of h
        abline(h=0.9,col="blue")
      })
      output$WGCNAPowerPlot2 <- renderPlot({
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1], sft$fitIndices[,5],
             xlab="Soft threshold (power)",ylab="Mean connectivity", type="n",
             main = paste("Mean connectivity"))
        text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
      })
      removeModal()
    }
  })

  observeEvent(input$action4_WGCNA,{
      if(is.null(finalExp())){
        smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
        return()
      }
      #WGCNA
      smartModal(error=F, title = "Using WGCNA to calculate merged clusters", content = "Calculating. This could take a while depend on number of genes. Please be patient.")
      finalExp.temp <- finalExp()
      row.names(finalExp.temp) <- finalSym()
      finalExp(finalExp.temp)
      datExpr <- t(finalExp()) # gene should be colnames, sample should be rownames
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
      t <- try(net <- blockwiseModules(datExpr, power = input$power,
                                      TOMType = "unsigned", minModuleSize = input$minModuleSize,     #30,
                                      reassignThreshold = input$reassignThreshold, mergeCutHeight =  input$mergeCutHeight,    # 0.25,
                                      numericLabels = TRUE, pamRespectsDendro = FALSE,
                                      saveTOMs = FALSE,
                                      saveTOMFileBase = "femaleMouseTOM"))
      if("try-error" %in% class(t)) {
        removeModal()
        smartModal(error=T, title = "Error in goodGenes",
                   content = sprintf("Too few genes with valid expression levels in the required number of samples."))
        return()
      }
      

      netcolors = net$colors
      matrixdata<- data.frame(cbind(finalSymChar(), netcolors))
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
                          main = sprintf("Power = %d, minModuleSize = %d, mergeCutHeight = %.4f",
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
      temp_eigengene <- matrix(0, nrow = length(unique(netcolors))-1, ncol = dim(finalExp())[2]) # Clusters * Samples
      for (i in 1: (length(unique(netcolors))-1) ){
        geneID <- which(matrixdata$netcolors == i)
        # ===== Calculate Eigengene Start
        X <- finalExp()[geneID,]
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
      geneCharVector_global(geneCharVector_withoutColor) # remove the color label
      text.final(temptext)
      colnames(temp_eigengene) <- sampleID()
      eigengene_matrix(temp_eigengene)
      output$eigengene_matrix_select_row_ui <- renderUI({
        selectInput(inputId="eigengene_matrix_select_row", label="Please select row:",
                    choices = 1:dim(eigengene_matrix())[1], selected = 1, multiple = FALSE,
                    selectize = TRUE, width = NULL, size = NULL)
      })
      
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
        return(eigengene_matrix())
      },rownames = TRUE, colnames = TRUE, na = "", bordered = TRUE, digits = 4)

      print('tab4')
      session$sendCustomMessage("download_cluster_ready","-")
      session$sendCustomMessage("myCallbackHandler", "tab4")
  })
  
  observeEvent(input$run_survival_analysis,{
    if(is.null(input$eigengene_matrix_select_row)){
      sendSweetAlert(session, title = "Error", text = "You haven't indicated which row of eigengene matrix will be dichotomized and analyzed.", type = "error",
                     btn_labels = "OK", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    event = input$survival_event
    time = input$survival_time
    event = as.numeric(unlist( regmatches(event, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", event)) ))
    time = as.numeric(unlist( regmatches(time, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", time)) ))
    
    if(length(event) != length(sampleID()) | length(time) != length(sampleID())){
      sendSweetAlert(session, title = "Error", text = sprintf("Number of OS/EFS events (%d) or number of OS/EFS times (%d) doesn't match number of samples.", length(event), length(time)),
                     type = "error", btn_labels = "OK", html = FALSE, closeOnClickOutside = TRUE)
      return()
    }
    
    row = as.numeric(input$eigengene_matrix_select_row)
    eigengene = eigengene_matrix()[row,]
    output$survival_analysis_results_ui <- renderUI({
      plotOutput("survival_plot", width = "100%", height = "400px", click = NULL,
                 dblclick = NULL, hover = NULL, hoverDelay = NULL,
                 hoverDelayType = NULL, brush = NULL, clickId = NULL, hoverId = NULL,
                 inline = FALSE)
    })
    mv = median(eigengene)
    group = cbind(length(time))
    for(j in 1:length(time)){
      if(eigengene[j] < mv){
        group[j] = 1
      }else{
        group[j] = 2
      }
    }
    
    mySurvTest = Surv(time, event)
    # logrank
    log1 = survdiff(mySurvTest ~ group)
    logrank.p = pchisq(log1$chisq, 1, lower.tail=FALSE)
    
    # KM Curve
    fit = survfit(mySurvTest ~ group)
    n1 = sum(group==1)
    leg1 = paste("Low risk(", n1, ")", sep = "")
    n2 = sum(group==2)
    leg2 = paste("High risk(", n2, ")", sep = "")
    
    output$survival_plot <- renderPlot({
      plot(fit, mark.time=TRUE, xlab = "Survival times", ylab = "Survival probability", lty = 1:2,
           col = 1:2, cex = 0.5)
      title(main = paste("Survival plot of eigengene module ", row, sep=""))
      grid()
      legend(x = "topright", legend = c(leg1, leg2), lty = 1:2,
             col = 1:2, cex = 0.65)
      text(10, 0.1, paste("p=", formatC(logrank.p, format="g", digits = 5), sep = ""),
           pos = 4, cex = 1)
    })
  })
  
  #   +------------------------------------------------------------+
  #   |
  #   |
  #   |                  GO Enrichment Analysis
  #   |
  #   |
  #   +--------------------------------
  observeEvent(input$action_finaldata4enrichr,{
    if(is.null(finalSym())){
      smartModal(error=T, title = "Operation Failed", content = "You have not selected any data. Please go to previous section.")
      return()
    }
    genes_str <- finalSym()
    genes_str <- unlist(strsplit(genes_str, " /// "))
    # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
    print("genes_str for enrich analysis: ")
    print(genes_str)
    final_genes_str(genes_str)
    updateTextAreaInput(session, "textareainput_GOEA",
                        label = paste(sprintf("Number of Genes: %d", length(final_genes_str())), input$controller),
                        value = paste(paste(final_genes_str(), collapse = '\n'), input$controller))
    
    enriched(enrichr(final_genes_str(), enrichr_dbs))
    
    Map(function(id) {
      dbres = enriched()[[enrichr_dbs[id]]]
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
      genes_str <- geneCharVector_global()[[cluster]]
      genes_str <- unlist(strsplit(genes_str, " /// "))
      # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
      print("genes_str for enrich analysis: ")
      print(genes_str[-1])
      final_genes_str(genes_str[-1])
      updateTextAreaInput(session, "textareainput_GOEA",
                          label = paste(sprintf("Number of Genes: %d", length(final_genes_str())), input$controller),
                          value = paste(paste(final_genes_str(), collapse = '\n'), input$controller))

      enriched(enrichr(final_genes_str(), enrichr_dbs))
      
      Map(function(id) {
        dbres = enriched()[[enrichr_dbs[id]]]
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
      genes_str <- geneCharVector_global()[[cluster]]
      genes_str <- unlist(strsplit(genes_str, " /// "))
      genes_str <- genes_str[-1]
      # genes_str <- c('PHF|14','RBM|3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
      print("genes_str for circos plot: ")
      print(genes_str)
      updateTextAreaInput(session, "textareainput_circos",
                          label = paste(sprintf("Number of Genes: %d", length(genes_str)), input$controller),
                          value = paste(paste(genes_str, collapse = '\n'), input$controller))
      
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
      output$circos_plot_ui_hg38 <- renderUI({
        plotOutput("circos_plot_component_hg38", width = input$circos_param_size, height = input$circos_param_size)
      })
      output$circos_plot_ui_hg19 <- renderUI({
        plotOutput("circos_plot_component_hg19", width = input$circos_param_size, height = input$circos_param_size)
      })
      output$circos_plot_component_hg38 <- renderPlot({
        factors_count = as.data.frame(hg38.ring.lengthsum)
        factors = factor(factors_count[,1], levels = factors_count[,1])
        xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
        rownames(xlim) = factors_count[,1]
        BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
        circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human genome (GRCh38/hg38) (genes)",
                         input$circos_param_size,
                         input$circos_param_genelink,
                         input$circos_param_genesymbol)
      })
      
      output$circos_plot_component_hg19 <- renderPlot({
        factors_count = as.data.frame(hg19.ring.lengthsum)
        factors = factor(factors_count[,1], levels = factors_count[,1])
        xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
        rownames(xlim) = factors_count[,1]
        BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
        circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                         input$circos_param_size,
                         input$circos_param_genelink,
                         input$circos_param_genesymbol)
      })
      removeModal()
      print('tab4_circos_plots')
      session$sendCustomMessage("myCallbackHandler", "tab4_circos_plots")
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
    output$circos_plot_ui_hg38 <- renderUI({
      plotOutput("circos_plot_component_hg38", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_ui_hg19 <- renderUI({
      plotOutput("circos_plot_component_hg19", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_component_hg38 <- renderPlot({
      factors_count = as.data.frame(hg38.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genome (GRCh38/hg38)",
                       input$circos_param_size,
                       input$circos_param_genelink,
                       input$circos_param_genesymbol)
    })
    
    output$circos_plot_component_hg19 <- renderPlot({
      factors_count = as.data.frame(hg19.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                       input$circos_param_size,
                       input$circos_param_genelink,
                       input$circos_param_genesymbol)
    })
    removeModal()
  })
  
  observeEvent(input$action_finaldata4circos,{
    if(is.null(data_final())){
      sendSweetAlert(session, title = "Insufficient Input Data", text = "Please finish previous steps.",
                     type = "error", btn_labels = "OK", html = F, closeOnClickOutside = T)
      return()
    }
    smartModal(error=F, title = "Processing", content = "We are working on your customized circos plot ...")
    # save(data_final(), file = "~/Desktop/datafinal.Rdata")
    genes_str <- levels(data_final()[,1])
    genes_str <- unlist(strsplit(genes_str, " /// "))
    # print("genes_str for circos plot: ")
    # print(genes_str)
    
    updateTextAreaInput(session, "textareainput_circos",
                        label = paste(sprintf("Number of Genes: %d", length(genes_str)), input$controller),
                        value = paste(paste(genes_str, collapse = '\n'), input$controller))
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
    output$circos_plot_ui_hg38 <- renderUI({
      plotOutput("circos_plot_component_hg38", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_ui_hg19 <- renderUI({
      plotOutput("circos_plot_component_hg19", width = input$circos_param_size, height = input$circos_param_size)
    })
    output$circos_plot_component_hg38 <- renderPlot({
      factors_count = as.data.frame(hg38.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg38.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg38", myTitle = "Human Genome (GRCh38/hg38)",
                       input$circos_param_size,
                       input$circos_param_genelink,
                       input$circos_param_genesymbol)
    })
    
    output$circos_plot_component_hg19 <- renderPlot({
      factors_count = as.data.frame(hg19.ring.lengthsum)
      factors = factor(factors_count[,1], levels = factors_count[,1])
      xlim = cbind(rep(0, dim(factors_count)[1]), factors_count[,2])
      rownames(xlim) = factors_count[,1]
      BED.data <- data.frame(hg19.matched[,c(4,6:7,10,14)])
      circlizeGenomics(BED.data, factors, xlim, mySpecies="hg19", myTitle = "Human Genome (GRCh37/hg19)",
                       input$circos_param_size,
                       input$circos_param_genelink,
                       input$circos_param_genesymbol)
    })
    removeModal()
    session$sendCustomMessage("myCallbackHandler", "tab4")
    session$sendCustomMessage("myCallbackHandler", "tab4_circos_plots")
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
        text_download = gsub(",", ",", noquote(text.final()))
        write.table(text_download, file, eol = "\r\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
      if (sep == 1){
        text_download = gsub(",", "\t", noquote(text.final()))
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
      write.table(eigengene_matrix(), file = file, append = FALSE, quote = TRUE, sep = separator,
                  eol = "\r\n", na = "NA", dec = ".", row.names = T,
                  col.names = NA, qmethod = c("escape", "double"),
                  fileEncoding = "")
    }
  )

  output$download_finaldata <- downloadHandler(
    filename = function() {
      name = "finaldata.csv"
    },
    content = function(file) {
      write.table(data_final(), file = file, append = FALSE, quote = TRUE, sep = ',',
                  eol = "\r\n", na = "NA", dec = ".", row.names = F,
                  col.names = T, qmethod = c("escape", "double"),
                  fileEncoding = "")
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename = 'GO_results.zip',
    content = function(fname_in) {
      separator <- switch(input$filetype3, "csv" = ',', "txt" = '\t')
      if(separator == ','){
        fs <- paste0(enrichr_dbs, '.csv')
      }
      else{
        fs <- paste0(enrichr_dbs, '.txt')
      }
      fs <- c(fs, 'genes_list.txt')
      for(i in 1:length(enrichr_dbs)){
        write.table(enriched()[[enrichr_dbs[i]]], file = fs[i], sep = separator, col.names = NA)
        print(fs[i])
      }
      if(length(final_genes_str()) > 0){
        write(final_genes_str(), file = 'genes_list.txt', sep = "\n")
      }
      else{
        write("You cannot believe I disabled this function. Haha.", file = 'genes_list.txt')
      }
      zip(zipfile=fname_in, files=fs)
      if(file.exists(paste0(fname_in, ".zip"))) {file.rename(paste0(fname_in, ".zip"), fname_in)}

    },
    contentType = "application/zip"
  )
}