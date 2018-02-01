# Zhi Huang 01/26/2017

navbarPage( theme = "style.css",
  "GeneCoExpression!",
  
  tabPanel("Main Page",
           titlePanel("Gene Co-Expression Analysis Site"),
           tags$head(tags$script('Shiny.addCustomMessageHandler("myCallbackHandler",
                                 function(typeMessage) {console.log(typeMessage)
                                 if(typeMessage == "tab1"){
                                 console.log("got here");
                                 $("a:contains(1. DataSet)").click();
                                 }
                                 if(typeMessage == "tab2"){
                                 $("a:contains(2. Data Cleaning)").click();
                                 }
                                 if(typeMessage == "tab3"){
                                 $("a:contains(3. Choose Method)").click();
                                 }
                                 if(typeMessage == "tab4"){
                                 $("a:contains(4. Summary and Download)").click();
                                 }
                                 });
                                 ')),
           # Horizontal line ----
           tags$hr(),
           actionButton("action1", "Proceed")
           ),
  tabPanel("1. DataSet",
           # App title ----
           titlePanel("Choose GSE Data or Upload Your Own Data"),
           singleton(
             tags$head(tags$script(src = "message-handler.js"))
           ),
           # Sidebar layout with input and output definitions ----
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               # Input: Select a file ----
               h4("File uploader"),
               fileInput("csvfile", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               
               # Include clarifying text ----
               helpText("Note: Maximum csv file size allowed for uploading is 300MB."),
               
               # Input: Checkbox if file has header ----
               checkboxInput("header", "Header", TRUE),
               
               fluidRow(
                 # Input: Select separator ----
                 column(6, radioButtons("sep", "Separator",
                                        choices = c(Comma = ",",
                                                    Semicolon = ";",
                                                    Tab = "\t"),
                                        selected = ",")),
                 # Input: Select quotes ----
                 column(6, radioButtons("quote", "Quote",
                                        choices = c(None = "",
                                                    "Double Quote" = '"',
                                                    "Single Quote" = "'"),
                                        selected = '"'))
               ),
               
               # Horizontal line ----
               tags$hr(),
               actionButton("action2", "Confirm when Complete")
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               h4("NCBI GEO Data"),
               tabsetPanel(
                 id = 'dataset',
                 tabPanel("Series", DT::dataTableOutput("mytable1"))
                 # tabPanel("Datasets", DT::dataTableOutput("mytable2")),
                 # tabPanel("Samples", DT::dataTableOutput("mytable3"))
               ),
               
               tags$script("$(document).on('click', '#mytable1 button', function () {
                        Shiny.onInputChange('lastClickId',this.id)
                       });")
             )
           )
           
         ),
  tabPanel("2. Data Cleaning",
           titlePanel("Verifying & Cleaning Data"),
           
           # Sidebar layout with input and output definitions ----
           sidebarLayout(
             
             # Sidebar panel for inputs ----
             sidebarPanel(
               # Input: Select a file ----
               h4("Choose Quicklook dimensions"),
               helpText("Quicklook starting from the beginning to specific rows and columns."),
               
               fluidRow(
                 column(6, numericInput("quicklook_row", "Number of rows:", 100, step = 1, min = 1)),
                 column(6, numericInput("quicklook_col", "Number of column:", 10, step = 1, min = 1))
               ),
               # Horizontal line ----
               tags$hr(),
               h4("Verify starting column and row of expression data"),
               helpText("This is assigner. Choose starting column and row for expression data. Index starting from 1."),
               
               fluidRow(
                 column(6, numericInput("starting_row", "starting row:", 1, step = 1, min = 1)),
                 column(6, numericInput("starting_col", "starting column:", 2, step = 1, min = 1))
               ),
               # Horizontal line ----
               tags$hr(),
               h4("Verify Gene ID"),
               helpText("We suppose Gene ID is in column 1."),
               numericInput("starting_gene_row", "starting row:", 1, step = 1, min = 1),
               # Horizontal line ----
               tags$hr(),
               h4("Remove Genes"),
               helpText("Remove data with lowest percentile absolute expression value shared by all samples. Then remove data with lowest percentile variance across samples."),
               fluidRow(
                 column(6, numericInput("absolute_expval", "Lowest Absolute Percentile (%) To Remove:", 20, step = 1, min = 0)),
                 column(6, numericInput("variance_expval", "Lowest Variance Percentile (%) To Remove:", 10, step = 1, min = 0))
               ),
               numericInput("max_gene_retain", "Maximum Number of Genes to Retain:", 20000, step = 5000, min = 0),
               checkboxInput("NAconverter", "Convert NA value to 0", TRUE),
               actionButton("action3", "Continue")
             ),
             
             # Main panel for displaying outputs ----
             mainPanel(
               h4("Data Quicklook"),
               tabsetPanel(
                 id = 'tabset2',
                 tabPanel("Original Data", DT::dataTableOutput("mytable4")),
                 tabPanel("Expression Value", DT::dataTableOutput("mytable5")),
                 tabPanel("Gene ID", tableOutput("mytable6"))
               )
             )
           )
  ),
  tabPanel("3. Choose Method",
           titlePanel("Select Method for Gene Co-Expression Analysis"),
           navlistPanel(
             "Method List",
             tabPanel("lmQCM",
                      h3("lmQCM"),
                      h4("An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph with Applications in Gene Co-Expression Module Discovery in Cancers"),
                      tags$hr(),
                      fluidRow(
                        column(6, numericInput("gamma", "gamma:", 0.55, step = 0.05)),
                        column(6, numericInput("lambda", "lambda", 1, step = 0.05))
                      ),
                      fluidRow(
                        column(6, numericInput("t", "t:", 1, step = 0.05)),
                        column(6, numericInput("beta", "beta:", 0.4, step = 0.05))
                      ),
                      numericInput("minClusterSize", "Minimum Cluster Size:", 10, step = 1, width = NULL),
                      
                      # Horizontal line ----
                      tags$hr(),
                      actionButton("action4_lmQCM", "Confirm and Run!")
             ),
             tabPanel("WGCNA",
                      h3("WGCNA"),
                      h4("An R package for weighted correlation network analysis"),
                      tags$hr(),
                      
                      h4("Step 1: Power Visualization"),
                      helpText("Visualize the data"),
                      fluidRow(
                        column(6, plotOutput('WGCNAPowerPlot1', height = "300px")),
                        column(6, plotOutput('WGCNAPowerPlot2', height = "300px"))
                      ),
                      actionButton("checkPower", "Check Power"),
                      
                      # Horizontal line ----
                      tags$hr(),
                      
                      h4("Step 2: Choose Parameters"),
                      helpText("Choose the power and remaining parameters"),
                      fluidRow(
                        column(6, numericInput("power", "power:", 6, step = 1, min = 1)),
                        column(6, numericInput("reassignThreshold", "Reassign Threshold", 0, step = 0.01))
                      ),
                      fluidRow(
                        column(6, numericInput("mergeCutHeight", "Merge Cut Height:", 0.25, step = 0.01)),
                        column(6, numericInput("verbose", "verbose:", 3, step = 1))
                      ),
                      numericInput("minModuleSize", "Minimum Module Size:", 10, step = 1, width = NULL),
                      
                      # Horizontal line ----
                      tags$hr(),
                      actionButton("action4_WGCNA", "Confirm and Run!"),
                      plotOutput('WGCNAresult', height = "500px")
             ),
             tabPanel("More",
                      h3("More...")
             )
           )
         ),
  tabPanel("4. Summary and Download",
              mainPanel(
                h4("Download Results"),
                radioButtons("filetype", "Choose file type and download processed data:",
                             choices = c("csv", "txt")),
                downloadButton('downloadData', 'Download'),
                # Horizontal line ----
                tags$hr(),
                h4("Preview"),
                tableOutput("clusterResult")
              )
           ),
  navbarMenu(
    "More",
    tabPanel("Login",
             DT::dataTableOutput("table"))
  )
)
