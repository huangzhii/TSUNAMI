# Zhi Huang 03/16/2018
library(plotly)
navbarPage( theme = "style.css",

            # title="bioinfo tool",
            # title=div(img(src="images/iulogo.png",
            #               width = 20,
            #               style = "margin:0px 0px;"), "BioInfo Suite"),

            title=div(img(src="images/tbi_tsunami_logo.png",
                          height = 28,
                          style = "margin:0px 0px; padding-bottom: 5px")),
            tabPanel("Analysis",
                     navlistPanel(widths = c(2, 10),
                                  #theme = "style.css",
                                  #"GeneCoExpression!",
                                  
                                  tabPanel("INFO",
                                           sidebarLayout(
                                             position = "right",
                                             sidebarPanel(
                                               h4("Overview", style="color: STEELBLUE"),
                                               tabsetPanel(
                                                 id = 'overview_tabset',
                                                 tabPanel("NCBI GEO Database",
                                                          uiOutput("url_GEO"),
                                                          helpText("GEO is a public functional genomics data repository supporting MIAME-compliant data submissions. Numerous array- and sequence-based data are available for downstream analysis."),
                                                          
                                                          h5("Statistics", style="color: STEELBLUE"),
                                                          br(),
                                                          # Output: Histogram ----
                                                          plotlyOutput(outputId = "NCBI_GEO_Sample_Histogram", height = "200px"),
                                                          br(),
                                                          br(),
                                                          plotlyOutput(outputId = "NCBI_GEO_Release_Date", height = "200px")
                                                 )
                                               )
                                             ),
                                             mainPanel(
                                               h2("TBI-TSUNAMI:", style="color: STEELBLUE; font-size: 22px"),
                                               h2("Translational BIoinformatics Tool SUite for Network Analysis and MIning", style="color: STEELBLUE; font-size: 20px; margin: 0px"),
                                               
                                               h5("Introduction", style="color: STEELBLUE; padding-top: 10px"),
                                               "Recent technological advances in gene co-expression module discovery have allowed researcher to explorer highly correlated genes among groups. Thus TBI-TSUNAMI is developed to cater such needs. It has several highlight features and advantages:",
                                               
                                               tags$ul(
                                                 tags$li("Easy-access and real-time performance based on web server;"),
                                                 tags$li("Supporting multiple data formats and data preprocessing interface is bundled together;"),
                                                 tags$li("Multiple co-expression analysis interfaces available with a high flexibility of variable selection;"),
                                                 tags$li("Integrated downstream Enrichr enrichment analysis;"),
                                                 tags$li("All results are downloadable with multiple formats (CSV, txt, etc.) are supported.")
                                               ),
                                               "All of which bring convenience to researchers for multiple purposes.",
                                               h5("Pipeline Flowchart", style="color: STEELBLUE; padding-top: 10px"),
                                               tags$div(
                                                 tags$img(src='images/flowchart.png',
                                                          height="500",
                                                          alt="TBI-TSUNAMI Flowchart", class="center"),
                                                 style="text-align: center; padding: 20px"
                                               ),
                                               "Figure above: TBI-TSUNAMI in flowchart. Blue rectangles represent operations; Transparent rounded rectangles represent optional processes; Pink rounded rectangles represent Download processes; Black arrows represent regular procedure; Dashed orange arrows represent optional processes.",
                                               h5("Definition of Co-expression", style="color: STEELBLUE; padding-top: 10px"),
                                               "A gene co-expression network (GCN) is an undirected graph, where each node corresponds to a gene, and a pair of nodes is connected with an edge if there is a significant co-expression relationship between them.[1] Having gene expression profiles of a number of genes for several samples or experimental conditions, a gene co-expression network can be constructed by looking for pairs of genes which show a similar expression pattern across samples, since the transcript levels of two co-expressed genes rise and fall together across samples. Gene co-expression networks are of biological interest since co-expressed genes are controlled by the same transcriptional regulatory program, functionally related, or members of the same pathway or protein complex.",
                                               
                                               tags$head(
                                                 tags$script(HTML("(function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){(i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
                                                                  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
                                                                  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
                                                                  
                                                                  ga('create', 'UA-113406500-2', 'auto');
                                                                  ga('send', 'pageview');"))
                                                 ),
                                               tags$head(tags$script(HTML("document.title = 'TSUNAMI';"))), # rename the title by JS
                                               tags$head(tags$script('Shiny.addCustomMessageHandler("myCallbackHandler",
                                                                     function(typeMessage) {console.log(typeMessage)
                                                                     if(typeMessage == "tab1"){
                                                                     console.log("got here");
                                                                     $("a:contains(1. DataSet)").click();
                                                                     }
                                                                     if(typeMessage == "tab2"){
                                                                     $("a:contains(2. Data Preprocessing)").click();
                                                                     }
                                                                     if(typeMessage == "tab3"){
                                                                     $("a:contains(3. Choose Method)").click();
                                                                     }
                                                                     if(typeMessage == "tab4"){
                                                                     $("a:contains(4. Result)").click();
                                                                     }
                                                                     if(typeMessage == "tab5"){
                                                                     $("a:contains(5. GO Enrichment Analysis)").click();
                                                                     }
                                                                     });
                                                                     
                                                                     // disable download at startup.
                                                                     $(document).ready(function() {
                                                                     $("#downloadData1").attr("disabled", "true").attr("onclick", "return false;");
                                                                     $("#downloadData2").attr("disabled", "true").attr("onclick", "return false;");
                                                                     $("#downloadData3").attr("disabled", "true").attr("onclick", "return false;");
                                                                     $("#download_finaldata").attr("disabled", "true").attr("onclick", "return false;");
                                                                     
                                                                     Shiny.addCustomMessageHandler("download_cluster_ready", function(message) {
                                                                       $("#downloadData1").removeAttr("disabled").removeAttr("onclick");
                                                                       $("#downloadData2").removeAttr("disabled").removeAttr("onclick");
                                                                     });
                                                                     Shiny.addCustomMessageHandler("download_go_ready", function(message) {
                                                                     $("#downloadData3").removeAttr("disabled").removeAttr("onclick");
                                                                     });
                                                                     Shiny.addCustomMessageHandler("download_finaldata_ready", function(message) {
                                                                     $("#download_finaldata").removeAttr("disabled").removeAttr("onclick");
                                                                     });
                                                                     })
                                                                     ')),
                                               # Horizontal line ----
                                               tags$hr(),
                                               actionButton("action1", "Proceed")
                                               ) # end of mainPanel
                                                 ) # end of sidebarLayout
                                                 ), # end of tabPanel
                                  tabPanel("1. DataSet",
                                           # App title ----
                                           # titlePanel("Choose GSE Data or Upload Your Own Data"),
                                           singleton(
                                             tags$head(tags$script(src = "message-handler.js"))
                                           ),
                                           # Sidebar layout with input and output definitions ----
                                           sidebarLayout(
                                             position = "right",
                                             
                                             # Sidebar panel for inputs ----
                                             sidebarPanel(
                                               # Input: Select a file ----
                                               h4("File Uploader", style="color: STEELBLUE"),
                                               fileInput("csvfile", "Choose File",
                                                         multiple = FALSE,
                                                         accept = c("text/csv",
                                                                    "text/comma-separated-values,text/plain",
                                                                    ".csv", ".xlsx", ".xls")),
                                               
                                               # Include clarifying text ----
                                               helpText("Note: Maximum file size allowed for uploading is 300MB. If uploaded data is with .xlsx or .xls, separater can be any value, but please make sure data are located in Sheet1."),
                                               
                                               # Input: Checkbox if file has header ----
                                               checkboxInput("header", "Header", TRUE),
                                               
                                               fluidRow(
                                                 # Input: Select separator ----
                                                 column(6, radioButtons("sep", "Separator",
                                                                        choices = c(Comma = ",",
                                                                                    Semicolon = ";",
                                                                                    Tab = "\t",
                                                                                    Space = " "),
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
                                               h5("NCBI GEO Database"),
                                               tabsetPanel(
                                                 id = 'dataset',
                                                 tabPanel("Series", DT::dataTableOutput("mytable1"))
                                               ),
                                               
                                               tags$script("$(document).on('click', '#mytable1 button', function () {
                                                           Shiny.onInputChange('dataset_lastClickId',this.id)
                                                           });")
                                  )
                                             )
                                  
                                             ),
                                  tabPanel("2. Data Preprocessing",
                                           # titlePanel("Verifying & Cleaning Data"),
                                           # Sidebar layout with input and output definitions ----
                                           sidebarLayout(
                                             position = "right",
                                             # Sidebar panel for inputs ----
                                             sidebarPanel(
                                               tabsetPanel(
                                                 tabPanel("Main",
                                                   # Input: Select a file ----
                                                   h5("Choose Preview dimensions"),
                                                   helpText("Preview starting from the beginning to specific rows and columns.", style="margin: 0px"),
                                                   helpText("Default value when leave it blank: # of rows = 100, # of columns = 10.", style="color: STEELBLUE; font-size: 12px"),
                                                   
                                                   fluidRow(
                                                     column(6, numericInput("quicklook_row", "# of rows:", 100, step = 1, min = 1)),
                                                     column(6, numericInput("quicklook_col", "# of columns:", 10, step = 1, min = 1))
                                                   ),
                                                   # Horizontal line ----
                                                   # tags$hr(),
                                                   h5("Verify starting column and row of expression data"),
                                                   helpText("Choose starting column and row for expression data.", style="margin: 0px"),
                                                   helpText("Default value when leave them blank: starting row = 1, starting column = 2.", style="color: STEELBLUE; font-size: 12px"),
                                                   
                                                   fluidRow(
                                                     column(6, numericInput("starting_row", "starting row:", 1, step = 1, min = 1)),
                                                     column(6, numericInput("starting_col", "starting column:", 2, step = 1, min = 1))
                                                   ),
                                                   # Horizontal line ----
                                                   # tags$hr(),
                                                   h5("Verify Gene Symbol"),
                                                   helpText("We suppose Gene Symbol is in column 1.", style="margin: 0px"),
                                                   helpText("Default value when leave it blank: 1.", style="color: STEELBLUE; font-size: 12px"),
                                                   numericInput("starting_gene_row", "starting row:", 1, step = 1, min = 1),
                                                   
                                                   # Horizontal line ----
                                                   tags$hr(),
                                                   h5("Remove Genes"),
                                                   helpText("Remove data with lowest percentile absolute expression value shared by all samples. Then remove data with lowest percentile variance across samples.", style="margin: 0px"),
                                                   helpText("Default value when leave them blank: 0.", style="color: STEELBLUE; font-size: 12px"),
                                                   fluidRow(
                                                     column(6, numericInput("absolute_expval", "Lowest Absolute Percentile (%) To Remove:", 20, step = 1, min = 0)),
                                                     column(6, numericInput("variance_expval", "Lowest Variance Percentile (%) To Remove:", 10, step = 1, min = 0))
                                                   ),
                                                   checkboxInput("checkbox_NA", "Convert NA value to 0 in Expression Data", TRUE),
                                                   checkboxInput("checkbox_logarithm", "Take the log (e) of Expression Data (Default: Unchecked)", FALSE),
                                                   checkboxInput("checkbox_empty", "Remove rows with empty Gene Symbol", TRUE),
                                                   checkboxInput("checkbox_duplicated", "Keep only one row with largest mean expression value when Gene Symbol is duplicated", TRUE),
                                                   numericInput("max_gene_retain", "Maximum Number of Genes to Retain (i.e. Top N genes sorted by mean expression values among all samples. Leave blank for keeping all data):", 10000, step = 1000, min = 0)
                                                   
                                                 ), # EOF tabpanel main
                                                 tabPanel("Advanced",
                                                    h5("Choose Advanced Processes"),
                                                    verbatimTextOutput("summary_advanced"),
                                                    h5("Sorting"),
                                                    checkboxInput("sorting_adv_checkbox", "Sort expression data ascending to learn OS / EFS", F),
                                                    selectizeInput(
                                                      'choose_OS_EFS', 'Specify OS or EFS:',
                                                      choices = c("OS", "EFS")),
                                                    helpText("If yes, please select objective row: row index and range of columns."),
                                                    numericInput("row_osefs_ind", "Row of OS_IND/EFS_IND:", 9, step = 1, width = NULL, min = 1),
                                                    numericInput("row_osefs", "Row of OS/EFS:", 10, step = 1, width = NULL, min = 1),
                                                    numericInput("sort_col_start", "Starting Col:", 2, step = 1, width = NULL, min = 1),
                                                    checkboxInput("select_pval_adv_checkbox", "Pick Expression Data only with satisfied P-value.", F),
                                                    helpText("Calculated by median and the (non-central) Chi-Squared Distribution."),
                                                    numericInput("advance_selection_pvalue", "P-value smaller than:", 0.05, step = 0.001, width = NULL, min = 0)
                                                 )
                                               ), # EOF tabsetPanel
                                               actionButton("action3", "Continue to Co-Expression Analysis",style="color: WHITE; background-color: DODGERBLUE")
                                             ), # EOF siderbarPanel
                                             
                                             # Main panel for displaying outputs ----
                                             mainPanel(
                                               h5("Data Summary"),
                                               verbatimTextOutput("summary"),
                                               helpText("Convert Probe ID to Gene Symbol with Platform GPL*** (Optional for self-uploaded data):"),
                                               tags$span(style="color:STEELBLUE", "Be sure to verify (modify) Gene Symbol starting row on the Sidebar Panel!"),
                                               fluidRow(
                                                 column(6, textInput("platform_text", NULL, value = "Unknown", width = NULL, placeholder = NULL)),
                                                 column(6, actionButton("action_platform", "Convert"))
                                               ),
                                               h5("Data Preview"),
                                               tabsetPanel(
                                                 id = 'tabset2',
                                                 tabPanel("Original Data", DT::dataTableOutput("mytable4")),
                                                 tabPanel("Expression Value", DT::dataTableOutput("mytable5")),
                                                 tabPanel("Gene Symbol", tableOutput("mytable6"))
                                               )
                                             )
                                           ) # EOF siderbarLayout
                                       ), # EOF tabPanel 2. Data Preprocessing
                                  tabPanel("3. Choose Method",
                                           # titlePanel("Select Method for Gene Co-Expression Analysis"),
                                           
                                           tabsetPanel(
                                             id = "MethodList",
                                             tabPanel("lmQCM",
                                                      h5("lmQCM: An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph with Applications in Gene Co-Expression Module Discovery in Cancers"),
                                                      helpText("Zhang, Jie, and Kun Huang. \"Normalized ImQCM: An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph with Applications in Gene Co-Expression Module Discovery in Cancers.\" Cancer informatics 13 (2014): CIN-S14021."),
                                                      tags$hr(),
                                                      h5("Parameter Choosing"),
                                                      helpText("Gamma (γ) (Default = 0.55) controls the threshold for the initiation of each new module, lambda (λ) (Default = 1) and t (Default = 1) define the adaptive
                                                               threshold of the module density to ensure proper stopping
                                                               criterion for the greedy search for each module (Usually λ and t won't change), and beta (β) (Default = 0.4) is the
                                                               threshold for overlapping ratio for merging"),
                                                      fluidRow(
                                                        column(6, numericInput("gamma", "gamma (γ):", 0.55, step = 0.05)),
                                                        column(6, numericInput("lambda", "lambda (λ)", 1, step = 0.05))
                                                      ),
                                                      fluidRow(
                                                        column(6, numericInput("t", "t:", 1, step = 0.05)),
                                                        column(6, numericInput("beta", "beta (β):", 0.4, step = 0.05))
                                                      ),
                                                      fluidRow(
                                                        column(6, numericInput("minClusterSize", "Minimum Cluster Size:", 10, step = 1, width = NULL)),
                                                        column(6, selectizeInput(
                                                          'massiveCC', 'Calculation of Correlation Coefficient',
                                                          choices = c("pearson", "spearman")))
                                                      ),
                                                      
                                                      
                                                      # Horizontal line ----
                                                      tags$hr(),
                                                      actionButton("action4_lmQCM", "Confirm and Run",
                                                                   style="color: WHITE; background-color: DODGERBLUE"),
                                                      br(),
                                                      br(),
                                                      br(),
                                                      br()
                                                               ),
                                             tabPanel("WGCNA",
                                                      h5("WGCNA: An R package for weighted correlation network analysis"),
                                                      helpText("The WGCNA as an analysis method is described in: Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834"),
                                                      helpText("The package implementation is described in the article: Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559"),
                                                      tags$hr(),
                                                      
                                                      h5("Step 1: Pick Soft Thresholding"),
                                                      helpText("The soft thresholding, is a value used to power the correlation of the genes to that threshold. The assumption on that by raising the correlation to a power will reduce the noise of the correlations in the adjacency matrix. To pick up one threshold use the pickSoftThreshold function, which calculates for each power if the network resembles to a scale-free graph. The power which produce a higher similarity with a scale-free network is the one you should use."),
                                                      
                                                      actionButton("checkPower", "Check Power (β)"),
                                                      uiOutput("WGCNAPowerPlot1and2"),
                                                      
                                                      # Horizontal line ----
                                                      tags$hr(),
                                                      
                                                      h5("Step 2: Choose Parameters"),
                                                      helpText("Choose the power and remaining parameters. Default are as showned."),
                                                      helpText("power (β, Default = 6): The soft thresholding. 6 is large enough so that the resulting network exhibited approximate scale free topology."),
                                                      helpText("reassignThreshold (Default = 0): P-value ratio threshold for reassigning genes between modules."),
                                                      helpText("mergeCutHeight (Default = 0.25): Dendrogram cut height for module merging."),
                                                      helpText("verbose (Default = 3): Integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose."),
                                                      helpText("minModuleSize (Default = 10): Minimum module size for module detection."),
                                                      fluidRow(
                                                        column(6, numericInput("power", "power (β):", 6, step = 1, min = 1)),
                                                        column(6, numericInput("reassignThreshold", "Reassign Threshold", 0, step = 0.01))
                                                      ),
                                                      fluidRow(
                                                        column(6, numericInput("mergeCutHeight", "Merge Cut Height:", 0.25, step = 0.01)),
                                                        column(6, numericInput("verbose", "verbose:", 3, step = 1))
                                                      ),
                                                      numericInput("minModuleSize", "Minimum Module Size:", 10, step = 1, width = NULL),
                                                      
                                                      # Horizontal line ----
                                                      tags$hr(),
                                                      actionButton("action4_WGCNA", "Confirm and Run",
                                                                   style="color: WHITE; background-color: DODGERBLUE"),
                                                      uiOutput("WGCNAresultUI"),
                                                      br(),
                                                      br(),
                                                      br()
                                             ), # EOF WGCNA TAB
                                             tabPanel("Verify Final Data",
                                                      h5("Final Incoming Data"),
                                                      
                                                      fluidRow(
                                                        column(6,
                                                               helpText("You can verify the final incoming data and also download it."),
                                                               downloadButton('download_finaldata', 'Download Final Data (CSV)')),
                                                        column(6,
                                                               helpText("GO Enrichment Analysis for following all Genes."),
                                                               actionButton("action_finaldata4enrichr", "GO Enrichment Analysis",
                                                                               style="color: WHITE; background-color: DODGERBLUE"),
                                                               helpText("Warning: Directly process large # of genes may cause very slow GO process. We suggest user perform Co-expression clustering and do GO analysis with small amount of genes.", style="color: STEELBLUE; font-size: 12px"))
                                                      ),
                                                      
                                                      
                                                      h5("Data Preview"),
                                                      DT::dataTableOutput("mytable_finaldata")
                                                      
                                             ) # EOF TAB Verify Final Data
                                             )
                                           ),
                                  tabPanel("4. Result",
                                           mainPanel(
                                             h4("Download Results", style="color: STEELBLUE; padding-top: 10px"),
                                             fluidRow(
                                               column(6,
                                                      radioButtons("filetype1", "Merged Clusters with Gene Symbol:",
                                                                   choices = c("csv", "txt")),
                                                      downloadButton('downloadData1', 'Download')
                                               ),
                                               column(6,
                                                      radioButtons("filetype2", "Eigengene Matrix:",
                                                                   choices = c("csv", "txt")),
                                                      downloadButton('downloadData2', 'Download')
                                               )
                                             ),
                                             
                                             h4("Preview", style="color: STEELBLUE; padding-top: 10px"),
                                             tabsetPanel(
                                               id = 'tabset',
                                               tabPanel("Merged Clusters", DT::dataTableOutput("clusterResult")),
                                               tabPanel("Eigengene Matrix", tableOutput("mytable7"))
                                             ),
                                             
                                             tags$script("$(document).on('click', '#clusterResult button', function () {
                                                         Shiny.onInputChange('go_lastClickId',this.id)
                                                         });")
                                           )
                                           ),
                                  tabPanel("5. GO Enrichment Analysis",
                                           mainPanel(
                                             
                                             tabsetPanel(
                                               id = 'tabset_GOEA',
                                               tabPanel("Enrichr",
                                                        h4("Enrichment Analysis - by Enrichr", style="color: STEELBLUE"),
                                                        h5("Adjusted P-value (q-value):"),
                                                        helpText("The q-value is an adjusted p-value using the Benjamini-Hochberg method for correction for multiple hypotheses testing. Users can read more about this method, and why it is needed here:"),
                                                        helpText("Yoav Benjamini and Yosef Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological)
                                                                 Vol. 57, No. 1 (1995), pp. 289-300"),
                                                        h5("Relationship between P-value, Z-score, and combined score:"),
                                                        helpText("The combined score is a combination of the p-value and z-score calculated by multiplying the two scores as follows:
                                                                 c = ln(p) * z
                                                                 Where c is the combined score, p is the p-value computed using Fisher's exact test, and z is the z-score computed to assess the deviation from the expected rank. The combined score provides a compromise between both methods and in several benchmarks we show that it reports the best rankings when compared with the other scoring schemes."),
                                                        
                                                        
                                                        tabsetPanel(
                                                          id = 'tabset',
                                                          tabPanel("GO_Biological_Process_2017b",
                                                                   DT::dataTableOutput("mytable_Enrichr_1")),
                                                          tabPanel("GO_Molecular_Function_2017b",
                                                                   DT::dataTableOutput("mytable_Enrichr_2")),
                                                          tabPanel("GO_Cellular_Component_2017b",
                                                                   DT::dataTableOutput("mytable_Enrichr_3")),
                                                          tabPanel("Jensen_DISEASES",
                                                                   DT::dataTableOutput("mytable_Enrichr_4")),
                                                          tabPanel("Reactome_2016",
                                                                   DT::dataTableOutput("mytable_Enrichr_5")),
                                                          tabPanel("KEGG_2016",
                                                                   DT::dataTableOutput("mytable_Enrichr_6")),
                                                          tabPanel("Transcription_Factor_PPIs",
                                                                   DT::dataTableOutput("mytable_Enrichr_7")),
                                                          tabPanel("Genome_Browser_PWMs",
                                                                   DT::dataTableOutput("mytable_Enrichr_8")),
                                                          tabPanel("TRANSFAC_and_JASPAR_PWMs",
                                                                   DT::dataTableOutput("mytable_Enrichr_9")),
                                                          tabPanel("ENCODE_TF_ChIP-seq_2015",
                                                                   DT::dataTableOutput("mytable_Enrichr_10")),
                                                          tabPanel("Chromosome_Location (Cytoband)",
                                                                   DT::dataTableOutput("mytable_Enrichr_11")),
                                                          tabPanel("miRTarBase_2017",
                                                                   DT::dataTableOutput("mytable_Enrichr_12")),
                                                          tabPanel("TargetScan_microRNA_2017",
                                                                   DT::dataTableOutput("mytable_Enrichr_13")),
                                                          tabPanel("ChEA_2016",
                                                                   DT::dataTableOutput("mytable_Enrichr_14"))
                                                        )
                                                        )#,
                                               #
                                               # tabPanel("Target Gene Symbols",
                                               #          h5("The target gene symbols allow users to copy and use in other GO analysis website."),
                                               #          textAreaInput("textareainput_GOEA", "Gene Symbols", value = "", width = '300px', height = '400px', placeholder = NULL)
                                               # )
                                                        )
                                             
                                           ),
                                           sidebarPanel(
                                             h4("Download Results & Further Analysis", style="color: STEELBLUE"),
                                             radioButtons("filetype3", "Choose file type:",
                                                          choices = c("csv", "txt")),
                                             downloadButton('downloadData3', 'Download'),
                                             br(),
                                             br(),
                                             textAreaInput("textareainput_GOEA", "Gene Symbols", value = "", width = 'auto', height = '300px', placeholder = NULL),
                                             helpText("The target gene symbols allow users to copy and use in other GO analysis website."),
                                             uiOutput("url_toppgene"),
                                             uiOutput("url_david"),
                                             uiOutput("url_enrichr"),
                                             uiOutput("url_gorilla")
                                           )
                                  )
                                  )  #, style='width: 80%'
                     
            ),
            tabPanel("Tutorial",
                     h3("Tutorial", style="color: STEELBLUE; padding-bottom: 20px"),
                     h4("Video Tutorial", style="color: STEELBLUE; padding-bottom: 20px"),
                     tags$div(
                       HTML('<iframe width="720" height="480" src="https://www.youtube.com/embed/rRIRMW_RRS4" frameborder="0" allowfullscreen></iframe>'),
                       style="text-align: center; padding: 20px"
                     ),
                     h4("Google Slides", style="color: STEELBLUE; padding-bottom: 20px"),
                     h4("Github", style="color: STEELBLUE; padding-bottom: 20px"),
                     h4("Report Bugs", style="color: STEELBLUE; padding-bottom: 20px")
            ),
            tabPanel("FAQ",
                     h3("Frequently Asked Questions", style="color: STEELBLUE; padding-bottom: 20px"),
                     h4("General Questions", style="color: STEELBLUE; padding-bottom: 20px"),
                     h5("What is the TBI-TSUNAMI website?", style="color: STEELBLUE"),
                     p("The TBI-TSUNAMI (Translational Bioinformatics Tool Suite for Network Analysis and Mining) was developed at Indiana University School of Medicine."),
                     h5("How do I get started?", style="color: STEELBLUE"),
                     p("Check our tutorial.")
            ),
            tabPanel("News",
                      h3("News", style="color: STEELBLUE; padding-bottom: 20px"),
                       h4("March 20, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
                       tags$ul(
                         tags$li("R package 'lmQCM' was released to CRAN."),
                         tags$li("Create flowchart.")
                       ),
                       h4("March 16, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
                       tags$ul(
                         tags$li("Renamed our website as TBI-TSUNAMI."),
                         tags$li("Various of bugs are fixed.")
                       ),
                      h4("March 02, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
                      tags$ul(
                        tags$li("First prototype platform has been founded.")
                      )
            ),
            tabPanel("About",
                      h3("About Us", style="color: STEELBLUE; padding-bottom: 20px"),
                      "The TBI-TSUNAMI (Translational Bioinformatics Tool Suite for Network Analysis and Mining) was developed at Indiana University School of Medicine.",
                     
                     tags$div(
                       tags$img(src='images/IUSM2.png',
                                height="100",
                                alt="TBI-TSUNAMI", class="center", style="padding: 30px"),
                       tags$img(src='images/regenstrief.png',
                                height="100",
                                alt="TBI-TSUNAMI", class="center", style="padding: 30px"),
                       style="text-align: center; padding: 20px"
                     ),
                     
                      h4("Development Team", style="color: STEELBLUE; padding-bottom: 20px"),
                      h5("Prof. Kun Huang's Laboratory", style="color: STEELBLUE"),
                      tags$ul(
                        tags$li("Zhi Huang"),
                        tags$li("Zhi Han"),
                        tags$li("Jie Zhang"),
                        tags$li("Kun Huang")
                      ),
                      h4("Publications", style="color: STEELBLUE; padding-bottom: 20px"),
                      "Nice paper",
                      h4("Funding for the TBI-TSUNAMI is or has been provided by:", style="color: STEELBLUE; padding-bottom: 20px"),
                      tags$ul(
                        tags$li("Good funding"),
                        tags$li("Nice funding")
                      )
                   
            )
            
            
            # navbarMenu(
            #   "More",
            #   tabPanel("Developer",
            #            h4("Author Information"),
            #            helpText("Indiana University School of Medicine"),
            #            h4("Publication"),
            #            helpText("Please cite ...")
            #            )
            # )
  )
