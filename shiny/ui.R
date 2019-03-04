# Zhi Huang 03/28/2018
library(plotly)
library(shinyWidgets)
navbarPage( theme = "style.css",
            
            # title="bioinfo tool",
            # title=div(img(src="images/iulogo.png",
            #               width = 20,
            #               style = "margin:0px 0px;"), "BioInfo Suite"),
            
            title=div(a(img(src="images/tsunami_logo.png",
                            height = 28,
                            style = "margin:0px 0px; padding-bottom: 5px"), href="https://apps.medgen.iupui.edu/rsc/tsunami")),
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
                                                          h5(a("Gene Expression Omnibus", href="https://www.ncbi.nlm.nih.gov/geo/", target="_blank")),
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
                                               h2("TSUNAMI:", style="color: STEELBLUE; font-size: 22px"),
                                               h2("Translational Bioinformatics Tool SUite for Network Analysis and MIning", style="color: STEELBLUE; font-size: 20px; margin: 0px"),
                                               
                                               h5("Introduction", style="color: STEELBLUE; padding-top: 10px"),
                                               "Gene co-expression network (GCN) mining aims to mine gene modules with highly correlated expression profiles across sample cohorts. It may help to reveal latent molecular interactions, identify novel gene functions, pathways and drug targets, as well as providing disease mechanistic insights on for biological researchers. TSUNAMI is developed to allow biological researchers with no programing background to perform GCN mining themselves. It has several highlight features and advantages:",
                                               
                                               tags$ul(
                                                 tags$li("User friendly interface, easy-access and real-time co-expression network mining based on web server;"),
                                                 tags$li("Direct access and search of GEO database as well as user-input expression matrix for network mining;"),
                                                 tags$li("Support multiple data formats and data preprocessing interface is bundled together;"),
                                                 tags$li("Multiple co-expression analysis tools available with a high flexibility of variable selection;"),
                                                 tags$li("Integrated downstream Enrichr GO enrichment analysis and link to other GO tools as well;"),
                                                 tags$li("All results can be downloaded with multiple formats (CSV, txt, etc.).")
                                               ),
                                               "All of which bring convenience to researchers for multiple purposes.",
                                               h5("Pipeline Flowchart", style="color: STEELBLUE; padding-top: 10px"),
                                               tags$div(
                                                 tags$img(src='images/flowchart_concise_flat4.png',
                                                          width="600",
                                                          alt="TSUNAMI Flowchart", class="center"),
                                                 style="text-align: center; padding: 20px"
                                               ),
                                               "Figure above: TSUNAMI in flowchart. Blue blocks represent operations; Pink rounded rectangles represent Download processes; Dashed arrow means optional process.",
                                               
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
                                                                     if(typeMessage == "tab4_circos_plots"){
                                                                     $("a:contains(Circos Plots)").click();
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
                                               actionButton("action1", "Proceed"),
                                               br(),
                                               br(),
                                               br()
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
                                                                    ".csv", ".xlsx")),
                                               
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
                                               h5("Database"),
                                               helpText("Example GSE Microarray Data: GSE17537; GSE88882; GSE98761; GSE40294; GSE73119; GSE31399; GSE21361; GSE13002; GSE4309; GSE61084; GSE61085.", style="font-size: 12px"),
                                               # helpText("Example RNA-seq Expression Data: ", a("TCGA-BLCA",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-BLCA.csv",target="_blank"), a("TCGA-BRCA",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-BRCA.csv",target="_blank"),
                                               #          a("TCGA-CESC",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-CESC.csv",target="_blank"),a("TCGA-ESCA",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-ESCA.csv",target="_blank"),
                                               #          a("TCGA-HNSC",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-HNSC.csv",target="_blank"),a("TCGA-KIRC",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-KIRC.csv",target="_blank"),
                                               #          a("TCGA-KIRP",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-KIRP.csv",target="_blank"),a("TCGA-LIHC",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/TCGA-LIHC.csv",target="_blank"),
                                               #          style="font-size: 12px"),
                                               helpText("Example Single-cell RNA-seq Data: ", a("GSE59739_DataTable",href="http://web.ics.purdue.edu/~huang898/TSUNAMI_data/GSE59739_DataTable.txt",target="_blank"), style="font-size: 12px"),
                                               tabsetPanel(
                                                 id = 'dataset',
                                                 tabPanel("TCGA mRNA-seq Data",
                                                          h5("Table of illuminahiseq rnaseqv2 RSEM genes normalized mRNA-seq data",
                                                             style="color: black; font-size: 14px; font-weight: bold"),
                                                          DT::dataTableOutput("mytable0")),
                                                 tabPanel("GEO Series Matrix", DT::dataTableOutput("mytable1"))
                                               ),
                                               
                                               tags$script("$(document).on('click', '#mytable1 button', function () {
                                                           Shiny.onInputChange('dataset_lastClickId_mytable1',this.id)
                                                           });
                                                           $(document).on('click', '#mytable0 button', function () {
                                                           Shiny.onInputChange('dataset_lastClickId_mytable0',this.id)
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
                                                 tabPanel("Basic",
                                                          # Input: Select a file ----
                                                          # h5("Choose Preview dimensions"),
                                                          # helpText("Preview starting from the beginning to specific rows and columns.", style="margin: 0px"),
                                                          # helpText("Default value when leave it blank: # of rows = 100, # of columns = 10.", style="color: STEELBLUE; font-size: 12px"),
                                                          # 
                                                          # fluidRow(
                                                          #   column(6, numericInput("quicklook_row", "# of rows:", 100, step = 1, min = 1)),
                                                          #   column(6, numericInput("quicklook_col", "# of columns:", 10, step = 1, min = 1))
                                                          # ),
                                                          # Horizontal line ----
                                                          # tags$hr(),
                                                          h5("Verify starting column and row of expression data"),
                                                          helpText("Choose starting column and row for expression data.", style="margin: 0px"),
                                                          helpText("Default value when leave them blank: starting row = 1, starting column = 2.", style="color: STEELBLUE; font-size: 12px"),
                                                          
                                                          fluidRow(
                                                            column(6, numericInput("starting_row", "Gene and Expression starting row:", 1, step = 1, min = 1)),
                                                            column(6, numericInput("starting_col", "Expression starting column:", 2, step = 1, min = 1))
                                                          ),
                                                          h5("Convert Probe ID to Gene Symbol"),
                                                          helpText("Convert Probe ID to Gene Symbol with Platform GPL*** (Optional for self-uploaded data):"),
                                                          tags$span(style="color:STEELBLUE", "Be sure to verify (modify) Gene Symbol."),
                                                          fluidRow(
                                                            column(8, textInput("platform_text", NULL, value = "Unknown", width = NULL, placeholder = NULL)),
                                                            column(4, actionButton("action_platform", "Convert"))
                                                          ),
                                                          # h5("Verify Gene Symbol"),
                                                          # helpText("We suppose Gene Symbol is in column 1.", style="margin: 0px"),
                                                          # helpText("Default value when leave it blank: 1.", style="color: STEELBLUE; font-size: 12px"),
                                                          # numericInput("starting_gene_row", "starting row:", 1, step = 1, min = 1),
                                                          
                                                          # Horizontal line ----
                                                          tags$hr(),
                                                          h5("Remove Genes"),
                                                          helpText("Remove rows with lowest percentile mean expression value shared by all samples. Then remove data with lowest percentile variance across samples.", style="margin: 0px"),
                                                          helpText("Default value when leave them blank: 0.", style="color: STEELBLUE; font-size: 12px"),
                                                          fluidRow(
                                                            column(6, numericInput("mean_expval", "Lowest Mean Percentile (%) To Remove:", 50, step = 1, min = 0, max = 100)),
                                                            column(6, numericInput("variance_expval", "Lowest Variance Percentile (%) To Remove:", 50, step = 1, min = 0, max = 100))
                                                          ),
                                                          checkboxInput("checkbox_NA", "Convert NA value to 0 in Expression Data", TRUE),
                                                          checkboxInput("checkbox_logarithm", "Take the log2(x+1) of Expression Data x (Default: Unchecked)", FALSE),
                                                          checkboxInput("checkbox_empty", "Remove rows with empty Gene Symbol", TRUE),
                                                          checkboxInput("checkbox_duplicated", "Keep only one row with largest mean expression value when Gene Symbol is duplicated", TRUE),
                                                          # numericInput("max_gene_retain", "Maximum Number of Genes to Retain (i.e. Top N genes sorted by mean expression values among all samples. Leave blank for keeping all data):", 10000, step = 1000, min = 0)
                                                          actionButton("action3", "Continue to Co-Expression Analysis",style="color: WHITE; background-color: DODGERBLUE")
                                                          
                                                 ),
                                                 
                                                 tabPanel("Advanced",
                                                          h5("Select Samples Subgroup", style="color: black; font-size: 14px; font-weight: bold"),
                                                          actionLink("data_sample_subgroup_selectall","Select/Deselect All"),
                                                          uiOutput("data_sample_subgroup_ui")
                                                 )
                                                 # EOF tabpanel main
                                                 # tabPanel("Advanced",
                                                 #          h5("Choose Advanced Processes"),
                                                 #          checkboxInput("sorting_adv_checkbox", "Sort expression data ascending to learn OS / EFS", F),
                                                 #          conditionalPanel(condition = "input.sorting_adv_checkbox == 1",
                                                 #                           selectizeInput(
                                                 #                             'choose_OS_EFS', 'Specify OS or EFS:',
                                                 #                             choices = c("OS", "EFS")),
                                                 #                           helpText("If yes, please select objective row: row index and range of columns. Please refer OS/EFS position from Original Data (not Verified Data)."),
                                                 #                           helpText("OS_IND/EFS_IND must either valued 0 or 1. OS/EFS must be numeric."),
                                                 #                           numericInput("row_osefs_ind", "Row of OS_IND/EFS_IND:", 9, step = 1, width = NULL, min = 1),
                                                 #                           numericInput("row_osefs", "Row of OS/EFS:", 10, step = 1, width = NULL, min = 1),
                                                 #                           numericInput("sort_col_start", "Starting Col:", 2, step = 1, width = NULL, min = 1)
                                                 #          ),
                                                 #          checkboxInput("select_pval_adv_checkbox", "Pick Expression Data only with satisfied P-value.", F),
                                                 #          conditionalPanel(condition = "input.select_pval_adv_checkbox == 1",
                                                 #                           helpText("Calculated by median and the (non-central) Chi-Squared Distribution."),
                                                 #                           numericInput("advance_selection_pvalue", "P-value smaller than:", 0.05, step = 0.001, width = NULL, min = 0)
                                                 #          )
                                                 #          
                                                 # )
                                               ) # EOF tabsetPanel
                                             ), # EOF siderbarPanel
                                             
                                             # Main panel for displaying outputs ----
                                             mainPanel(
                                               h5("Data Summary"), 
                                               verbatimTextOutput("summary"),
                                               h5("Data Preview"),
                                               tabsetPanel(
                                                 id = 'tabset2',
                                                 tabPanel("Original Data", DT::dataTableOutput("mytable4")),
                                                 tabPanel("Verified Data", DT::dataTableOutput("mytable5"))
                                               )
                                             )
                                           ) # EOF siderbarLayout
                                  ), # EOF tabPanel 2. Data Preprocessing
                                  tabPanel("3. Choose Method",
                                           # titlePanel("Select Method for Gene Co-Expression Analysis"),
                                           
                                           tabsetPanel(
                                             id = "MethodList",
                                             tabPanel("lmQCM",
                                                      h4("lmQCM: An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph", style="color: STEELBLUE; padding-top: 10px"),
                                                      tags$div(
                                                        a(tags$img(src='images/lmQCM_logo.png',
                                                                   height="60",
                                                                   alt="lmQCM", class="center", style="padding: 5px"), href="https://CRAN.R-project.org/package=lmQCM", target="_blank"),
                                                        style="text-align: center; padding: 5px"
                                                      ),
                                                      helpText("If you benefit from the results, please cite:"),
                                                      a("Zhang, Jie, and Kun Huang. \"Normalized ImQCM: An Algorithm for Detecting Weak Quasi-Cliques in Weighted Graph with Applications in Gene Co-Expression Module Discovery in Cancers.\" Cancer informatics 13 (2014): CIN-S14021.",href="http://journals.sagepub.com/doi/abs/10.4137/CIN.S14021",target="_blank"),
                                                      tags$hr(),
                                                      h5("Parameter Choosing"),
                                                      helpText("Gamma (γ) (Default = 0.7, Recommend: 0.70 ~ 0.75) controls the threshold for the initiation of each new module, lambda (λ) (Default = 1) and t (Default = 1) define the adaptive
                                                               threshold of the module density to ensure proper stopping
                                                               criterion for the greedy search for each module (Usually λ and t won't change), and beta (β) (Default = 0.4) is the
                                                               threshold for overlapping ratio for merging"),
                                                      helpText("Weight Normalization is to normalize the correlation matrix (default: Not selected). However we recommend to check it while the expression data comes from microarray."),
                                                      prettyCheckbox(inputId = "lmQCM_weight_normalization", label = "Weight Normalization",
                                                                     value = F, status = "default", icon = icon("check")),
                                                      fluidRow(
                                                        column(6, numericInput("gamma", "gamma (γ):", 0.7, step = 0.05)),
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
                                                      h4("WGCNA: An R package for weighted correlation network analysis", style="color: STEELBLUE; padding-top: 10px"),
                                                      helpText("If you benefit from the results, please cite:"),
                                                      helpText("The WGCNA as an analysis method is described in:",a("Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834",href="https://econpapers.repec.org/article/bpjsagmbi/v_3a4_3ay_3a2005_3ai_3a1_3an_3a17.htm", target="_blank")),
                                                      helpText("The package implementation is described in the article:", a("Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559",href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559", target="_blank")),
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
                                                      # helpText("verbose (Default = 3): Integer level of verbosity. Zero means silent, higher values make the output progressively more and more verbose."),
                                                      helpText("minModuleSize (Default = 10): Minimum module size for module detection."),
                                                      fluidRow(
                                                        column(6, numericInput("power", "power (β):", 6, step = 1, min = 1)),
                                                        column(6, numericInput("reassignThreshold", "Reassign Threshold", 0, step = 0.01))
                                                      ),
                                                      fluidRow(
                                                        column(6, numericInput("mergeCutHeight", "Merge Cut Height:", 0.25, step = 0.01)),
                                                        # column(6, numericInput("verbose", "verbose:", 3, step = 1))
                                                        column(6, numericInput("minModuleSize", "Minimum Module Size:", 10, step = 1, width = NULL))
                                                      ),
                                                      
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
                                                      h4("Final Incoming Data", style="color: STEELBLUE; padding-top: 10px"),
                                                      
                                                      fluidRow(
                                                        column(4,
                                                               helpText("You can verify the final incoming data and also download it."),
                                                               downloadButton('download_finaldata', 'Download Final Data (CSV)')),
                                                        column(4,
                                                               helpText("GO Enrichment Analysis for following all Genes."),
                                                               actionButton("action_finaldata4enrichr", "GO Enrichment Analysis",
                                                                            style="color: WHITE; background-color: DODGERBLUE"),
                                                               helpText("Warning: Directly process large # of genes may cause very slow GO process. We suggest user perform Co-expression clustering and do GO analysis with small amount of genes.", style="color: STEELBLUE; font-size: 12px")),
                                                        column(4,
                                                               helpText("Circos Plot"),
                                                               actionButton("action_finaldata4circos", "Circos Plot for All Genes",
                                                                            style="color: WHITE; background-color: #FFC300"),
                                                               helpText("When finished, go to 4. Result Circos plots section.", style="color: #D5A200; font-size: 12px"),
                                                               helpText("We strongly recomment user clean the genes first through our Data Preprocessing section. If genes are not get cleaned, such as RBM|123 cannot be found which RBM is supposed to be in hg38 database.", style="color: #D5A200; font-size: 12px")
                                                        )
                                                      ),
                                                      
                                                      
                                                      h4("Data Preview", style="color: STEELBLUE; padding-top: 10px"),
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
                                               tabPanel("Eigengene Matrix",
                                                        tableOutput("mytable7"),
                                                        
                                                        h4("Survival Analysis", style="color: STEELBLUE; padding-top: 10px"),
                                                        helpText("Please select which row of the eigengene matrix would be applied for survival analysis. Groups are dichotomized by its median value."),
                                                        uiOutput("eigengene_matrix_select_row_ui"),
                                                        helpText("Please copy and paste following information in the correct order with regard to sample IDs (column names in above eigengene matrix), Note: separator (space, comma, new line, tab, or semicolon) will be identified automatically."),
                                                        textAreaInput("survival_event", "OS/EFS Events (1: deceased; 0: censored)", value = "", width = '100%', height = "20%", placeholder = "e.g., 0, 1, 0, ..."),
                                                        textAreaInput("survival_time", "OS/EFS Times", value = "", width = '100%', height = "20%", placeholder = "e.g., 12.3, 10.6, 5.0, ..."),
                                                        actionButton("run_survival_analysis", "Confirm and Run"),
                                                        uiOutput("survival_analysis_results_ui")
                                               ),
                                               tabPanel("Circos Plots",
                                                        
                                                        h4("Circos Plot", style="color: STEELBLUE"),
                                                        uiOutput("circos_plot_ui_hg38"),
                                                        uiOutput("circos_plot_ui_hg19"),
                                                        # plotOutput("circos_plot_component_hg38", width = "500px"),
                                                        # plotOutput("circos_plot_component_hg19", width = "500px"),
                                                        
                                                        h4("Parameters of Circos Plot", style="color: STEELBLUE"),
                                                        fluidRow(
                                                          column(6,
                                                                 textAreaInput("textareainput_circos", "Gene Symbols",
                                                                               value = "NKX2-5\nMEF2A\nGATA4\nHAND1\nHAND2\nTBX5\nSRF",
                                                                               width = 'auto', height = '300px', placeholder = NULL)
                                                          ),
                                                          column(6,
                                                                 sliderInput("circos_param_size", "Plot Size:",
                                                                             min = 100, max = 2000,
                                                                             value = 500),
                                                                 fluidRow(
                                                                   column(6, checkboxInput("circos_param_genelink", "Show Gene Links", TRUE)),
                                                                   column(6, checkboxInput("circos_param_genesymbol", "Show Gene Symbols", TRUE))
                                                                 )
                                                          )
                                                        ),
                                                        actionButton("circos_button_update_1", "Update Plots",style="color: WHITE; background-color: DODGERBLUE"),
                                                        helpText("You can directly use your own data here without any previous operation.")
                                                        
                                               )
                                             ),
                                             
                                             tags$script("$(document).on('click', '#clusterResult button', function () {
                                                         Shiny.onInputChange('button_lastClickId',this.id)
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
                                                        a("Yoav Benjamini and Yosef Hochberg. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society. Series B (Methodological)
                                                          Vol. 57, No. 1 (1995), pp. 289-300", href="http://www.jstor.org/stable/2346101", target="_blank"),
                                                        h5("Relationship between P-value, Z-score, and combined score:"),
                                                        helpText("The combined score is a combination of the p-value and z-score calculated by multiplying the two scores as follows:
                                                                 c = ln(p) * z
                                                                 Where c is the combined score, p is the p-value computed using Fisher's exact test, and z is the z-score computed to assess the deviation from the expected rank. The combined score provides a compromise between both methods and in several benchmarks we show that it reports the best rankings when compared with the other scoring schemes."),
                                                        
                                                        
                                                        tabsetPanel(
                                                          id = 'tabset',
                                                          tabPanel("GO_Biological_Process_2018",
                                                                   DT::dataTableOutput("mytable_Enrichr_1")),
                                                          tabPanel("GO_Molecular_Function_2018",
                                                                   DT::dataTableOutput("mytable_Enrichr_2")),
                                                          tabPanel("GO_Cellular_Component_2018",
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
                                             
                                             h5(a("ToppGene", href="https://toppgene.cchmc.org/enrichment.jsp", target="_blank")),
                                             h5(a("DAVID", href="https://david.ncifcrf.gov/tools.jsp", target="_blank")),
                                             h5(a("Enrichr", href="http://amp.pharm.mssm.edu/Enrichr/", target="_blank"))
                                             
                                           )
                                  )
)  #, style='width: 80%'

),
tabPanel("Tutorial",
         h3("Tutorial", style="color: STEELBLUE; padding-bottom: 20px"),
         h4("Google Slides", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
         tags$div(
           HTML("<iframe src=\"https://docs.google.com/presentation/d/e/2PACX-1vRoTd-UJbIyX6JEraVA1jZw_BWunHm_3qiq3qtxFE3y2DXjopc-SsnGmPXIalOof4iRX5d2eMfqGbji/embed?start=false&loop=false&delayms=3000\" frameborder=\"0\" width=\"960\" height=\"600\" allowfullscreen=\"true\" mozallowfullscreen=\"true\" webkitallowfullscreen=\"true\"></iframe>"),
           style="text-align: center; padding: 20px"
         ),
         h3("Presentation", style="color: STEELBLUE; padding-bottom: 20px"),
         h4("Google Slides", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
         tags$div(
           HTML("<iframe src=\"https://docs.google.com/presentation/d/e/2PACX-1vSv-c_P5P2dFCo9oP67JeBWRIrjZkxLgEkytC6edxUps7l4udMdWHqZx9kiltOwlIoWyWgJH-yDPqJY/embed?start=false&loop=false&delayms=3000\" frameborder=\"0\" width=\"960\" height=\"749\" allowfullscreen=\"true\" mozallowfullscreen=\"true\" webkitallowfullscreen=\"true\"></iframe>"),
           style="text-align: center; padding: 20px"
         ),
         # h4("Video Tutorial", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
         # tags$div(
         #   HTML('<iframe width="720" height="480" src="https://www.youtube.com/embed/d3eDKqm5yA0" frameborder="0" allowfullscreen></iframe>'),
         #   style="text-align: center; padding: 20px"
         # ),
         h4("Github", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
         tags$div(
           a("https://github.com/huangzhii/TSUNAMI", href="https://github.com/huangzhii/TSUNAMI"),
           style="text-align: center; padding: 0px"
         ),
         h4("Report Bugs", style="text-align: center; color: STEELBLUE; padding-bottom: 20px"),
         
         tags$div(
           a("https://github.com/huangzhii/TSUNAMI/issues/", href="https://github.com/huangzhii/TSUNAMI/issues/"),
           style="text-align: center; padding: 0px"
         )
),
tabPanel("FAQ",
         h3("Frequently Asked Questions", style="color: STEELBLUE; padding-bottom: 20px"),
         h4("General Questions", style="color: STEELBLUE; padding-bottom: 20px"),
         h5("What is the TSUNAMI website?", style="color: STEELBLUE"),
         p("The TSUNAMI (Translational Bioinformatics Tool Suite for Network Analysis and Mining) was developed at Indiana University School of Medicine."),
         h5("How do I get started?", style="color: STEELBLUE"),
         p("Please refer to our tutorial."),
         h5("Can I use TSUNAMI to analysis my data from my mobile devices?", style="color: STEELBLUE"),
         p("Yes you can!"),
         tags$div(
           tags$img(src='images/iphonex.png',
                    width="200",
                    alt="TSUNAMI_iPhoneX", style="padding: 0px"),
           style="padding: 10px"
         ),
         p("TSUNAMI website adopted responsive web design and is compatible with any mobile terminal. Every process, analysis, and computation is performed on the server behind your mobile browser. File uploading system would still work even on the phone when you are waiting a bus.")
),
tabPanel("News",
         h3("News", style="color: STEELBLUE; padding-bottom: 20px"),
         h4("April 24, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("Updated GEO offline data list to date 04/24/2018."),
           tags$li("Fixed a bug when percentiles are 0 or NULL."),
           tags$li("Moved platform converter to the right siderbar."),
           tags$li("Add conditional Panel on Advanced Data Preprocessing.")
         ),
         h4("March 27, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("Texts are modified."),
           tags$li("Footer added."),
           tags$li("Update pipeline flowchart."),
           tags$li("Update funding information."),
           tags$li("Added a Google Slides tutorial talk.")
         ),
         h4("March 20, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("R package 'lmQCM' was released to CRAN."),
           tags$li("Create flowchart.")
         ),
         h4("March 16, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("Renamed our website as TSUNAMI."),
           tags$li("Various of bugs are fixed.")
         ),
         h4("March 02, 2018", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("First prototype platform has been deployed.")
         )
),
tabPanel("About",
         h3("About Us", style="color: STEELBLUE; padding-bottom: 20px"),
         "The TSUNAMI (Translational Bioinformatics Tool Suite for Network Analysis and Mining) was developed at Indiana University School of Medicine.",
         "The design of such user-friendly implementations of our TSUNAMI pipeline provides a comprehensive analysis tool suite for users to study gene interaction from raw transcriptomic data level all the way to the gene ontology level with just simple button clicks.",
         tags$div(
           tags$img(src='images/IUSM2.png',
                    height="100",
                    alt="TSUNAMI", class="center", style="padding: 30px"),
           tags$img(src='images/regenstrief.png',
                    height="100",
                    alt="TSUNAMI", class="center", style="padding: 30px"),
           style="text-align: center; padding: 20px"
         ),
         h4("Our Other Softwares", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$div(
           a(tags$img(src='images/lmQCM_logo.png',
                      height="60",
                      alt="lmQCM", class="center", style="padding: 5px"), href="https://CRAN.R-project.org/package=lmQCM", target="_blank"),
           br(),a("R package: lmQCM", href="https://CRAN.R-project.org/package=lmQCM", target="_blank"),
           br(),br(),
           a(tags$img(src='images/annoPeak_logo.png',
                      height="40",
                      alt="annoPeak", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/content/19/", target="_blank"),
           br(),a("annoPeakR: a web-tool to annotate, visualize and compare peak sets from ChIP-seq/ChIP-exo", href="https://apps.medgen.iupui.edu/rsc/content/19/", target="_blank"),
           style="text-align: center; padding: 5px"
         ),
         br(),
         tags$div(
           a(tags$img(src='images/iGenomicsR_logo2.png',
                      height="40",
                      alt="iGenomicsR", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/content/27/", target="_blank"),
           br(),a("iGenomicsR: An integrative platform to explore, visualize, and analyze multidimensional genomics data for disease", href="https://apps.medgen.iupui.edu/rsc/content/27/", target="_blank"),
           br(),br(),
           a(tags$img(src='https://apps.medgen.iupui.edu/rsc/content/25/images/circosviewer_logo.png',
                      height="50",
                      alt="Circos Viewer", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/content/25/", target="_blank"),
           br(),a("Circos Viewer: A Circos Plot Viewer.", href="https://apps.medgen.iupui.edu/rsc/content/25/", target="_blank"),
           br(),br(),
           a(tags$img(src='https://apps.medgen.iupui.edu/rsc/content/23/_w_f03c09b9/images/iGPSeplus_logo.png',
                      height="40",
                      alt="iGPSe Plus", class="center", style="padding: 5px"), href="https://apps.medgen.iupui.edu/rsc/content/23/", target="_blank"),
           br(),a("iGPSe Plus: Integrative Genomic based Canser Patient Stratification", href="https://apps.medgen.iupui.edu/rsc/content/23/", target="_blank"),
           style="text-align: center; padding: 5px"
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
         tags$ul(
           tags$li("-")
         ),
         h4("Funding for the TSUNAMI is or has been provided by:", style="color: STEELBLUE; padding-bottom: 20px"),
         tags$ul(
           tags$li("Partially supported by IUSM startup fund, the NCI ITCR U01 (CA188547)."),
           tags$li("Data Science and Bioinformatics Program for Precision Health Initiative, Indiana University.")
         )
         
),


# navbarMenu(
#   "More",
#   tabPanel("Developer",
#            h4("Author Information"),
#            helpText("Indiana University School of Medicine"),
#            h4("Publication"),
#            helpText("Please cite ...")
#            )
# )
tags$div(
  p(a("TSUNAMI", href="https://apps.medgen.iupui.edu/rsc/tsunami"), "Version v1.8 | ", a("IUSM",href="https://medicine.iu.edu/", target="_blank"), " | ", a("RI",href="http://www.regenstrief.org/", target="_blank"), style="color: grey; font-size: 12px"), 
  p("Questions and feedback: zhihuan@iu.edu | ", a("Report Issue", href="https://github.com/huangzhii/TSUNAMI/issues", target="_blank"), " | ", a("Github", href="https://github.com/huangzhii/TSUNAMI/", target="_blank"), style="color: grey; font-size: 12px"),
  style="text-align: center; padding-top: 40px"
)
)
