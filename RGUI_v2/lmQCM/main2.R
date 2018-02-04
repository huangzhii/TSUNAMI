library(shiny)
library(rsconnect)
library(plyr)
library(data.table)


options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB
options(shiny.sanitize.errors = FALSE)
setwd("/Users/zhi/Desktop/GeneCoexpression/shiny"); #mac
source("GeneCoExpressionAnalysis.R")

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("lmQCM for mining the locally dense structures"),
  
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
      # Input: Select a method ----
      selectInput("method", "Choose a method:",
                  choices = c("lmQCM", "WGCNA")),
      
      # Input: Specify the number of observations to view ----
      fluidRow(
        column(6, numericInput("gamma", "gamma:", 0.55, step = 0.05)),
        column(6, numericInput("lambda", "lambda", 1, step = 0.05))
      ),
      fluidRow(
        column(6, numericInput("t", "t:", 1, step = 0.05)),
        column(6, numericInput("beta", "beta:", 0.4, step = 0.05))
      ),
      numericInput("minClusterSize", "Minimum Cluster Size:", 10, step = 1, width = NULL),
      
      
      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change). This is useful if
      # the computations required to render output are inordinately
      # time-consuming.
      actionButton("update", "Update View"),
      # Horizontal line ----
      tags$hr(),
      radioButtons("filetype", "Choose file type and download processed data:",
                   choices = c("csv", "txt")),
      downloadButton('downloadData', 'Download')
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  # tabPanel("Plot", plotOutput("plot")),
                  tabPanel("Log", verbatimTextOutput("log")),
                  tabPanel("Results", tableOutput("results"))
      )
      
    )
    
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  
  method <- eventReactive(input$update, {
    switch(input$method,
           "lmQCM" = "lmQCM",
           "WGCNA" = "WGCNA")
  }, ignoreNULL = FALSE)
  
  gamma <- eventReactive(input$update, {input$gamma}, ignoreNULL = FALSE)
  lambda <- eventReactive(input$update, {input$lambda}, ignoreNULL = FALSE)
  t <- eventReactive(input$update, {input$t}, ignoreNULL = FALSE)
  beta <- eventReactive(input$update, {input$beta}, ignoreNULL = FALSE)
  minClusterSize <- eventReactive(input$update, {input$minClusterSize}, ignoreNULL = FALSE)
  y <- reactive( if (dataset()) return(1) else return(0) )
  
  output$log <- renderText({
    return("Please upload a CSV file.")
  })
  dataset <- reactive({
    req(input$csvfile)
    print("CSV file Uploaded.")
    dataset <- read.csv(input$csvfile$datapath,
                     header = input$header,
                     sep = input$sep,
                     quote = input$quote)
    print("CSV file Readed.")
    return(dataset)
  })
  text <- reactive({
    req(dataset())
    text <- lmQCM(step1, dataset(), gamma(), lambda(), t(), beta(), minClusterSize())
    return(text)
  })
  
  output$results <- renderTable({
    text_multiline <- strsplit(text(), " ")
    out <- transpose(data.frame(lapply(text_multiline, `length<-`, max(lengths(text_multiline)))))
  },rownames = TRUE, colnames = FALSE, na = "", bordered = TRUE)
  # Show the first "n" observations ----
  # The use of isolate() is necessary because we don't want the table
  # to update whenever input$obs changes (only when the user clicks
  # the action button)
  
  # obs = 10
  # output$view <- renderTable({
  #   head(datasetInput(), n = isolate(input$obs))
  # })
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      name = paste("lmQCMresult","gamma",gamma(),"lambda",lambda(),"t",t(),"beta",beta(),"minClusterSize",minClusterSize(), sep = "_", collapse = NULL)
      paste(name, input$filetype, sep = ".")
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      sep <- switch(input$filetype, "csv" = 0, "txt" = 1)
      text = text()
      if (sep == 0){
        text = gsub(" ", ",", noquote(text))
        write.table(text, file, eol = "\r\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
      if (sep == 1){
        text = gsub(" ", "\t", noquote(text))
        write.table(text, file, eol = "\r\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)