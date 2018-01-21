library(shiny)
library(rsconnect)

options(shiny.maxRequestSize=300*1024^2) # to the top of server.R would increase the limit to 300MB


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("lmQCM for mining the locally dense structures"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      fileInput("csvfile", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Include clarifying text ----
      helpText("Note: Maximum csv file size allowed for uploading is 300MB."),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      # Input: Select a method ----
      selectInput("method", "Choose a method:",
                  choices = c("lmQCM", "WGCNA")),
      
      # Input: Specify the number of observations to view ----
      numericInput("gamma", "gamma:", 0.55),
      numericInput("lambda", "lambda", 1),
      numericInput("t", "t:", 1),
      numericInput("beta", "beta:", 0.4),
      numericInput("minClusterSize", "Minimum Cluster Size:", 10),
      
      
      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change). This is useful if
      # the computations required to render output are inordinately
      # time-consuming.
      actionButton("update", "Update View")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  # tabPanel("Plot", plotOutput("plot")),
                  tabPanel("Summary", verbatimTextOutput("summary")),
                  tabPanel("Table", tableOutput("table"))
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
  
  # data <- reactive({
  #   # file <- input$csvfile
  #   if(is.null(input$csvfile)) {return(0)}
  #   read.csv(input$csvfile$datapath,
  #                    header = input$header,
  #                    sep = input$sep,
  #                    quote = input$quote)
  # })
  
  gamma <- eventReactive(input$update, {
    input$gamma
  }, ignoreNULL = FALSE)
  lambda <- eventReactive(input$update, {
    input$lambda
  }, ignoreNULL = FALSE)
  t <- eventReactive(input$update, {
    input$t
  }, ignoreNULL = FALSE)
  beta <- eventReactive(input$update, {
    input$beta
  }, ignoreNULL = FALSE)
  minClusterSize <- eventReactive(input$update, {
    input$minClusterSize
  }, ignoreNULL = FALSE)
  
  # Generate a summary of the dataset ----
  output$summary <- renderPrint({
    req(input$csvfile)
    
    data <- read.csv(input$csvfile$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    print("CSV file Uploaded.")
    if (method() == "lmQCM"){
      source("GeneCoExpressionAnalysis.R")
      # Check if is a new data, then determine to perform lmQCM_1 or not.
      step1 <- 0
      if (length(data) >0){
        print("Processing CSV file.")
        output <- lmQCM(step1, data, gamma(), lambda(), t(), beta(), minClusterSize())
        # Run GeneCoExpressionAnaylsis
        # print(output)
        for (i in 1:(length(output))) {
          cat(as.matrix(output[[i]]), sep=' ')
          cat('\n')
        }
      }
      else {
        print("Please upload a CSV file.")
      }
    } else {
      method()
      # Run WGCNA
    }
  })
  output$table <- renderTable({
    if(is.null(data())) {return()}
    input$csvfile
  })
  # Show the first "n" observations ----
  # The use of isolate() is necessary because we don't want the table
  # to update whenever input$obs changes (only when the user clicks
  # the action button)
  
  # obs = 10
  # output$view <- renderTable({
  #   head(datasetInput(), n = isolate(input$obs))
  # })
  
}

# Create Shiny app ----
shinyApp(ui, server)