library(shiny)
library(ggplot2)  # for the diamonds dataset

ui <- fluidPage(
  title = "Examples of DataTables",
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.dataset === "diamonds"',
        checkboxGroupInput("show_vars", "Columns in diamonds to show:",
                           names(diamonds), selected = names(diamonds))
      ),
      conditionalPanel(
        'input.dataset === "mtcars"',
        helpText("Click the column header to sort a column.")
      ),
      conditionalPanel(
        'input.dataset === "iris"',
        helpText("Display 5 records by default.")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("diamonds", DT::dataTableOutput("mytable1")),
        tabPanel("mtcars", DT::dataTableOutput("mytable2")),
        tabPanel("iris", DT::dataTableOutput("mytable3"))
      )
    )
  )
)

server <- function(input, output) {

  # choose columns to display
  diamonds2 = diamonds[sample(nrow(diamonds), 1000), ]

  load(file='~/Desktop/enrichment_result.Rdata')

  for (i in 1:length(enrichment_result)){
    enrichment_result[[i]][[6]] = paste(enrichment_result[[i]][[6]], collapse = ',')
  }

  enrichment_result2 = unlist(enrichment_result)

  enrichment_result2 = data.frame(t(matrix(enrichment_result2, nrow = length(enrichment_result[[1]]), byrow = F)))
  print(enrichment_result2)
  output$mytable1 <- DT::renderDataTable({DT::datatable( enrichment_result2, selection="none", escape=FALSE,
                                       options = list(paging = F, searching = F, autoWidth = TRUE, dom='t',ordering=F),
                                       rownames = F)
    })

  # sorted columns are colored now because CSS are attached to them
  output$mytable2 <- DT::renderDataTable({
    DT::datatable(mtcars, options = list(orderClasses = TRUE))
  })

  # customize the length drop-down menu; display 5 rows per page by default
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(iris, options = list(lengthMenu = c(5, 30, 50), pageLength = 5))
  })

}

shinyApp(ui, server)
