library(shiny)
library(ggplot2)  # for the diamonds dataset
library(enrichR)
library(DT)
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
        tabPanel("mtcars", DT::dataTableOutput("mytable_GO_1")),
        tabPanel("mtcars", DT::dataTableOutput("mytable_GO_2")),
        tabPanel("mtcars", DT::dataTableOutput("mytable_GO_3")),
        tabPanel("iris", DT::dataTableOutput("mytable3"))
      )
    )
  )
)

server <- function(input, output) {

  # choose columns to display
  diamonds2 = diamonds[sample(nrow(diamonds), 1000), ]

  enrichr_dbs <- c("GO_Biological_Process_2017b", "GO_Molecular_Function_2017b", "GO_Cellular_Component_2017b")
  # print(geneCharVector_global[[cluster]])
  # genes_str <- geneCharVector_global[[cluster]]
  genes_str <- c('PHF14','RBM3','Nlrx1','MSL1','PHF21A','ARL10','INSR')
  enriched <- enrichr(genes_str[-1], enrichr_dbs)
  print(genes_str)

  output$mytable1 <- DT::renderDataTable({DT::datatable( enriched[["GO_Biological_Process_2017b"]], selection="none", escape=FALSE,
                                       options = list(paging = F, searching = F, autoWidth = TRUE, dom='t',ordering=F),
                                       rownames = F)
    })

  # sorted columns are colored now because CSS are attached to them
  output$mytable2 <- DT::renderDataTable({
    DT::datatable(mtcars, options = list(orderClasses = TRUE)) %>%
      formatStyle('mpg', backgroundColor = 'yellow') %>% formatRound(colnames(mtcars), digits=3)
  })
  cname <- letters[1:3]
  Map(function(id) {
    if (id==1){db = mtcars}
    if (id==2){db = mpg}
    if (id==3){db = iris}
    output[[paste("mytable_GO",id,sep="_")]]<- DT::renderDataTable({
      DT::datatable(db, options = list(orderClasses = TRUE))
    })
  },1:3)

  # customize the length drop-down menu; display 5 rows per page by default
  output$mytable3<- DT::renderDataTable({
    DT::datatable(data_final[, 1:10],
                  extensions = 'Responsive', escape=F, selection = 'none')
  })

}

shinyApp(ui, server)
