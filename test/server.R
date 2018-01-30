library(shiny)
shinyServer(function(input, output,session) {
  observe({
    if(input$action > 0){
      print('1')
      session$sendCustomMessage("myCallbackHandler", "1")
    }
  })
  observe({
    if(input$action1 > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler", "2")
    }
  })
}
)