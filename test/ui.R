navbarPage("My Application",
                   
                   tabPanel
                   (
                     "Select Data range",
                     sidebarLayout
                     (
                       sidebarPanel
                       (tags$head(tags$script('
                                              Shiny.addCustomMessageHandler("myCallbackHandler",
                                              function(typeMessage) {console.log(typeMessage)
                                              if(typeMessage == 1){
                                              console.log("got here");
                                              $("a:contains(Select Resolution)").click();
                                              }
                                              if(typeMessage == 2){
                                              $("a:contains(Select Data range)").click();
                                              }
                                              });
                                              ')),
                         h3("Select Data Range"),
                         selectInput("select", label = h3("Select Sector"),choices = list("Sector 1" = 1, "Sector 2" = 2,"Sector 3" = 3), selected = 1),br(),
                         dateRangeInput("dates", label = h3("Select Date range")),br(),
                         actionButton("action", label = "Proceed to select resolution")
                       ),
                       mainPanel("Output")
                       )
                     ),
                   
                   tabPanel
                   (
                     "Select Resolution",
                     sidebarLayout
                     (
                       sidebarPanel
                       (
                         h3("Select Resolution"),
                         numericInput("num1", label = h3("Select X-Grid Size"), value = 2),br(),
                         numericInput("num2", label = h3("Select Y-Grid Size"), value = 2),br(),
                         numericInput("num3", label = h3("Outlier Removal"), value = 2),br(),
                         numericInput("num4", label = h3("Frequency"), value = 2),br(),
                         actionButton("action1", label = "Proceed to Service Parameters")
                         
                       ),
                       mainPanel("Output"),
                       
                     )
                   )
)