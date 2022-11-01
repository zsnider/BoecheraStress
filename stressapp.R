

# Load R packages
library(shiny)
library(shinythemes)


# Read in the RF model
model <- readRDS("solmod.rds")


####################################
# User interface                   #
####################################

ui <- pageWithSidebar(
  
  # Page header
  headerPanel('Boechera Stress Predictor'),
  
  # Input values
  sidebarPanel(
    #HTML("<h3>Input parameters</h3>"),
    tags$label(h3('Input parameters')),
    textInput("ID", "Sample ID:", "1A"),
    numericInput("solclear", 
                 label = "Solar Radiation", 
                 value = 20),
    actionButton("submitbutton", "Submit", 
                 class = "btn btn-primary")
  ),
  
  mainPanel(
    tags$label(h3('Status/Output')), # Status/Output Text Box
    verbatimTextOutput('contents'),
    textOutput('newdata'),
    tags$head(tags$style("#newdata{color: Black;
                                 font-size: 25px;
                                 font-style: bold;
                                 }"
    )
    ),
    textOutput('info')
  )
)

####################################
# Server                           #
####################################

server<- function(input, output, session) {
  
  # Input Data
  datasetInput <- reactive({  
    
    df <- data.frame(
      ID = as.character(c(input$ID)),
      solclear = as.numeric(c(input$solclear)),
      stringsAsFactors = FALSE)
    
    #Species <- 0
    #df <- rbind(df, Species)
    #input <- transpose(df)
    #write.table(input,"input.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #test <- read.csv(paste("input", ".csv", sep=""), header = TRUE)
    
    #Output <- data.frame(Prediction=predict(model,test), round(predict(model,test,type="prob"), 3))
    Output <- predict(model, df)
    print(Output)
    
  })
  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate("Calculation complete.") 
    } else {
      return("Server is ready for calculation.")
    }
  })
  
  #Info Text Box
  output$info <- renderPrint({
    if (input$submitbutton>0) { 
      isolate("This result is given in YII (quantum yield). 
              Plants that are more stressed will exhibit a lower YII.") 
    }
  })

  # Prediction results table
  output$newdata <- renderText({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } 
  })
  
}

####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)