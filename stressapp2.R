# Load R packages
library(shiny)
library(shinythemes)
library(tableHTML)


# Read in the RF model
model <- readRDS("testmod.rds")


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
    textInput("ID", "Sample ID:", "1a"),
    numericInput("Elevation", 
                 label = "Elevation (in meters)", 
                 value = 1500),
    numericInput("tmean", 
                 label = "Mean Annual Temperature (°C)", 
                 value = 12),
    numericInput("PPT", 
                 label = "Mean Annual Precipitation (mm)", 
                 value = 600),
    numericInput("tdmean", 
                 label = "Mean Annual Dew Point Temperature (°C)", 
                 value = 0.4),
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
      Elevation = as.numeric(c(input$Elevation)),
      tmean = as.numeric(c(input$tmean)),
      PPT = as.numeric(c(input$PPT)),
      tdmean = as.numeric(c(input$tdmean)),
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
  output$info <- renderText({
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

rsconnect::setAccountInfo(name='zsnider',
                          token='79499AE0FBC07D122C3DC9299B930833',
                          secret='amR8Y6VGycgey1DDtyLWQAXR2iWukaAXz9xp8alM')

####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)