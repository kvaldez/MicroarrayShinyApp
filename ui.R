library(shiny)
library(DT)
library(Biobase)
library(shinyjs)

shinyUI(
  fluidPage(
    useShinyjs(),
    
    
    headerPanel("GSE Table"),
    headerPanel(h5("made on shinyR")),
    textInput("gseid", label= h3("Accession Code"), value="8 digit GSE code"),
    helpText("*Make sure to input values correctly or the page will refresh."),
    helpText("*Data takes about 10 seconds to load"),
    actionButton(inputId="button", label="Display"),
    
    br(),
    br(),
    

    shinyjs::hidden(
    div(id= "hide1",
        fluidRow(
      
      column(3,
        numericInput("number", "Group amount:",
                     c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
        )
      )),
    
      column(3, 
        uiOutput("ui")
      ),
      
      br(),

    actionButton("button2", "Define")
    )),
  

    
    br(),
    shinyjs::hidden(
      div(id= "hide2",
        
    
    br(),
    br(),
    
   
          fluidRow(column(10, wellPanel(DT:: dataTableOutput("mytable"))))
        )),
    
    shinyjs:: hidden(
      div(id= "hide3",
    
    sidebarPanel(
      actionButton("test","Contrast")),
    mainPanel(
      uiOutput("value"))
        ))
  ))
  
