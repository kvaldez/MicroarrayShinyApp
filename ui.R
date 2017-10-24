library(shiny)
library(DT)
library(Biobase)
library(shinyjs)

shinyUI(
  fluidPage(
    tags$head(
      tags$style(
        HTML(".shiny-notification {
             height: 100px;
             width: 800px;
             position:fixed;
             font-weight: 500;
             font-size: 20px;
             background-color: #c1cdcd;
             top: calc(40% - 50px);;
             left: calc(50% - 400px);;
             }"
        )
      )
    ),
    useShinyjs(),
    
    titlePanel("CCBR Microarray analysis workflow", windowTitle="CCBR Microarray analysis workflow"),
    h5("(For Affymetrix human and mouse data)"),
    
    sidebarPanel(
      selectInput('analysisType', "Choose type of analysis",
                  c('Upload CEL files' = "CEL", 'Analyze GEO data' = 'GEO', ""),
                  selected = ''
      )
    ),
      
    conditionalPanel(
      condition = "input.analysisType == 'CEL'",
      fluidRow(),
      column(2,
        textInput("ProjectID", label=h6("Project ID:"), value="CCBR", width="150px")
      ),
      column(3,
        fileInput("Indir", label=h6("Select CEL files"),multiple =T)
      )
    ),
    conditionalPanel(
      condition = "input.analysisType == 'GEO'",
      fluidRow(),
      column(2,
        textInput("ProjectID", label=h6("Project ID:"), value="CCBR", width="150px")
      ),
      column(2,
        textInput("gseid", label= h6("Accession Code"), value="8 digit GSE code", width="150px")
        ),
      br(),
      br(),
      actionButton(inputId="button", label="Display")
    ),
    
    
    br(),
    br(),
    

    shinyjs::hidden(
      div(id= "hide1",
          fluidRow(
          column(3,
                 selectInput("number", "Number of groups:",
                             c("1"="1", "2"="2", "3"="3", "4"="4", "5"="5", "6"="6", "7"="7", "8"="8", "9"="9", "10"="10", ""),
                             selected=""
                 )
                 #actionButton("button3", "Enter")
          )),
          column(3, 
                 uiOutput("ui")
          ),
          br(), 
          fluidRow(
            column(3,
          actionButton("button2", "Define")
            )
      ))
    ),
    br(),
    
    
    #shinyjs::hidden(
    #  div(id='hide1',
      # conditionalPanel(
      # condition = 'input.number == 1',
      #     fluidRow(),
      #       column(2,
      #         #lapply(1:input$number, function(i) {     
      #         selectInput("group1", "Please select a group",
      #                     #choices = paste0('Group', i),
      #                     choices = c("Group 1" = "Group 1"),
      #                     selected = "Group 1"
      #         #)}
      #       ),
      #       actionButton("button2", "Define")
      # )
      # ),
      # 
    #   conditionalPanel(
    #   condition = 'input.number == 2',
    #   fluidRow(),
    #   column(2,
    #          #lapply(1:input$number, function(i) {     
    #          selectInput("group1", "Please select a group",
    #                      #choices = paste0('Group', i),
    #                      choices = c("Group 1" = "Group 1", "Group 2" = "Group 2"),
    #                      selected = "Group 1"
    #                      #)}
    #          ),
    #          actionButton("button2", "Define")
    #   )
    # ),
    
      
    br(),
    shinyjs::hidden(
      div(id= "hide2",
        
    
    br(),
    br(),

          fluidRow(column(10, wellPanel(DT:: dataTableOutput("mytable"))))
        )),
    
    shinyjs:: hidden(
      div(id= "hide3",
    mainPanel(
      actionButton("test","Contrast")),
    mainPanel(
      uiOutput("value"))
        )),
    
    shinyjs::hidden(
      div(id='hideAnalysis',
        mainPanel(
          actionButton(inputId="analyze",label="Start")
        )
    )),
    
    shinyjs::hidden(
    div(id='hideResults',
        mainPanel(
          navbarPage(title = "Results",
                     navbarMenu (title="Pre-normalization QC plots",
                                 tabPanel("Histogram",plotOutput("rawhist")),
                                 tabPanel("Maplots", uiOutput("rawmaplot")),
                                 tabPanel("Boxplots", plotOutput("rawbox")),
                                 tabPanel("RLE",plotOutput("rle")),
                                 tabPanel("NUSE",plotOutput("nuse"))
                     ),
                    navbarMenu (title="Post-normalization plots",
                                tabPanel("Histogram",plotOutput("rmahist")),
                                tabPanel("Maplots",uiOutput("normaplot")),
                                tabPanel("Boxplots",plotOutput("rmabox")),
                                tabPanel("3D-PCA",rglwidgetOutput("pca3d")),
                                tabPanel("Interactive Heatmap",plotlyOutput("heatmap"))
                    ),
                    navbarMenu (title="DEG-Enrichments-tables",
                                tabPanel("Differentially Expressed Genes",DT::dataTableOutput("deg"))
                    )
          )
        )
    )
    )
  )
)
 
 

  
  
