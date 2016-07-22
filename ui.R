library(shiny)

shinyUI(fluidPage(
  titlePanel("CCBR Microarray analysis workflow", windowTitle="CCBR Microarray analysis workflow"),
  h5("(For human and mouse data)"),
  sidebarLayout(
    sidebarPanel(
      # tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
                width=12,
                 fluidRow(align="Top",
                          
                          column(3,
                                 textInput("ProjectID", label=h6("Project ID:"), value="ccbr", width="150px"),
                                 radioButtons("Platform", label=h6("Select platform"), choices=c("hgu133plus2" = "h133p2", "Mouse.gene.2.0.st" = "mst2"),selected="mst2")
                          ),
                          column(3,
                                 # textInput("Indir", label=h6("Path to input directory"), value="/Users/elloumif/Documents/cels6/cels",width="300px")
                                 fileInput("Indir", label=h6("Select CEL files"),multiple =T)
                                 
                          ),
                          
                          column(3,
                                 fileInput("pheno", label=h6("Choose phenotype file"))
                          ),
                          column(3,
                                 fileInput("const", label=h6("Choose contrast file")),
                                 numericInput("NumContrasts", label=h6("Which contrast to show"),value="1", width="150px"),
                                 numericInput("pval", label=h6("Contrast Pvalue threshold for Pathway/GO analysis"),value="0.05", width="150px")
                          )
                          #       numericInput("NumContrasts", label=h6("Which contrast to show"),value="1", width="150px"),
                           #      numericInput("fdr", label=h6("FDR for Pathway enrichment"),value="0.05", width="150px")
                          #)
                 ),
                
            
                # submitButton(text="Click to assign samples to groups and create contrasts")
                actionButton(inputId="go",label="Start"),
                actionButton(inputId="rep",label="Generate Report"),
                downloadButton('downloadReport', label = 'Download Report'),
                downloadButton('downloadTables', label = 'Download Tables')
    ),
    mainPanel(
      navbarPage( title = "Microarray",
                  tabPanel("Results"," "),
                  navbarMenu (title="Pre-normalization QC plots",
                    tabPanel("Histogram",plotOutput("rawhist")),
                    tabPanel("Maplots",uiOutput("rawmaplot")),
                    tabPanel("Boxplots", plotOutput("rawbox")),
                    tabPanel("RLE",plotOutput("rle")),
                    tabPanel("NUSE",plotOutput("nuse"))
                  ),
                  navbarMenu (title="Post-normalization plots",
                              tabPanel("Histogram",plotOutput("rmahist")),
                              tabPanel("Maplots",uiOutput("normaplot")),
                              tabPanel("Boxplots",plotOutput("rmabox")),
                              tabPanel("PCA2D",plotOutput("pca2d")),
                              tabPanel("Heatmap",plotOutput("heatmap"))
                              
                  ),
                  navbarMenu (title="DEG-Enrichments-tables",
                              tabPanel("DEG",DT::dataTableOutput("deg")),
                              tabPanel("KEGG",DT::dataTableOutput("kegg")),
                              tabPanel("GO",DT::dataTableOutput("go"))
                  ),
                  navbarMenu (title="Help",
                              tabPanel("Manual",htmlOutput("manu"))
                               
                  )
                  
      )
      # end solution2
    )
  )
  ## div(style="display:inline-block",submitButton("Generate PDF report"), style="float:center")
))
