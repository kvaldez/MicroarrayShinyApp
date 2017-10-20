library(shiny)
library(shinyjs)
library(GEOquery)
library(pd.mogene.2.0.st)
library(mogene20sttranscriptcluster.db)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(pd.hugene.2.0.st)
library(hugene20sttranscriptcluster.db)
library(pd.clariom.s.human.ht)
library(clariomshumanhttranscriptcluster.db)
library(pd.clariom.s.human)
library(clariomshumantranscriptcluster.db)
library(pd.clariom.s.mouse.ht)
library(clariomsmousehttranscriptcluster.db)
library(pd.clariom.s.mouse)
library(clariomsmousetranscriptcluster.db)
library(pd.mouse430.2)
library(mouse4302.db)
library(pd.hg.u133a)
library(hgu133a.db)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(pd.hg.u133a.2)
library(hgu133a2.db)
library(pd.huex.1.0.st.v2)
library(huex10sttranscriptcluster.db)
library(pd.hg.u219)
library(hgu219.db)
library(pd.mg.u74av2)
library(mgu74av2.db)
library(pd.mouse430a.2)
library(mouse430a2.db)
library(pd.moe430a)
library(moe430a.db)
library(pd.hg.u95av2)
library(hgu95av2.db)
library(pd.hta.2.0)
library(hta20transcriptcluster.db)
library(pd.moex.1.0.st.v1)
library(moex10sttranscriptcluster.db)
library(pd.hg.u133b)
library(hgu133b.db)
#library(GSEA)
library(limma)
library(oligo)
library(gplots)
library(geneplotter)
library(multtest)
library(rgl)
library(rglwidget)
library(DT)
library(getopt)
library(annotate)
library(knitr)
library(reshape)
library(RColorBrewer)
library(mixOmics)
library(calibrate)
library(rmarkdown)
library(ggplot2)
library(ggfortify)
library(shinyRGL)
library(plotly)
library(htmltools)
library(heatmaply)
library(Biobase)

shinyServer(function(input, output, session) {
  
  v <- reactiveValues(data = NULL, platform=NULL)
  
  observeEvent(
    input$button, 
    isolate({
      shinyjs::show("hide1")
    })
  )
  
  observeEvent(
    input$button,
    isolate({
      shinyjs::show("hide2")
    })
  )
  
  observeEvent(
    input$button2,
    isolate({
      shinyjs:: show("hide3")
    })
  )
  
  observeEvent(
    input$test,
    isolate({
      shinyjs:: show("hideAnalysis")
    })
  )
  
  
  observeEvent(
    input$button, 
    isolate({
      withProgress(message = 'Loading files...', value = 0.25, {
        if (input$gseid=="8 digit GSE code") return()
        id=input$gseid
        id = gsub(" ","",id,fixed=TRUE)       
        
        incProgress(0.25)
        
        #gds <- getGEO(input$gseid, GSEMatrix = F,getGPL=T,AnnotGPL=T)
        gds <- getGEO(id, GSEMatrix = F,getGPL=T,AnnotGPL=T)
        
        gds
        
        mytable=matrix("",length(GSMList(gds)),3)
        colnames(mytable)=c("gsm","title","description")
        
        for (k in 1:length(GSMList(gds)))
        {
          mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession, Meta(GSMList(gds)[[k]])$title, Meta(GSMList(gds)[[k]])$description)
        }
        
        mytable <- data.frame(mytable)
        
        mytable$group <- ""
        v$data <- mytable
        
        incProgress(0.25)
        output$mytable = DT::renderDataTable({
          if (is.null(v$data)) return()
          if (is.null(v$platform)) warning()
          DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8), pageLength = 8))
      })
      })
    })
  )
  
  observeEvent(
    input$button2,
    isolate({
      
      if (is.null(v$data)) return()
      if (is.null(input$group1)) return()
      
      #print(v$data)
      #print( "Error??2")

      v$data[input$mytable_rows_selected, "group" ] <- input$group1
      #print( "Error??3")
      
      output$mytable = DT::renderDataTable({
        
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8,10), pageLength = 8))
      })
    })
  )
  observeEvent(input$number, {
               isolate({
                 output$ui <- renderUI({
                   if (is.null(input$number))
                     return()

                   switch(input$number,
                          "1" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1"),
                                            selected = "Group 1"
                          ),
                          "2" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2"),
                                            selected = "Group 1"
                          ),
                          "3" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3"),
                                            selected = "Group 1"
                          ),
                          "4" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4"),
                                            selected = "Group 1"
                          ),
                          "5" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5"),
                                            selected = "Group 1"
                          ),
                          "6" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5",
                                                        "Group 6" = "Group 6"),
                                            selected = "Group 1"
                          ),
                          "7" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5",
                                                        "Group 6" = "Group 6",
                                                        "Group 7" = "Group 7"),
                                            selected = "Group 1"
                          ),
                          "8" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5",
                                                        "Group 6" = "Group 6",
                                                        "Group 7" = "Group 7",
                                                        "Group 8" = "Group 8"),
                                            selected = "Group 1"
                          ),
                          "9" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5",
                                                        "Group 6" = "Group 6",
                                                        "Group 7" = "Group 7",
                                                        "Group 8" = "Group 8",
                                                        "Group 9" = "Group 9"),
                                            selected = "Group 1"
                          ),
                          "10" = selectInput("group1", "Please select a group",
                                            choices = c("Group 1" = "Group 1",
                                                        "Group 2" = "Group 2",
                                                        "Group 3" = "Group 3",
                                                        "Group 4" = "Group 4",
                                                        "Group 5" = "Group 5",
                                                        "Group 6" = "Group 6",
                                                        "Group 7" = "Group 7",
                                                        "Group 8" = "Group 8",
                                                        "Group 9" = "Group 9",
                                                        "Group 10" = "Group 10"),
                                            selected = "Group 1"
                          )
                   )
                 })
               })
         })

  observe({
    if (input$test == 0)
      return()
    isolate({
      output$value <-renderTable({
        side1 <- paste0(switch(input$number,
                               "1" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1"),
                                                 selected = "Group 1"
                               ),
                               "2" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2"),
                                                 selected = "Group 1"
                               ),
                               "3" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3"),
                                                 selected = "Group 1"
                               ),
                               "4" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4"),
                                                 selected = "Group 1"
                               ),
                               "5" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4",
                                                             "Group 5" = "Group 5"),
                                                 selected = "Group 1"
                               ),
                               "6" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4",
                                                             "Group 5" = "Group 5",
                                                             "Group 6" = "Group 6"),
                                                 selected = "Group 1"
                               ),
                               "7" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4",
                                                             "Group 5" = "Group 5",
                                                             "Group 6" = "Group 6",
                                                             "Group 7" = "Group 7"),
                                                 selected = "Group 1"
                               ),
                               "8" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4",
                                                             "Group 5" = "Group 5",
                                                             "Group 6" = "Group 6",
                                                             "Group 7" = "Group 7",
                                                             "Group 8" = "Group 8"),
                                                 selected = "Group 1"
                               ),
                               "9" = selectInput("group1", NULL,
                                                 choices = c("Group 1" = "Group 1",
                                                             "Group 2" = "Group 2",
                                                             "Group 3" = "Group 3",
                                                             "Group 4" = "Group 4",
                                                             "Group 5" = "Group 5",
                                                             "Group 6" = "Group 6",
                                                             "Group 7" = "Group 7",
                                                             "Group 8" = "Group 8",
                                                             "Group 9" = "Group 9"),
                                                 selected = "Group 1"
                               ),
                               "10" = selectInput("group1", NULL,
                                                  choices = c("Group 1" = "Group 1",
                                                              "Group 2" = "Group 2",
                                                              "Group 3" = "Group 3",
                                                              "Group 4" = "Group 4",
                                                              "Group 5" = "Group 5",
                                                              "Group 6" = "Group 6",
                                                              "Group 7" = "Group 7",
                                                              "Group 8" = "Group 8",
                                                              "Group 9" = "Group 9",
                                                              "Group 10" = "Group 10"),
                                                  selected = "Group 1"
                               
        )
                               
        ), TRUE:input$test)
        vs.side2 <- paste0(switch(input$number,
                                  "1" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1"),
                                                    selected = "Group 1"
                                  ),
                                  "2" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2"),
                                                    selected = "Group 1"
                                  ),
                                  "3" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3"),
                                                    selected = "Group 1"
                                  ),
                                  "4" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4"),
                                                    selected = "Group 1"
                                  ),
                                  "5" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4",
                                                                "Group 5" = "Group 5"),
                                                    selected = "Group 1"
                                  ),
                                  "6" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4",
                                                                "Group 5" = "Group 5",
                                                                "Group 6" = "Group 6"),
                                                    selected = "Group 1"
                                  ),
                                  "7" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4",
                                                                "Group 5" = "Group 5",
                                                                "Group 6" = "Group 6",
                                                                "Group 7" = "Group 7"),
                                                    selected = "Group 1"
                                  ),
                                  "8" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4",
                                                                "Group 5" = "Group 5",
                                                                "Group 6" = "Group 6",
                                                                "Group 7" = "Group 7",
                                                                "Group 8" = "Group 8"),
                                                    selected = "Group 1"
                                  ),
                                  "9" = selectInput("group1", NULL,
                                                    choices = c("Group 1" = "Group 1",
                                                                "Group 2" = "Group 2",
                                                                "Group 3" = "Group 3",
                                                                "Group 4" = "Group 4",
                                                                "Group 5" = "Group 5",
                                                                "Group 6" = "Group 6",
                                                                "Group 7" = "Group 7",
                                                                "Group 8" = "Group 8",
                                                                "Group 9" = "Group 9"),
                                                    selected = "Group 1"
                                  ),
                                  "10" = selectInput("group1", NULL,
                                                     choices = c("Group 1" = "Group 1",
                                                                 "Group 2" = "Group 2",
                                                                 "Group 3" = "Group 3",
                                                                 "Group 4" = "Group 4",
                                                                 "Group 5" = "Group 5",
                                                                 "Group 6" = "Group 6",
                                                                 "Group 7" = "Group 7",
                                                                 "Group 8" = "Group 8",
                                                                 "Group 9" = "Group 9",
                                                                 "Group 10" = "Group 10"),
                                                     selected = "Group 1"
        )
        
        ), 
        TRUE:input$test)
        data.frame(side1, vs.side2)
        
      }, sanitize.text.function = function(x) x)
    })
  })
  observeEvent(input$analyze, {
      raw=reactive(input$analyze,
        {
      withProgress(message = 'Loading files...', value = 0.25, {
      id=input$gseid
      id = gsub(" ","",id,fixed=TRUE) 
      
      supp = getGEOSuppFiles(id, makeDirectory = T, baseDir = getwd())
      fileID = paste0(id, '_RAW.tar')
      #system(paste0('tar -xvf', fileID))
      untar(paste0(getwd(),'/',id,'/',fileID))
      incProgress(0.25)
      
      Pheno = v$data
      cels = paste0(Pheno$gsm,'_',Pheno$title,'.CEL.gz')   #adds filename
      rownames(Pheno) = cels
      incProgress(0.25)
      
      pd = AnnotatedDataFrame(Pheno)
      celfiles = read.celfiles(cels, phenoData = pd)
      cat(celfiles@annotation,file="annotation.txt")
      
      if (celfiles@annotation!="pd.hg.u133.plus.2" & celfiles@annotation!="pd.mogene.2.0.st" & celfiles@annotation!="pd.hugene.2.0.st" & celfiles@annotation!="pd.clariom.s.human.ht" & celfiles@annotation!="pd.clariom.s.human" & celfiles@annotation!="pd.clariom.s.mouse.ht" & celfiles@annotation!="pd.clariom.s.mouse" & celfiles@annotation!='pd.mouse430.2' & celfiles@annotation!='pd.hg.u133a' & celfiles@annotation!='pd.hugene.1.0.st.v1' & celfiles@annotation!='pd.mogene.1.0.st.v1' & celfiles@annotation!='pd.hg.u133a.2' & celfiles@annotation!='pd.huex.1.0.st.v2' & celfiles@annotation!='pd.hg.u219' & celfiles@annotation!='pd.mg.u74av2' & celfiles@annotation!='pd.mouse430a.2' & celfiles@annotation!='pd.moe430a' & celfiles@annotation!='pd.hg.u95av2' & celfiles@annotation!='pd.hta.2.0' & celfiles@annotation!='pd.moex.1.0.st.v1' & celfiles@annotation!='pd.hg.u133b') {
        #cat("Please sort your phenotype on sample name and upload it again. \n")
        info(paste0("Affymetrix platform: ",celfiles@annotation," NOT supported. Leaving..."))
        stopApp(-1)
      }
      incProgress(0.25)
      
      celfiles
      })
})
  })
  
  norm=reactive(
    {
      withProgress(message = 'Normalization', detail = 'starting ...', value = 0, {
        if (raw()@annotation=="pd.hg.u133.plus.2" | raw()@annotation=="pd.clariom.s.human.ht" | raw()@annotation=="pd.clariom.s.human" | raw()@annotation=="pd.clariom.s.mouse.ht" | raw()@annotation=="pd.clariom.s.mouse" | raw()@annotation=='pd.mouse430.2' | raw()@annotation=='pd.hg.u133a' | raw()@annotation=='pd.hg.u133a.2' | raw()@annotation=='pd.hg.u219' | raw()@annotation=='pd.mg.u74av2' | raw()@annotation=='pd.mouse430a.2' | raw()@annotation=='pd.moe430a' | raw()@annotation=='pd.hg.u95av2' | raw()@annotation=='pd.hg.u133b') {
          incProgress(0.5)
          celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL)
        } else {
          incProgress(0.5)
          celfiles.rma =rma(raw(), background=TRUE, normalize=TRUE, subset=NULL, target="core")
        }
      })
    })
  qc=reactive(
    {
      withProgress(message = 'Fitting probe level model', detail = 'starting ...', value = 0.75, {
        validate(
          need(raw()@annotation!= "pd.mogene.1.0.st.v1", 'NUSE and RLE plots unavailable for this platform.')
        )
        celfiles.qc=fitProbeLevelModel(raw())
      })
    }
  )
})


