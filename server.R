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

setwd('/Users/valdezkm/Documents/MicroarrayPipeline_VERSION2/MicroarrayShinyApp')

shinyServer(function(input, output) {
  
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
    input$analyze,
    isolate({
      shinyjs:: show("hideResults")
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
  observeEvent(
    input$analyze, {
      isolate({
      raw <- reactive(
        {
      withProgress(message = 'Loading files...', value = 0.25, {
      id=input$gseid
      id = gsub(" ","",id,fixed=TRUE) 
      
      system(paste0('rm *.CEL.gz'))        #removes previous CEL files
      getGEOSuppFiles(id, makeDirectory = T, baseDir = getwd())
      fileID = paste0(id, '_RAW.tar')
      #system(paste0('tar -xvf', fileID))
      untar(paste0(getwd(),'/',id,'/',fileID))
      incProgress(0.25)
 
      #cels = paste0(Pheno$gsm,'_',Pheno$title,'.CEL.gz')   #adds filename
      Pheno = v$data
      system(paste0('ls ',getwd(),'/*.gz > SampleName.txt'))    #list contents of new directory with zipped CEL files
      SampleName = read.delim('SampleName.txt', sep='\n', header = F)
      SampleName = basename(as.character(unlist(SampleName)))
      rownames(Pheno) = SampleName
      cels = SampleName
      
      if (length(grep('*CEL*',SampleName,ignore.case = T)) == 0) {
        info("Raw files must be CEL files")
      }
      
      incProgress(0.25)
      
      pd = AnnotatedDataFrame(Pheno)
      celfiles = read.celfiles(cels, phenoData = pd)
      colnames(pData(celfiles))[2] = 'SampleID'     
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
      # list of DEG
      deg=reactive(
        {
          ##-------------
          withProgress(message = 'Computing differentially expressed genes', value = 0, {
            facs <- factor(pData(raw())$SampleGroup)
            labfacs=levels(facs)
            nbfacs=length(labfacs)
            file1=input$const

            contra=read.delim(file1$datapath)
            nb=dim(contra)[1]
            cons=c()
            #validate(
            #  need((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs), "One of the groups in contrast file does not match a group in phenotype file. Make sure names match and upload again. 
            #Once correct file is entered, 'Computing differentially expressed genes' message will display.")
            #)
            
            for (k in 1:nb) {
              if ((contra[k,1] %in% labfacs) & (contra[k,2] %in% labfacs) )
              {
                cons=c(cons,paste(contra[k,1],"-",contra[k,2],sep=""))
              } else {
                #cat("One of the groups in contrasts file at line :",k+1,"does not match a group in phenotype file..Quitting!!!\n")
                info('One of the groups in contrast file does not match a group in phenotype file. Make sure names match and upload again.')
                print( contra )
                stopApp(-1)
              }
            }
            
            
            myfactor <- factor(pData(norm())$SampleGroup)
            design1 <- model.matrix(~0+myfactor)
            colnames(design1) <- levels(myfactor)
            
            fit1 <- lmFit(norm(),design1)
            contrast.matrix <- makeContrasts(contrasts=cons,levels=design1)
            
            fit2 <- contrasts.fit(fit1, contrast.matrix)
            ebayes.fit2=eBayes(fit2) # smooths the std error
            incProgress(0.25, detail = 'Limma model fitted')
            # #EXTRACTING ALL GENES FOR EACH CONTRAST
            
            ##ANNOTATE PROBESET IDS FROM ANNOTATION PACKAGE FROM BIOCONDUCTOR
            ## load libraries as sources of annotation
            
            #library(mogene20sttranscriptcluster.db)
            #if (input$Platform=="mst2") {
            #if (raw()@annotation=="pd.mogene.2.0.st") {  
            #  Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
            #} else {
            # if (input$Platform=="h133p2") {
            #   if (raw()@annotation=="pd.hg.u133.plus.2") {
            #     Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
            #  } 
            #} 
            
            if (raw()@annotation=="pd.mogene.2.0.st") {  
              Annot <- data.frame(ACCNUM=sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene20sttranscriptclusterGENENAME), paste, collapse=", "))
            } else {
              # if (input$Platform=="h133p2") {
              if (raw()@annotation=="pd.hg.u133.plus.2") {
                Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
              } else {
                if (raw()@annotation=="pd.hugene.2.0.st") {
                  Annot <- data.frame(ACCNUM=sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse=", "))
                } else {
                  if (raw()@annotation=="pd.clariom.s.human.ht") {
                    Annot <- data.frame(ACCNUM=sapply(contents(clariomshumanhttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumanhttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumanhttranscriptclusterGENENAME), paste, collapse=", "))
                  } else {
                    if (raw()@annotation=="pd.clariom.s.mouse.ht") {
                      Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousehttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousehttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousehttranscriptclusterGENENAME), paste, collapse=", "))
                    } else {
                      if (raw()@annotation=="pd.clariom.s.mouse") {
                        Annot <- data.frame(ACCNUM=sapply(contents(clariomsmousetranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomsmousetranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomsmousetranscriptclusterGENENAME), paste, collapse=", "))
                      } else {
                        if (raw()@annotation=="pd.clariom.s.human") {
                          Annot <- data.frame(ACCNUM=sapply(contents(clariomshumantranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(clariomshumantranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(clariomshumantranscriptclusterGENENAME), paste, collapse=", "))
                        } else {
                          if (raw()@annotation=="pd.mouse430.2") {
                            Annot <- data.frame(ACCNUM=sapply(contents(mouse4302ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse4302SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse4302GENENAME), paste, collapse=", "))
                          } else {
                            if (raw()@annotation=='pd.hg.u133a') {
                              Annot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))
                            } else {
                              if (raw()@annotation=='pd.hugene.1.0.st.v1') {
                                Annot <- data.frame(ACCNUM=sapply(contents(hugene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hugene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hugene10sttranscriptclusterGENENAME), paste, collapse=", "))
                              } else {
                                if (raw()@annotation=='pd.mogene.1.0.st.v1') {
                                  Annot <- data.frame(ACCNUM=sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mogene10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(mogene10sttranscriptclusterGENENAME), paste, collapse=", "))
                                } else {
                                  if (raw()@annotation=='pd.hg.u133a.2') {
                                    Annot <- data.frame(ACCNUM=sapply(contents(hgu133a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133a2GENENAME), paste, collapse=", "))
                                  } else {
                                    if (raw()@annotation=='pd.huex.1.0.st.v2') {
                                      Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "))
                                    } else {
                                      if (raw()@annotation=='pd.hg.u219') {
                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu219ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu219SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu219GENENAME), paste, collapse=", "))
                                      } else {
                                        if (raw()@annotation=='pd.ht.hg.u133.plus.pm') {
                                          Annot <- data.frame(ACCNUM=sapply(contents(hgu133plus2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133plus2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133plus2GENENAME), paste, collapse=", "))
                                        } else {
                                          if (raw()@annotation=='pd.mg.u74av2') {
                                            Annot <- data.frame(ACCNUM=sapply(contents(mgu74av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mgu74av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mgu74av2GENENAME), paste, collapse=", "))
                                          } else {
                                            if (raw()@annotation=='pd.mouse430a.2') {
                                              Annot <- data.frame(ACCNUM=sapply(contents(mouse430a2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(mouse430a2SYMBOL), paste, collapse=", "), DESC=sapply(contents(mouse430a2GENENAME), paste, collapse=", "))
                                            } else {
                                              if (raw()@annotation=='pd.moe430a') {
                                                Annot <- data.frame(ACCNUM=sapply(contents(moe430aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moe430aSYMBOL), paste, collapse=", "), DESC=sapply(contents(moe430aGENENAME), paste, collapse=", "))
                                              } else {
                                                if (raw()@annotation=='pd.hg.u95av2') {
                                                  Annot <- data.frame(ACCNUM=sapply(contents(hgu95av2ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu95av2SYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu95av2GENENAME), paste, collapse=", "))
                                                } else {
                                                  if (raw()@annotation=='pd.hta.2.0') {
                                                    Annot <- data.frame(ACCNUM=sapply(contents(hta20transcriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hta20transcriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(hta20transcriptclusterGENENAME), paste, collapse=", "))
                                                  } else {
                                                    if (raw()@annotation=='pd.moex.1.0.st.v1') {
                                                      Annot <- data.frame(ACCNUM=sapply(contents(moex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(moex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(moex10sttranscriptclusterGENENAME), paste, collapse=", "))
                                                    } else {
                                                      if (raw()@annotation=='pd.hg.u133b') {
                                                        Annot <- data.frame(ACCNUM=sapply(contents(hgu133bACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133bSYMBOL), paste, collapse=", "), DESC=sapply(contents(hgu133bGENENAME), paste, collapse=", "))
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            
            incProgress(0.25, detail = 'preparing for pathway analysis ...')
            mylist=vector("list",nb)
            
            for (i in 1:nb)
            {
              
              all.genes.con = topTable(ebayes.fit2, coef = i, number=nrow(ebayes.fit2))
              
              # Merge data frames together (like a database table join)
              
              all <- merge(all.genes.con, Annot,by.x=0, by.y=0, all.x=T)
              all=all[order(all$P.Value),]
              colnames(all)[1]="probsetID"
              
              #add fold change and rearrange columns
              all$FC = ifelse(all$logFC<0, -1/(2^all$logFC), 2^all$logFC)
              all = all[,c(9,1,8,10,11,2,5,6,3,4,7)]
              
              # Write out to a file:
              write.table(all,file=paste(input$ProjectID,"_",cons[i],"_all_genes.txt",sep=""),sep="\t",row.names=F)
              # cat("Contrast: ",i," done \n")
              
              mylist[[i]]=all
              ## end for
            }
            all <- merge(exprs(norm()), Annot,by.x=0, by.y=0, all.x=T)
            write.table(all,file=paste(input$ProjectID,"_normalized_data.txt",sep=""),sep="\t",row.names=F)
            #  
            names(mylist)=cons
            
            incProgress(0.5, detail = 'DEG done')
            
            #mylist
            list(mylist=mylist)
          })
          ##-------------
        }
      )
      
  #Processing all outputs
  
  ####creates a list of colors specific to each group
  fs = factor(pData(raw())$group)
  #fs = factor(pData(raw())$SampleGroup)
  lFs=levels(fs)
  numFs=length(lFs)
  colors = list()
  for (i in 1:numFs){
    colors[which(fs==lFs[i])] = i*5
  }
  colors = unlist(colors)
  ####end
  
  output$projectid=renderText({paste("Project ID: ",input$ProjectID)})
  output$rawhist=renderPlot(
    {
      hist(raw(),which="all", main =" Raw Samples distribution")
    }
  )
  
  ###Beginning raw maplot###
  output$rawmaplot=renderUI({
    #facs <- pData(raw())$SampleID
    facs <- pData(raw())$SampleID
    nbfacs=length(facs)
    plot_output_list <- lapply(1:nbfacs, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname, height = 400, width = 600)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  facs <- pData(raw())$SampleID
  nbfacs=length(facs)

  for (i in 1:nbfacs) {
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        withProgress(message = 'Generating Raw Maplot', detail = paste0('Plot ', my_i, ' ...'), value = (my_i/nbfacs), {
          MAplot(raw(),which=my_i,plotFun=smoothScatter,refSamples=c(1:nbfacs), main='', cex=2)
        })
      })
    })
  }
  ###end raw maplot###
  
  output$rawbox=renderPlot(
    {
      boxplot(raw(), col=colors, which="all", main="Boxplots before normalization",las=2,names=pData(raw())$SampleID)
    }
  )
  output$rle=renderPlot(
    {
      RLE(qc(), main="RLE plot",names=pData(raw())$SampleID, col=colors)
    }
  )
  output$nuse=renderPlot(
    {
      NUSE(qc(), main="NUSE plot",names=pData(raw())$SampleID, col=colors)
    }
  )
  output$rmahist=renderPlot(
    {
      hist(norm(), main ="Distribution after Normalization")
    }
  )
  
  ## MVAplot after normalization
  output$normaplot=renderUI({
    facs2 <- pData(norm())$SampleID
    nbfacs2=length(facs2)
    plot_output_list2 <- lapply(1:nbfacs2, function(i) {
      plotname2 <- paste("plota", i, sep="")
      plotOutput(plotname2, height = 400, width = 600)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list2)
  })
  
  facs2 <- pData(norm())$SampleID
  nbfacs2=length(facs2)
  
  for (i in 1:nbfacs2) {
    local({
      my_i <- i
      plotname2 <- paste("plota", my_i, sep="")
      # MA plots are then used to visualize intensity-dependent ratio for each group
      output[[plotname2]] <- renderPlot({
        withProgress(message = 'Generating Normalized Maplot', detail = paste0('Plot ', my_i, ' ...'), value = my_i/nbfacs2, {
          MAplot(norm(),which=my_i,plotFun=smoothScatter,refSamples=c(1:nbfacs2),main='', cex=2)
        })
      })
    })
  }
  output$rmabox=renderPlot(
    {
      boxplot(norm(),col=colors, main="Boxplots after RMA normalization",las=2,names=pData(norm())$SampleID)
    }
  )
  ## end mvaplat after normalization
  
  ## pca 3D
  output$pca3d=renderRglwidget(
    {
      withProgress(message = 'Generating PCA', detail = 'starting ...', value = 0.5, {
        tedf= t(exprs(norm()))
        
        #removes zero  variances (issue with small sample sizes)
        if (length(which(apply(tedf, 2, var)==0)) >= 0){
          tedf = tedf[ , apply(tedf, 2, var) != 0]
        }
        
        pca=prcomp(tedf, scale. = T)
        incProgress(amount = 0.25, detail = 'determining variance ...')
        rgl.open(useNULL=T)
        bg3d('white')
        plot3d(pca$x[,1:3],col=colors, type='s',size=2)
        group.v=as.vector(pData(norm())$SampleID)
        text3d(pca$x, pca$y, pca$z, group.v, cex=0.6, adj=2)
        par3d(mouseMode = "trackball")
        rglwidget()
      })
    }
  )
  output$heatmap=renderPlotly(
    {
      mat=as.matrix(dist(t(exprs(norm()))))
      rownames(mat)=pData(norm())$SampleID
      colnames(mat)=rownames(mat)
      heatmaply(mat,margins = c(80,120,60,40),colorRampPalette(colors = c("red", "yellow")))
    }
  )
  output$deg=DT::renderDataTable(DT::datatable(
    {
      dat = deg()$mylist[[input$NumContrasts]]
      dat = dat[,-6]
      dat[,6:7] = format(dat[,6:7], scientific = TRUE)
      
      if (is.na(input$pval) & is.na(input$fc)) {   
        dat
      } else if (is.na(input$pval))  {
        dat = dat[(abs(as.numeric(dat[,5])) >= input$fc),]
      } else if (is.na(input$fc)) {
        dat = dat[(as.numeric(dat[,6]) <= input$pval),]
      } else {
        dat = dat[(as.numeric(dat[,6]) <= input$pval & abs(as.numeric(dat[,5])) >= input$fc),]
        dat
      }
      # deg()[[1]]
    }, caption =paste0("contrast: ",names(deg()$mylist)[input$NumContrasts])
  )
  )
  })
})
})


