library(shiny)
library(DT)
library(GEOquery)

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
    input$button, 
    isolate({
      
      if (input$gseid=="8 digit GSE code") return()
      
      gds <- getGEO(input$gseid, GSEMatrix = F,getGPL=T,AnnotGPL=T)
      
      mytable=matrix("",length(GSMList(gds)),3)
      colnames(mytable)=c("gsm","title","description")
      for (k in 1:length(GSMList(gds)))
      {
        mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession, Meta(GSMList(gds)[[k]])$title, Meta(GSMList(gds)[[k]])$description)
      }
      
      mytable <- data.frame(mytable)
      
      mytable$group <- ""
      v$data <- mytable
      
      
      output$mytable = DT::renderDataTable({
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8), pageLength = 8))
        
        
      })
      
    })
  )
  
  observeEvent(
    input$button2,
    isolate({
      print( input$group1 )
      print( v$data )
      
      if (is.null(v$data)) return()
      if (is.null(input$group1)) return()
      
      print(v$data)
      print( "Error??2")

      v$data[input$mytable_rows_selected, "group" ] <- input$group1
      print( "Error??3")
      
      output$mytable = DT::renderDataTable({
        
        if (is.null(v$data)) return()
        if (is.null(v$platform)) warning()
        DT::datatable(v$data, options = list(lengthMenu = c(2,4,6,8,10), pageLength = 8))
      })
    })
  )
  
  observeEvent( input$number,
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
  )
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
        
        ), TRUE:input$test)
        data.frame(side1, vs.side2)
      }, sanitize.text.function = function(x) x)
    })
  })
})


#Thanks for everything Fathi! This is my code with the GSE box,
# the displayed table, the group assignment, and the contrast. Hopefully you can connect this with your code 
# and get this application up and running. Onto my poster! :) 