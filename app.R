
### Library  and functions ----


library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(ggplot2)
library(googleVis)
library(colourpicker)
library(reshape2)
library(dplyr)
library(tidyr)
library(Hmisc)
library(profvis)
library(data.table)
library(GGally)
library(plyr)
library(UpSetR)
library(yarrr)
library(reticulate)
library(stringr)
library(protr)


#f_source(f_formatDIANN.R)
options(shiny.maxRequestSize=100*1024^5)
options(encoding="UTF-8")

source("f_formatDIANN.R", encoding = 'UTF-8')





### UI ----



ui <- dashboardPage(
  dashboardHeader(title = "Process Data DIANN"),
  dashboardSidebar( width = 120,
                    sidebarMenu(id = "tabs",
                                menuItem("Read data", tabName = "readData", icon = icon("readme")),
                                menuItem("Tables", tabName = "tables", icon = icon("table")),
                                menuItem("Graphs", tabName = "graphs", icon = icon("poll"))
                    )
  ),
  dashboardBody(
    
    tabItems(
      #### Read data ----
      tabItem(tabName = "readData",
              h1("DIA-NN Data"),
              fileInput("dataFile",label = NULL,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"),
              h2("Fasta file (if you want %coverage of protein)"),
              fileInput("fastaFile",label = NULL,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"
                        ),
              
              
              fluidRow(
                column(12,
                       h3("Parameters"), helpText("To check if your file is readable")),
                
                # Input: Checkbox if file has header
                column(4,
                       radioButtons(inputId = "header", 
                                    label = "Header",
                                    choices = c("Yes" = TRUE,
                                                "No" = FALSE),
                                    selected = TRUE, inline=T)),
                
                
                column(4,      # Input: Select separator 
                       radioButtons(inputId = "sep", 
                                    label = "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = "\t", inline=T)),
                
                # Input: Select quotes 
                column(4,
                       radioButtons(inputId = "quote", 
                                    label= "Quote",
                                    choices = c(None = "",
                                                "Double Quote" = '"',
                                                "Single Quote" = "'"),
                                    selected = "", inline=T))
                ,
                column(8,
                       h3("File preview"),
                       dataTableOutput(outputId = "preview")
                )
              ), 
              tags$br(),
              
              div(actionButton(inputId = "actBtnVisualisation", label = "Integration",icon = icon("play") ), align = "center")
              
              
              
      ),
      
      
      #### Display tables ----
      tabItem(tabName = "tables",
              h2("Proteins"),
              dataTableOutput('dataTableProteins'),
              downloadButton("downloadData", "Download"),
              h2("Peptide (input data)"),
              dataTableOutput('dataTablePeptides')),
              
      
      
      #### Display plot ----
      tabItem(tabName = "graphs",
              
              tabBox( width = 0,
                      
                      ##### Boxplot ----
                      tabPanel("Boxplot",  
                               textInput("titleboxplot", "Title of Boxplot", "Title"),
                               textInput("xboxplot", "X label", "Conditions"),
                               textInput("yboxplot", "Y label", "Log2Intensity"),
                               plotlyOutput('boxplot',width = "100%", height = "800"),
                               h2("DATA"), 
                               #dataTableOutput('dataTableBoxplot'),
                               downloadButton("downloadDataBoxplot", "Download")),
                      
                      ##### Barplot ----
                      tabPanel("Barplot",
                               textInput("titlebarplot", "Title of Barplot", "Title"),
                               textInput("xbarplot", "X label", "Conditions"),
                               textInput("ybarplot", "Y label", "nbProt"),
                               plotlyOutput('barplot',width = "100%", height = "800"),
                               h2("DATA"),
                               #dataTableOutput('dataTable4')),
                               downloadButton("downloadDataBarplot", "Download")),
                      
                      ##### Multiscatterplot ----
                      tabPanel("Multiscatterplot",
                               helpText("This may take up to two minutes...."),
                               plotOutput('scatterplot',width = "100%", height = "800")),
                      
                      ##### UpsetR Plot ----
                      tabPanel("UpsetR",
                               helpText("This plot is not possible if you have only one sample"),
                               uiOutput('orderselect'),
                               uiOutput('creasing'),
                               plotOutput('upsetR',width = "100%", height = "800"),
                               downloadButton("downloadDataUpset", "Download"),
                               helpText("This plot is not possible if you have only one sample")),
                      
                      ##### Pirateplot (violin plot) ----
                      tabPanel("Distribution CV",
                               helpText("This plot is not possible if you have no technical replicat"),
                               plotlyOutput('piratplot',width = "100%", height = "800"))
                      
                      
                     
                      
              )              
              
      )
    )
  )
)


#### Server ----


server <- function(input, output, session) {
  
  data <- reactiveValues()
  data2 <- reactiveValues()
  
  ## Preview ----
  
  output$preview <-  renderDataTable({
    
    req(input$dataFile)
    
    df <- read.csv(input$dataFile$datapath,
                   header = as.logical(input$header),
                   sep = input$sep,
                   quote = input$quote,
                   nrows=5
                   
                   
    )
  })
  
  
  
  
  ## Integration ----
  
  observeEvent(input$actBtnVisualisation, {
    
    if(!is.null(input$dataFile$datapath)){
      data$table = read.csv(input$dataFile$datapath,
                            header = as.logical(input$header),
                            sep = input$sep,
                            quote = input$quote)
      nbline <- nrow(data$table)
      nbcol <- ncol(data$table)
      confirmation <- paste("Your file has",nbline,"lines and",nbcol,"columns")
      sendSweetAlert(
        session = session,
        title = "Done !",
        text = confirmation,
        type = "success"
      )
      
     
      
      
    }
    if(!is.null(input$fastaFile$datapath)){
      data2$table = readFASTA(input$fastaFile$datapath)
      cat("test si ça marche")
      nbseq = length(data2$table)
      confirmation <- paste("Your fasta file has",nbseq,"sequence(s)")
      sendSweetAlert(
        session = session,
        title = "Done !",
        text = confirmation,
        type = "success"
      )  
    } 
    updateTabItems(session, "tabs", selected = "tables")
  })
  
 # observeEvent(input$test, {
    
  #  if(!is.null(input$fastaFile$datapath)){
  #    data2$table = readFASTA(input$fastaFile$datapath)
   #   nbseq = length(data2$table)
  #    confirmation <- paste("Your fasta file has",nbseq,"sequence(s)")
   #  sendSweetAlert(
  #     session = session,
  #     title = "Done !",
  #     text = confirmation,
  #     type = "success"
  #   )
      
  #   updateTabItems(session, "tabs", selected = "tables")
  #   }
    
#}) 
  
  ## Reading and exploration ----
  
  output$dataTablePeptides <- DT::renderDataTable({
    DT::datatable(data$table, options = list(orderClasses = TRUE,scrollX = TRUE))
  })
  
  ## Format data for datatables and plots ----
  
  ### Datatable "Proteins" -----
  
  
  cover<- reactive({
    
    protcov(data$table,data2$table)
  })
  
  output$dataTablecover <- DT::renderDataTable({
    DT::datatable(cover())
  })
 
  
  fileformat <- reactive({
    formatDIANN(data$table,cover())
    
  })
  
 
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("report_formated.csv", sep = "")
    },
    content = function(file) {
      write.csv(fileformat(), file, row.names = FALSE)
    })
  
  ### Boxplot data ----- 
  fileformat2 <- reactive({
    
   boxplotformat(fileformat()) 
    
  })
  
  ### Barplot data ----- 
    fileformat3a <- reactive({
    Barformat(fileformat())
    
    
  })
  
    fileformat3 <- reactive({
    data_summary(fileformat3a(),varname="nbProt",groupnames = "groupe")
    })  
   
  
    ### UpsetR data ----- 
    fileformat4 <- reactive({
    formatUpset(fileformat())
    
    
  })
  
  dataintersection <- reactive({
    
    formatListUpsetData(fileformat4())
  })
    
    ### Violin data ----- 
    fileformat5 <- reactive({
    formatpirateplot(fileformat())
    
    
  })

  
 
   
  ## Datatable for each plot ----
  
  output$dataTableProteins <- DT::renderDataTable({
    
    DT::datatable(fileformat(), options = list(orderClasses = TRUE,scrollX = TRUE))
  })
  
  output$dataTableBoxplot <- DT::renderDataTable({
    
    DT::datatable(fileformat2())
  }) 
  
  output$dataTableBarplot <- DT::renderDataTable({
    
    DT::datatable(fileformat3())
  }) 
  
  output$dataTableUpset <- DT::renderDataTable({
    
    DT::datatable(fileformat4())
  }) 
  
  output$dataPiratplot <- DT::renderDataTable({
    
    DT::datatable(fileformat5())
  }) 
  
  
  
  
  # Plots ----
  
  #### Boxplot ----
  
  boxplot_output <- reactive({
    df2 <- fileformat2()
    p<-ggplot(fileformat2(),aes(x=name,y=value,col=name))+geom_boxplot()+ggtitle(req(input$titleboxplot))+labs(x=input$xboxplot,y=input$yboxplot,colour=input$xboxplot)
    q <-ggplotly(p)
    q
    
  })
 
  output$boxplot <- renderPlotly(boxplot_output())
  
  
  ### to download data of boxplot
  
  
  output$downloadDataBoxplot <- downloadHandler(
    filename = function() {
      paste("dataBoxPlot.csv", sep = "")
    },
    content = function(file) {
      write.csv(fileformat2(), file, row.names = FALSE)
    })
  
  #### Barplot ----
  
  barplot_output <- reactive({
    
    df3 <- fileformat3()
    ngroup <- nrow(df3)
    bp <-ggplot(df3,aes(x=groupe,y=nbProt,fill=groupe))+geom_bar(stat="identity")+scale_fill_brewer()+ geom_errorbar(aes(ymin=nbProt-sd,ymax=nbProt+sd),width=0.5)+ggtitle(req(input$titlebarplot))+labs(x=input$xbarplot,y=input$ybarplot)
    bar <-ggplotly(bp)
    bar
  
  }) 
  
  output$barplot <-  renderPlotly(barplot_output())
  
  ###Multiscatterplot ----
  
  output$scatterplot <- renderPlot({
    
    df <- ggpairsformat(fileformat())    
    ggpairs(df) 
    
    
  })
  
  
  ###UpsetR----
  
  
  output$upsetR <- renderPlot({
    
    df4 <- fileformat4()
    
    p <- UpSetR::upset(data = df4,
                       text.scale = 3,
                       order.by = orderdat(),
                       decreasing = c(decrease()))
    p
    
  })
  
  
  
  output$orderselect <- renderUI({
    
    
    selectizeInput(inputId="order",
                   label="Order by",
                   choices =list("Degree" = "degree","Frequency" = "freq"),
                   width = '250px',
                   selected = "1"
                   
    )})
  
  output$creasing <- renderUI({
    
    
    selectizeInput(inputId="crease",
                   label="Increasing/Decreasing",
                   choices =list("Increasing" = "inc","Decreasing" = "dec"),
                   width = '250px',
                   selected = "1"
                   
    )})  
  
  
  
  orderdat <- reactive({
    orderdat <- as.character(input$order)
    if(orderdat == "degree"){
      orderdat <- c("degree")
    }
    else if(orderdat == "freq"){
      orderdat <- "freq"
    }
    return(orderdat)
  })
  
  decrease <- reactive({
    decrease <- as.character(input$crease)
    if(decrease == "inc"){
      decrease <- FALSE
    }
    else if(decrease == "dec"){
      decrease <- TRUE
    }
    return(decrease)
  })
  
  
  
  output$downloadDataUpset <- downloadHandler(
    filename = function() {
      paste("dataUpsetIntersection.csv", sep = "")
    },
    content = function(file) {
      write.csv(dataintersection(), file, row.names = FALSE)
    })
  
  ###Pirateplot(violin plot)----
  
   
    
  piratplot_output <- reactive({
    
    df5 <- fileformat5()
    
    jack<- ggplot(df5,aes(x=condition,y=cv,fill=condition))+geom_violin()+geom_point( color="black",size=0.1, position = position_jitter(w=0.01))+geom_crossbar(stat="summary", fun.y=mean, fun.ymax=mean, fun.ymin=mean,fatten=2, width=.5)+scale_fill_brewer(palette="Spectral")+theme_bw()
    sparrow <-ggplotly(jack)
    sparrow
  })
  
  output$piratplot <-  renderPlotly(piratplot_output())
  
  
  
 
  
 
  
  
  
}





shinyApp(ui, server)




