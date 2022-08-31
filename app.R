
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
      # Read data ----
      tabItem(tabName = "readData",
              h1("DIA-NN Data"),
              fileInput("dataFile",label = NULL,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"),
              
              
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
                
                
                column(4,      # Input: Select separator ----
                       radioButtons(inputId = "sep", 
                                    label = "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = "\t", inline=T)),
                
                # Input: Select quotes ----
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
      
      
      # tables
      tabItem(tabName = "tables",
              h2("Proteins"),
              dataTableOutput('dataTable2'),
              downloadButton("downloadData", "Download"),
              h2("Peptide (input data)"),
              dataTableOutput('dataTable')),
      # Button
      
      
      
      tabItem(tabName = "graphs",
              
              tabBox( width = 0,
                      
                      
                      tabPanel("Boxplot",  
                               textInput("titleboxplot", "Title of Boxplot", "Title"),
                               textInput("xboxplot", "X label", "Conditions"),
                               textInput("yboxplot", "Y label", "Log2Intensity"),
                               plotlyOutput('boxplot',width = "100%", height = "800"),
                               h2("DATA"), 
                               dataTableOutput('dataTable3'),
                               downloadButton("downloadDataBoxplot", "Download")),
                      tabPanel("Barplot",
                               textInput("titlebarplot", "Title of Barplot", "Title"),
                               textInput("xbarplot", "X label", "Conditions"),
                               textInput("ybarplot", "Y label", "Log2Intensity"),
                               plotOutput('barplot',width = "100%", height = "800"),
                               h2("DATA"),
                               dataTableOutput('dataTable4')),
                      tabPanel("Multiscatterplot",
                               helpText("This may take up to two minutes...."),
                               plotOutput('scatterplot',width = "100%", height = "800")),
                      tabPanel("UpsetR",
                               helpText("This plot is not possible if you have only one sample"),
                               uiOutput('orderselect'),
                               uiOutput('creasing'),
                               plotOutput('upsetR',width = "100%", height = "800")),
                      
                      tabPanel("Distribution CV",
                               helpText("This plot is not possible if you have no technical replicat"),
                               plotOutput('piratplot',width = "100%", height = "800"))
                      
                      
                      #downloadButton("downloadData", "Download")
                      
              )              
              
      )
    )
  )
)

################################################################################
# Server
################################################################################

server <- function(input, output, session) {
  
  data <- reactiveValues()
  
  #=============================================================================
  # Preview
  #=============================================================================
  output$preview <-  renderDataTable({
    
    req(input$dataFile)
    
    df <- read.csv(input$dataFile$datapath,
                   header = as.logical(input$header),
                   sep = input$sep,
                   quote = input$quote,
                   nrows=5
                   
                   
    )
  })
  
  
  
  #=============================================================================
  # Lecture
  #=============================================================================
  observeEvent(input$actBtnVisualisation, {
    
    if(!is.null(input$dataFile$datapath)){
      data$table = read.csv(input$dataFile$datapath,
                            header = as.logical(input$header),
                            sep = input$sep,
                            quote = input$quote)
      sendSweetAlert(
        session = session,
        title = "Done !",
        text = "Le fichier a bien été lu !",
        type = "success"
      )
      
      updateTabItems(session, "tabs", selected = "tables")
    }
    
  })
  
  #=============================================================================
  # Exploration du tableau
  #=============================================================================
  
  output$dataTable <- DT::renderDataTable({
    DT::datatable(data$table, options = list(orderClasses = TRUE,scrollX = TRUE))
  })
  
  #=============================================================================
  # Tableau format
  #=============================================================================
  
  
  
  
  fileformat <- reactive({
    formatDIANN(data$table)
    
  })
  
  fileformat2 <- reactive({
    Barformat(fileformat())
    
    
  })
  
  fileformat3 <- reactive({
    data_summary(fileformat2(),varname="nbProt",groupnames = "groupe")
    
    
  })
  
  fileformat4 <- reactive({
    formatUpset(fileformat())
    
    
  })
  
  fileformat5 <- reactive({
    formatpirateplot(fileformat())
    
    
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
  
  
  output$dataTable2 <- DT::renderDataTable({
    
    DT::datatable(fileformat(), options = list(orderClasses = TRUE,scrollX = TRUE))
  })
  
  output$dataTable3 <- DT::renderDataTable({
    
    DT::datatable(fileformat2())
  }) 
  
  output$dataTable4 <- DT::renderDataTable({
    
    DT::datatable(fileformat3())
  }) 
  
  output$dataTableUpset <- DT::renderDataTable({
    
    DT::datatable(fileformat4())
  }) 
  
  output$dataPiratplot <- DT::renderDataTable({
    
    DT::datatable(fileformat5())
  }) 
  
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
  #output$accessionselect <- renderUI({
  
  # accession <- fileformat()$Protein.Ids
  #selectizeInput(inputId="AccessionProtein",
  #               label="Accession Choice",
  #               choices = accession,
  #               width = '500px',
  #               selected = "1"
  
  #)
  
  
  # })
  
  #  output$accessionselect2 <- renderUI({
  
  #    accession <- fileformat()$Protein.Ids
  #   selectizeInput(inputId="AccessionProtein",
  #                   label="Accession Choice",
  #                   choices = accession,
  #                  width = '500px',
  #                 selected = "54"
  
  # )
  
  
  #})
  
  
  
  
  
  
  
  
  
  
  
  
  #=============================================================================
  # Graphiques
  #=============================================================================
  
  
  
  
  #Multiscatterplot
  
  output$scatterplot <- renderPlot({
    
    df <- ggpairsformat(fileformat())    
    ggpairs(df) 
    
    
  })
  
 
  #### Boxplot ----
  
  plot_output <- reactive({
    df2 <- boxplotformat(fileformat())
    p<-ggplot(df2,aes(x=name,y=value,col=name))+geom_boxplot()+ggtitle(req(input$titleboxplot))+labs(x=input$xboxplot,y=input$yboxplot,colour=input$xboxplot)
    q <-ggplotly(p)
    q
    
  })
  
  
  
   output$boxplot <- renderPlotly(plot_output())
    
   
    
     #### save file svg ----
     output$DownloadBoxPlot <- downloadHandler(
       filename = function() {
         paste("boxplot", input$extension, sep = ".")
       },
       content = function(file){
        
          #ggsave(file, plot_output(), device = input$extension)
        plotly_IMAGE(plot_output(),format="png",out_file="test.svg")
       })
     
   output$downloadDataBoxplot <- downloadHandler(
     filename = function() {
       paste("dataBoxPlot.csv", sep = "")
     },
     content = function(file) {
       write.csv(fileformat2(), file, row.names = FALSE)
     })
  
   
   
   
  output$barplot <- renderPlot({
    
    df3 <- fileformat3()
    ggplot(df3,aes(x=groupe,y=nbProt))+ geom_col(fill="lightgray",color="black")+ geom_errorbar(aes(ymin=nbProt-sd,ymax=nbProt+sd),width=0.5)+coord_flip()
  })
  
  
  
  
  #outfile <- tempfile(fileext='.png')
  output$upsetR <- renderPlot({
    
    df4 <- fileformat4()
    
    p <- UpSetR::upset(data = df4,
                       order.by = orderdat(),
                       decreasing = c(decrease()))
    p
    
  })
  
  output$piratplot <- renderPlot({
    
    df5 <- fileformat5()
    P<- pirateplot(formula=cv~condition,data=df5,theme=0,pal = "southpark",bean.f.o = .6, point.o = .3,inf.b.o = .8,avg.line.o=1, avg.line.col = "black")
    p
  })
  
  
  
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("report_formated.csv", sep = "")
    },
    content = function(file) {
      write.csv(fileformat(), file, row.names = FALSE)
    })
  
 
  
  
  
}



#p<-profvis::profvis(runApp(shinyApp(ui, server)))
#htmlwidgets::saveWidget(p, "profile.html")
# Can open in browser from R
#browseURL("profile.html")

shinyApp(ui, server)




