
################################################################################
# Library
################################################################################

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
library(rsvg)
#f_source(f_formatDIANN.R)
options(shiny.maxRequestSize=100*1024^5)
options(encoding="UTF-8")

source("f_formatDIANN.R", encoding = 'UTF-8')


################################################################################
# UI
################################################################################


ui <- dashboardPage(
  dashboardHeader(title = "Format DIANN"),
  dashboardSidebar( width = 120,
    sidebarMenu(id = "tabs",
                menuItem("Read data", tabName = "readData", icon = icon("readme")),
                menuItem("Tables", tabName = "tables", icon = icon("table")),
                menuItem("Graphs", tabName = "graphs", icon = icon("poll"))
    )
  ),
  dashboardBody(
    
    tabItems(
      # Read data
      tabItem(tabName = "readData",
              h1("Lecture des données"),
              fileInput("dataFile",label = NULL,
                        buttonLabel = "Browse...",
                        placeholder = "No file selected"),
                        
              
              fluidRow(
                column(5,
                       h3("Parameters"),
                       
                       # Input: Checkbox if file has header
                       radioButtons(inputId = "header", 
                                    label = "Header",
                                    choices = c("Yes" = TRUE,
                                                "No" = FALSE),
                                    selected = TRUE, inline=T),
                       
                       # Input: Select separator ----
                       radioButtons(inputId = "sep", 
                                    label = "Separator",
                                    choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"),
                                    selected = "\t", inline=T),
                       
                       # Input: Select quotes ----
                       radioButtons(inputId = "quote", 
                                    label= "Quote",
                                    choices = c(None = "",
                                                "Double Quote" = '"',
                                                "Single Quote" = "'"),
                                    selected = "", inline=T)
                ),
                column(9,
                       h3("File preview"),
                       dataTableOutput(outputId = "preview")
                )
              ), 
              tags$br(),
              
              div(actionButton(inputId = "actBtnVisualisation", label = "Integration",icon = icon("play") ), align = "center")
              
              
              
      ),
      
      # tables
      tabItem(tabName = "tables",
              h2("input exploration"),
              dataTableOutput('dataTable'),
              h2("Proteins"),
              dataTableOutput('dataTable2'),
              # Button
              downloadButton("downloadData", "Download")),
              
              
     tabItem(tabName = "graphs",
                      
                tabBox( width = 0,
                       
                     
                      tabPanel("Boxplot", downloadButton("downloadPlot", "Download"),
                               plotOutput('boxplot',width = "100%", height = "800"),
                               h2("DATA"), 
                               dataTableOutput('dataTable3')),
                      tabPanel("Barplot",
                               plotOutput('barplot',width = "100%", height = "800"),
                               h2("DATA"),
                               dataTableOutput('dataTable4')),
                      tabPanel("Multiscatterplot",
                               helpText("This may take up to two minutes...."),
                               plotOutput('scatterplot',width = "100%", height = "800")),
                      tabPanel("UpsetR",
                               helpText("This plot is not possible if you have only one sample"),
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
                   nrows=10
                   
                   
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
  
  output$boxplot <- renderPlot({
    
    df2 <- boxplotformat(fileformat())    
    boxplot(df2$value~df2$name,ylab = "Abondance",xlab="Conditions",main="All proteins/conditions") 
    
    
  }) 
  
  output$barplot <- renderPlot({
    
    df3 <- fileformat3()
    ggplot(df3,aes(x=groupe,y=nbProt))+ geom_col(fill="lightgray",color="black")+ geom_errorbar(aes(ymin=nbProt-sd,ymax=nbProt+sd),width=0.5)+coord_flip()
  })

  
  #outfile <- tempfile(fileext='.png')
  output$upsetR <- renderPlot({
    
    df4 <- fileformat4()
    
   p <- UpSetR::upset(data = df4)
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


