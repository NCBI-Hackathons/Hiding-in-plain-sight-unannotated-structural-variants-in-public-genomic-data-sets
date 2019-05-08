library(shiny)
library(ggplot2)
library(dplyr)
library(shinyWidgets)

sv <- read.table(file="deletionsBED.txt", sep="\t", header=FALSE,stringsAsFactors = FALSE)
colnames(sv) <- c("chr", "start","stop","ID","score","type")
start_min <- min(data[,"start"])
stop_max <- max(data[,"stop"])

ui <- fluidPage(
  titlePanel("Hiding in plane sight"),
  sidebarLayout(
    sidebarPanel(
      #sliderInput("posInput", "Position", start_min, stop_max, c(5, 50)),
      searchInput(
        inputId = "id", 
        label = "Enter a start position :", 
        placeholder = "20508559", 
        btnSearch = icon("search"), 
        btnReset = icon("remove"), 
        width = "100%"
      ),
      sliderInput("bitInput", "Score", 30, 200, c(60, 100)),
      checkboxGroupInput("typeInput", "SV type",
                   choices = c("Inversions", "deletion"),
                   selected = "Inversions"),
      #radioButtons("gnomADInput", "Variants Type",
      #             choices = c("Only seen in gnomAD", "Only seen here", "All variants"),
      #             selected = "All variants"),
     
      
     
      
     
      uiOutput("ChrOutput"),
      downloadBttn(
        outputId = "downloadData",
        style = "bordered",
        color = "primary"
      )
    ),
   
    
    mainPanel(
   
      tableOutput("results")
    )
  )
)

server <- function(input, output) {
  output$ChrOutput <- renderUI({
    selectInput("chrInput", "Chr",
                sort(unique(sv$chr)),
                selected = "chr6")
  })  
  
 

  filtered <- reactive({
    if (is.null(input$chrInput)) {
      return(NULL)
    }    

    sv %>%
      filter(
            #start >= input$posInput[1],
            #stop <= input$posInput[2],
             score >= input$bitInput[1],
            score <= input$bitInput[2],
             chr == input$chrInput,
            type == input$typeInput,
            start == input$id
            #ID == input$gnomADInput

      )
    
  })
  

  
  output$results <- renderTable({
    filtered()
  })
  
  output$downloadData <- downloadHandler(
    filename = function(file) {
      paste('data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(filtered(), con)
    }
  )
  
}

shinyApp(ui = ui, server = server)