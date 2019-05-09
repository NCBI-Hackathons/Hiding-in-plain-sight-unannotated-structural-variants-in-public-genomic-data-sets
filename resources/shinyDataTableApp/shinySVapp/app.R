library("shiny")
library("shinyWidgets")
library("stringr")

setwd("/Users/nhansen/WomenHackathon/shinySVapp")

sv <- read.table("chr21_deletions_annotated.bed", sep="\t", header=TRUE)
sv$mean_bitscore <- as.integer(sv$mean_bitscore)
sv$start <- as.integer(sv$start)
sv$end <- as.integer(sv$end)

ui <- fluidPage (
  titlePanel("Hiding in Plain Sight"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId="chr", 
                  label="Chromosome",
                  choices=unique(sv$chr)),
      textInput(inputId="poslowhigh", label="Minimum/maximum position (e.g., 10000-20000)", value = "", width = NULL,
                placeholder = NULL),
      sliderInput("scoreMinMax", "Score", 30, 300, c(60, 100)),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=sv$SV_type[1]),
      checkboxInput("gnomad", "Gnomad SV", value=TRUE),
      checkboxInput("thousandgenome", "ThousandGenomes SV", value=TRUE),
      radioButtons("display", "", choiceNames=c("Table", "Figure", "Download"), choiceValues=c("Table", "Figure", "Download"), selected="Figure", inline=TRUE)
      
    ),
    mainPanel( 
      dataTableOutput("datatable"),
      textOutput("testtext")
    )
  )
)


server <- function(input, output) {
  
  filteredsvs <- reactive( {
    if (input$display=="Table") {
      outputvars <- sv
      if (is.integer(input$scoreMinMax[1])) {
        outputvars <- outputvars[outputvars$mean_bitscore >= input$scoreMinMax[1], ]
      }
      if (is.integer(input$scoreMinMax[2])) {
        outputvars <- outputvars[outputvars$mean_bitscore <= input$scoreMinMax[2], ]
      }
      if (input$chr != "") {
        outputvars <- outputvars[outputvars$chr==input$chr,]
      }
      if (str_extract(mypos, "\\S*:{0,1}\\d+\\-\\d+") !="") {
        #start <- sub('(\\S+:){0,1}([0-9]+)\\-([0-9]+)', '\\1', input$poslowhigh)
        #end <- sub('(\\S+:){0,1}([0-9]+)\\-([0-9]+)', '\\2', input$poslowhigh)
        #outputvars <- outputvars[(!(outputvars$end < start)) & (!(outputvars$start > end)),]
      }
      if (!input$gnomad) {
        outputvars <- outputvars[outputvars$gnomad_SV_counts != 1,]
      }
      if (!input$thousandgenome) {
        outputvars <- outputvars[outputvars$X1kg_SV_counts != 1,]
      }
      outputvars <- outputvars[outputvars$SV_type %in% input$typeInput, ]
      
      outputvars
    }
    else if (input$display=="Figure") {
      outputvars <- data.frame()
    }
    
  })
  
  descriptivetext <- reactive({
    if (input$display=="Figure") {
      outputtext <- "This is a lovely figure!"
    }
    else if (input$display=="Table") {
      outputtext <- "This is when we display the table!"
    }
    else if (input$display=="Download") {
      outputtext <- "We will need to download in this case!"
    }
  })
  
  downloadbutton <- reactive({
    if (input$display=="Table") {
      widget <- downloadBttn(download, label = "Download", style = "unite"
                  )
        
        #downloadHandler(
        #filename=function() {"testfile.csv"},
        #content=function(con) {write.csv(filteredsvs(), con)})
      widget
    }
    else {
      widget <- NULL;
    }
    widget
  })

  output$testtext <- renderText(descriptivetext()) 
  output$datatable <- renderDataTable(filteredsvs())
  #output$downloaddata <- downloadbutton()
}

shinyApp(ui, server)

