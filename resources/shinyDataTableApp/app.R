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
      checkboxInput("thousandgenome", "ThousandGenomes SV", value=TRUE)
      
    ),
    mainPanel( 
      dataTableOutput("datatable") 
    )
  )
)


server <- function(input, output) {
  
  filteredsvs <- reactive( {
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
    #if (str_extract(mypos, "\\S+:{0,1}\\d+\\-\\d+") !="") {
      #start <- sub('\\S+:([0-9]+)\\-([0-9]+)', '\\1', input$poslowhigh)
      #end <- sub('\\S+:([0-9]+)\\-([0-9]+)', '\\2', input$poslowhigh)
      #outputvars <- outputvars[(!(outputvars$end < start)) & (!(outputvars$start > end)),]
    #}
    if (input$poslowhigh != "") {
      #m <- regexpr("(\d+)\-(\d+)", x, perl=TRUE)
      #c(low, high) <- 
      #outputvars <- outputvars[(!(outputvars$end < as.integer(input$poslow)) & !(outputvars$start > as.integer(input$poshigh))), ]
    }
    if (!input$gnomad) {
      outputvars <- outputvars[outputvars$gnomad_SV_counts != 1,]
    }
    if (!input$thousandgenome) {
      outputvars <- outputvars[outputvars$X1kg_SV_counts != 1,]
    }
    outputvars <- outputvars[outputvars$SV_type %in% input$typeInput, ]
    
    outputvars
  })
  
  output$datatable <- renderDataTable(filteredsvs())
}

shinyApp(ui, server)

