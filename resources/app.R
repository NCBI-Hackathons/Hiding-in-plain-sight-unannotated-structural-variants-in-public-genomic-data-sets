library("shiny")
library("shinyWidgets")
library("stringr")


setwd("/Users/rtorene9887/projects/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/resources")

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
      sliderInput("scoreMinMax", "Score", min(sv$mean_bitscore, na.rm=T), max(sv$mean_bitscore, na.rm=T), c(50, max(sv$mean_bitscore, na.rm=T))),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=unique(sv$SV_type)),
      checkboxGroupInput("sv_seen", "Annotation",
        choices=names(sv)[8:ncol(sv)],
        selected=names(sv)[8:ncol(sv)] )
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
    if (input$poslowhigh != "") {
      #m <- regexpr("(\d+)\-(\d+)", x, perl=TRUE)
      #c(low, high) <- 
      #outputvars <- outputvars[(!(outputvars$end < as.integer(input$poslow)) & !(outputvars$start > as.integer(input$poshigh))), ]
    }
    gnomad = "gnomAD" %in% input$sv_seen
    if (!gnomad) {
      outputvars <- outputvars[outputvars$gnomad_SV_counts == 0,]
    }
    kg = "1000_Genomes" %in% input$sv_seen
    if (!kg) {
      outputvars <- outputvars[outputvars$X1kg_SV_counts == 0,]
    }
    
    outputvars
  })
  
  outputvars <- sv
  output$displaytext <- renderText({paste("THIS:", lowscore) })
  
  output$selected_var <- renderText({ 
    input$typeInput
  })
  output$datatable <- renderDataTable(filteredsvs())
}

shinyApp(ui, server)

