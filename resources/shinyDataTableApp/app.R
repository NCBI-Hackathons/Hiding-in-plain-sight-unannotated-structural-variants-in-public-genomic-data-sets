library("shiny")
library("shinyWidgets")
library("stringr")


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
                         selected=sv$SV_type[1]),
      checkboxGroupInput("sv_seen", "Annotation",
                         choices=names(sv)[8:ncol(sv)],
                         selected=names(sv)[8:ncol(sv)]),
      radioButtons("display", "", choiceNames=c("Table", "Figure"), choiceValues=c("Table", "Figure"), selected="Figure", inline=TRUE),
      downloadButton("downloadData", "Download")
      
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
      if (!is.na(str_extract(input$poslowhigh, "(\\S*:){0,1}\\d+\\-\\d+"))) {
        start <- as.integer(sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\2', input$poslowhigh))
        end <- as.integer(sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\3', input$poslowhigh))
        outputvars <- outputvars[!((outputvars$end < start) | (outputvars$start > end)),]
      }
      gnomad = "gnomad_SV_counts" %in% input$sv_seen
      if (!gnomad) {
        outputvars <- outputvars[outputvars$gnomad_SV_counts == 0,]
      }
      kg = "X1kg_SV_counts" %in% input$sv_seen
      if (!kg) {
        outputvars <- outputvars[outputvars$X1kg_SV_counts == 0,]
      }
      simpleseq = "simple_sequence_flag" %in% input$sv_seen
      if (!simpleseq) {
        outputvars <- outputvars[outputvars$simple_sequence_flag == 0,]
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
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      "svdataset.csv"
    },
    content = function(file) {
      write.csv(filteredsvs(), file, row.names = FALSE)
    }
  )

  output$testtext <- renderText(descriptivetext()) 
  output$datatable <- renderDataTable(filteredsvs())

}

shinyApp(ui, server)

