# clear workspace 
rm(list=ls())

# load libraries NRIP1
library("ggplot2")
library("Gviz")
library("GenomicRanges")
library("tidyverse")
library("shiny")
library("shinyWidgets")
library("stringr")
# 28477797-33448354

sv <- read_delim("annotated_svs.bed", delim = "\t", col_names = TRUE) %>%
  rename("chr" = "#chr") %>%
  mutate(mean_bitscore=as.integer(mean_bitscore)) %>%
  mutate(width = (end-start))

annotation_options = names(sv)[8:ncol(sv)]
annotation_options = annotation_options[-match("width",annotation_options)]

min_bitscore<-min(sv$mean_bitscore)
max_bitscore<-max(sv$mean_bitscore)

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
      textInput(inputId="start", label="Start position (e.g., 14358791)", value = "14358791", width = NULL,
                placeholder = NULL),
      textInput(inputId="end", label="Stop position (e.g., 14366188)", value = "14366188", width = NULL,placeholder = NULL),
      sliderInput("scoreMinMax", "Score", min_bitscore, max_bitscore, c(60, 100)),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=sv$SV_type[1]),
      checkboxGroupInput("annotations", "Annotations",
                         choices=annotation_options,
                         selected=c())
    ),
    mainPanel(
      tabsetPanel(
      
      tabPanel("Table",dataTableOutput("datatable")),
      
      tabPanel("Image", imageOutput("myImage"))
      )
    )
  )
)



server <- function(input, output) {

  filtered_svs <- reactive({
    outputvars = sv %>%
      filter(chr == as.character(input$chr)) %>%
      filter(SV_type%in%input$typeInput) %>%
      filter(start >= as.numeric(input$start)) %>%
      filter(end <= as.numeric(input$end)) %>%
      filter(mean_bitscore >= as.numeric(input$scoreMinMax[1])) %>%
      filter(mean_bitscore <= as.numeric(input$scoreMinMax[2]))

    for (annottype in annotation_options) {
      if (annottype %in% input$annotations) {
        outputvars <- outputvars[ outputvars[,annottype] > 0, ]
      }
    }
    return(outputvars)
  })  
  
  
  deletiontrack <- reactive({
    outputvars <- filtered_svs()
    if (length(outputvars$start) > 0) {
      deletion.track <- AnnotationTrack(data=outputvars,
                                        start=outputvars$start, 
                                        width = outputvars$width, 
                                        genome="hg19", 
                                        type="l", 
                                        name="SV", 
                                        window=10,
                                        #group = rep(grouplist),
                                        # strand = "+", "-",
                                        chromosome=input$chr)
      feature(deletion.track) <- outputvars$SV_type
      return(deletion.track)
    }
    else {
      return(NULL)
    }
    
  })
  
  output$myImage <- renderImage({
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.pdf')
    
    func_chr <- function(chr, start, stop) {
      
      # build chromosome model 
      ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
      # build gene models 
      axisTrack <- GenomeAxisTrack()
      data(geneModels)
      biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = stop, name = "Gene Model", transcriptAnnotation = "symbol")
      # plot tracks
      mydeletiontrack <- deletiontrack()
      if (is.null(mydeletiontrack)) {
        plotTracks(list(ideoTrack,axisTrack,biomTrack), from = start, to = stop, labelPos = "below")
      }
      else {
        plotTracks(list(ideoTrack,axisTrack,biomTrack,mydeletiontrack), from = start, to = stop, labelPos = "below", Deletion="darkred", Insertion="darkgreen", Inversion="purple", TandemDup = "blue", MEI = "yellow", legend=TRUE)
      }
    }
    # generate png
    png(outfile, width = 1200, height = 1000, res = 210)
    func_chr(input$chr, as.numeric(input$start), as.numeric(input$end))
    dev.off()
    
    list(src = outfile,
         contentType = 'image/pdf',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = TRUE)

  output$datatable <- renderDataTable(filtered_svs())
}

shinyApp(ui, server)

