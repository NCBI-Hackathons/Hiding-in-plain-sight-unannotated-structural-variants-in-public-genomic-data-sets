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

setwd("/Users/rtorene9887/projects/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/resources/MergedApps/")
sv <- read_delim("/Users/rtorene9887/projects/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/output/annotated_svs.bed", delim = "\t", col_names = TRUE)[1:50000,] %>%
  rename("chr" = "#chr") %>%
  mutate(mean_bitscore=as.integer(mean_bitscore)) %>%
  mutate(width = (end-start))

min_bitscore<-min(sv$mean_bitscore, na.rm=T)
max_bitscore<-max(sv$mean_bitscore, na.rm=T)

sv$mean_bitscore <- as.integer(sv$mean_bitscore)
sv$start <- as.integer(sv$start)
sv$end <- as.integer(sv$end)


deletion.track <- AnnotationTrack(data=sv,
                                  start=sv$start, 
                                  width = sv$width, 
                                  genome="hg19", 
                                  type="l", 
                                  name="SV", 
                                  window=10,
                                  #group = rep(grouplist),
                                  # strand = "+", "-",
                                  chromosome=sv$chr)
# build deletion track 
annotation_options = names(sv)[8:ncol(sv)]
annotation_options = annotation_options[-match("width",annotation_options)]
# annotation_options = annotation_options[-match("pseudogenes",annotation_options)]

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
      sliderInput("scoreMinMax", "Score", min_bitscore, max_bitscore, c(60, max_bitscore)),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=sv$SV_type[1]),
      checkboxGroupInput("sv_seen", "Annotation",
                         choices=annotation_options,
                         selected=annotation_options),
      radioButtons("display", "", choiceNames=c("Table", "Figure"), choiceValues=c("Table", "Figure"), selected="Table", inline=TRUE),
      downloadButton("downloadData", "Download")
      
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
  

    filtered <- reactive( {
         outputvars <- sv
         if (is.integer(input$scoreMinMax[1])) {
           outputvars <- outputvars[outputvars$mean_bitscore >= input$scoreMinMax[1], ]
         }
         if (is.integer(input$scoreMinMax[2])) {
           outputvars <- outputvars[outputvars$mean_bitscore <= input$scoreMinMax[2], ]
         }
        # if (input$chr != "") {
         #  outputvars <- outputvars[outputvars$chr==input$chr,]
         #}
         #if (!is.na(str_extract(input$poslowhigh, "(\\S*:){0,1}\\d+\\-\\d+"))) {
          # outputvars <- outputvars[c(),]
           #start <- sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\1', input$poslowhigh)
           #end <- sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\2', input$poslowhigh)
           #outputvars <- outputvars[(!(outputvars$end < start) & !(outputvars$start > end)),]
       #  if (input$start != "") {
      #     outputvars <- outputvars[outputvars$start==input$start,]
      #   }
      #   if (input$end != "") {
      #     outputvars <- outputvars[outputvars$end==input$end,]
      #   }
           
           
         
         # gnomad = "gnomad_SV_counts" %in% input$sv_seen
         # if (!gnomad) {
         #   outputvars <- outputvars[outputvars$gnomad_SV_counts == 0,]
         # }
         # kg = "X1kg_SV_counts" %in% input$sv_seen
         # if (!kg) {
         #   outputvars <- outputvars[outputvars$X1kg_SV_counts == 0,]
         # }
         # simpleseq = "simple_sequence_flag" %in% input$sv_seen
         # if (!simpleseq) {
         #   outputvars <- outputvars[outputvars$simple_sequence_flag == 0,]
         # }
        outputvars <- outputvars[outputvars$SV_type %in% input$typeInput, ]
        
        
        outputvars=sv %>%
          filter(chr == as.character(input$chr)) %>%
          filter(start >= as.numeric(input$start)) %>%
          filter(end <= as.numeric(input$end))
        
        outputvars
      
  })
  
  filteredsvs1 <- reactive({
   
    sv <- sv
    
    outputvars=sv %>%
      #filter(SV_type == "Insertion") %>%
      filter(chr == as.character(input$chr)) %>%
      filter(SV_type%in%input$typeInput) %>%
      filter(start >= as.numeric(input$start)) %>%
      filter(end <= as.numeric(input$end)) %>%
      mutate(width = (end-start)+1) 
    #   filter(mean_bitscore >= as.numeric(input$scoreMinMax[1])) %>%
    #  filter(mean_bitscore <= as.numeric(input$input$scoreMinMax[2]))
    
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
    return(deletion.track)
    
  
  })
  output$results <- renderTable({
    filtered()
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
      plotTracks(list(ideoTrack,axisTrack,biomTrack,deletion.track), from = start, to = stop, labelPos = "below")
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
  
  
  
  
  
  output$datatable <- renderDataTable(filtered())
 
 
    
  
  
}

shinyApp(ui, server)

