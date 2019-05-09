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

#setwd("/Users/nhansen/WomenHackathon/shinySVapp")
sv <- read_delim("/Users/ariel/Dropbox/JHU/hackathon/Hiding-in-plain-sight-unannotated-structural-variants-in-public-genomic-data-sets/resources/shiny_graphics/test.bed", delim = "\t", col_names = TRUE) %>%
  rename("chr" = "#chr") %>%
  mutate(mean_bitscore=as.integer(mean_bitscore)) %>%
  mutate(width = (end-start))

deletion.track <- AnnotationTrack(data=sv,
                                  start=sv$start, 
                                  width = sv$width, 
                                  genome="hg19", 
                                  type="l", 
                                  name="SV", 
                                  window=10,
                                  #group = rep(grouplist),
                                  # strand = "+", "-",
                                  chromosome="chr21")
# build deletion track 

ui <- fluidPage (
  titlePanel("Hiding in Plain Sight"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId="chr", 
                  label="Chromosome",
                  choices=unique(sv$chr)),
      textInput(inputId="start", label="Minimum position (e.g., 10000-20000)", value = "14961234", width = NULL,
                placeholder = NULL),
      textInput(inputId="end", label="Maximum position (e.g., 10000-20000)", value = "16065902", width = NULL,placeholder = NULL),
      sliderInput("scoreMinMax", "Score", 30, 200, c(60, 100)),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=sv$SV_type[1]),
      checkboxInput("gnomad", "Gnomad SV", value=TRUE),
      checkboxInput("thousandgenome", "ThousandGenomes SV", value=TRUE)
      
    ),
    mainPanel(
      imageOutput("myImage")
    )
  )
)

 

server <- function(input, output) {
  sv <- sv
  filteredsvs <- reactive({
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

 
}

shinyApp(ui, server)
  
 