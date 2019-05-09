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
sv <- read_delim("test.bed", delim = "\t", col_names = TRUE) %>%
  rename("chr" = "#chr") %>%
  mutate(mean_bitscore=as.integer(mean_bitscore))
# build deletion track 
# define function
func_chr <- function(chr, start, stop) {
  
  # build chromosome model 
  ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  # build gene models 
  axisTrack <- GenomeAxisTrack()
  data(geneModels)
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = stop, name = "Gene Model", transcriptAnnotation = "symbol")
  # plot tracks
  plotTracks(list(ideoTrack,axisTrack,biomTrack, deletion.track), from = start, to = stop, labelPos = "below")
}

#deletion.track <- AnnotationTrack(data = newsv,
#                                  start=startlist, 
#                                  width = widthlist, 
#                                  genome="hg19", type="l", 
#                                  name="deletion", 
#                                  window=10,
#                                  group = rep(grouplist),
#                                  # strand = "+", "-",
#                                  chromosome=inchr)
#
##sv <- read.table("chr21_deletions_annotated.bed", sep="\t", header#=TRUE)
#sv$mean_bitscore <- as.integer(sv$mean_bitscore)
#sv$start <- as.integer(sv$start)
#sv$end <- as.integer(sv$end)

ui <- fluidPage (
  titlePanel("Hiding in Plain Sight"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId="chr", 
                  label="Chromosome",
                  choices=unique(sv$chr)),
      textInput(inputId="start", label="Minimum position (e.g., 10000-20000)", value = "28477797", width = NULL,
                placeholder = NULL),
      textInput(inputId="end", label="Maximum position (e.g., 10000-20000)", value = "29448354", width = NULL,placeholder = NULL),
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
  filteredsvs <- reactive({
    outputvars= sv %>%
    filter(chr == as.character(input$chr)) %>%
    #filter(SV_type%in%input$SV_type) %>%
    filter(SV_type == "Deletion") %>%
    filter(start >= as.numeric(input$start)) %>%
    filter(end <= as.numeric(input$end)) %>%
    mutate(width = (end-start)+1) #%>%
##    filter(mean_bitscore >= as.numeric(input$scoreMinMax[1])) %>%
##    filter(mean_bitscore <= as.numeric(input$input$scoreMinMax[2]))
  startlist <- outputvars$start
  widthlist <- outputvars$width
  outputvars})
#  grouplist <- outputvars$SV_type
  

#  deletion.track <- AnnotationTrack(data = outputvars,
#                                    start=startlist, 
#                                    width = widthlist, 
#                                    genome="hg19", type="l", 
#                                    name="deletion", 
#                                    window=10,
#                                    #group = rep(grouplist),
#                                    # strand = "+", "-",
#                                    chromosome=input$chr)})
#  
  output$myImage <- renderImage({
    # This file will be removed later by renderImage
    outfile <- tempfile(fileext = '.pdf')
    # generate png
    png(outfile, width = 1200, height = 1000, res = 210)
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr21")
    # build gene models 
    func_chr <- function(chr, start, stop) {
      
      # build chromosome model 
      ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
      # build gene models 
      axisTrack <- GenomeAxisTrack()
      data(geneModels)
      biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = stop, name = "Gene Model", transcriptAnnotation = "symbol")
      # build sv track 
#      deletion.track <- AnnotationTrack(data=dat,
#                                        start=start,
#                                        width = stop,
#                                        genome="hg19", type="l",
#                                        name="deletion",
#                                        window=10,
#                                       #group = rep(grouplist),
#                                        # strand = "+", "-",
#                                       chromosome=chr)
#      # plot tracks
      plotTracks(list(ideoTrack,axisTrack,biomTrack), from = start, to = stop, labelPos = "below")
    }

    func_chr(filteredsvs$chr, filteredsvs$start, filteredsvs$end)
    dev.off()
    
    list(src = outfile,
         contentType = 'image/pdf',
         width = 400,
         height = 300,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
}


shinyApp(ui, server)
  
 