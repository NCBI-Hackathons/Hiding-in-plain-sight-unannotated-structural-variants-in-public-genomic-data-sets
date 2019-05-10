library("shiny")
library("shinyWidgets")
library("stringr")
library("ggplot2")
library("Gviz")
library("GenomicRanges")
library("tidyverse")

sv <- read.table("chr21_deletions_annotated.bed", sep="\t", header=TRUE)
sv$mean_bitscore <- as.integer(sv$mean_bitscore)
sv$start <- as.integer(sv$start)
sv$end <- as.integer(sv$end)
sv$width <- sv$end - sv$start + 1

ui <- fluidPage (
  titlePanel("Hiding in Plain Sight"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId="chr", 
                  label="Chromosome",
                  choices=unique(sv$chr)),
      textInput(inputId="posstartstop", label="Minimum/maximum position (e.g., 10000-20000)", value = "", width = NULL,
                placeholder = NULL),
      sliderInput("scoreMinMax", "Score", min(sv$mean_bitscore, na.rm=T), max(sv$mean_bitscore, na.rm=T), c(50, max(sv$mean_bitscore, na.rm=T))),
      checkboxGroupInput("typeInput", "SV type",
                         choices=unique(sv$SV_type),
                         selected=sv$SV_type[1]),
      checkboxGroupInput("sv_seen", "Annotation",
                         choices=names(sv)[8:ncol(sv)],
                         selected=names(sv)[8:ncol(sv)]),
      radioButtons("display", "", choiceNames=c("Table", "Figure"), choiceValues=c("Table", "Figure"), selected="Table", inline=TRUE),
      downloadButton("downloadData", "Download")
      
    ),
    mainPanel( 
      dataTableOutput("datatable"),
      imageOutput("browserimage"),
      textOutput("testtext")
    )
  )
)

server <- function(input, output) {
  filteredsvs <- reactive( {
    if (input$display=="Table") {
      outputvars <- sv
      if (is.integer(input$scoreMinMax[1]) && is.integer(input$scoreMinMax[2])) {
        outputvars <- outputvars[outputvars$mean_bitscore >= input$scoreMinMax[1] & outputvars$mean_bitscore <= input$scoreMinMax[2], ]
      }
      if (input$chr != "") {
        outputvars <- outputvars[outputvars$chr==input$chr,]
      }
      startstop <- get_boundaries()
      if (!is.na(startstop[1]) && !is.na(startstop[2])) {
        outputvars <- outputvars[!((outputvars$end < startstop[1]) | (outputvars$start > startstop[2])),]
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
    }
    else if (input$display=="Figure") {
      outputvars <- data.frame()
    }
    
  })
  
  tidyfilteredsvs <- reactive({
    
    startstop <- get_boundaries()
    if (is.na(startstop[1])) {
      start <- 45000000
    }
    else {
      start <- startstop[1]
    }
    
    if (is.na(startstop[2])) {
      stop <- 47000000
    }
    else {
      stop <- startstop[2]
    }
    outputvars=sv %>%
        #filter(SV_type == "Insertion") %>%
        filter(chr == as.character(input$chr)) %>%
        filter(SV_type%in%input$typeInput) %>%
        filter(start >= as.numeric(start)) %>%
        filter(end <= as.numeric(stop)) %>%
        mutate(width = (end-start)+1)
      #   filter(mean_bitscore >= as.numeric(input$scoreMinMax[1])) %>%
      #  filter(mean_bitscore <= as.numeric(input$input$scoreMinMax[2]))
  })
  
  descriptivetext <- reactive({
    if (input$display=="Figure") {
      startstop <- get_boundaries()
      outputtext <- paste("This is a lovely figure!", startstop[1], startstop[2], sep=" ")
    }
    else if (input$display=="Table") {
      outputtext <- ""
    }
  })

  get_boundaries <- function(){
    if (!is.na(str_extract(input$posstartstop, "(\\S*:){0,1}\\d+\\-\\d+"))) {
      start <- as.integer(sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\2', input$posstartstop))
      end <- as.integer(sub('(\\S*:){0,1}([0-9]+)\\-([0-9]+)', '\\3', input$posstartstop))
    }
    else {
      start = NA
      end = NA
    }
    return(list(start, end))
  }
  
  getimage <- reactive({
    if (input$display=="Figure") { # display browser image
      outfile <- tempfile(fileext = '.pdf')
      png(outfile, width = 1200, height = 1000, res = 210)
      startstop <- get_boundaries()
      if (is.na(startstop[1])) {
        startstop[1] = 47000000
      }
      if (is.na(startstop[2])) {
        startstop[2] = 48000000
      }
      draw_chr(input$chr, as.integer(startstop[1]), as.integer(startstop[2]))
      dev.off()
      list(src = outfile,
           contentType='image/pdf',
           width = 400,
           height = 300,
           alt = "This is alternate text")
    }
    else { # make the image just a tiny empty one
      outfile <- tempfile(fileext = '.pdf')
      png(outfile, width = 1, height = 1)
      dev.off()
      list(src = outfile,
           contentType='image/pdf',
           width = 1,
           height = 1)
    }
  })
  
  draw_chr <- function(chr, start, stop) {
    # build chromosome model 
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
    # build gene models 
    axisTrack <- GenomeAxisTrack()
    data(geneModels)
    biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = stop, name = "Gene Model", transcriptAnnotation = "symbol")
    deletionTrack <- filterandcreatetrack(chr, start, stop)
    
    # plot tracks
    #plotTracks(list(ideoTrack,axisTrack,biomTrack), from = start, to = stop, labelPos = "below")
    plotTracks(list(ideoTrack, axisTrack,biomTrack, deletionTrack), from=start, to=stop, labelPos = "below")
  }
  
  filterandcreatetrack <- function(chr, start, stop) {
    outputvars <- tidyfilteredsvs()
    deletion.track <- AnnotationTrack(data=outputvars,
                                      #start=outputvars$start,
                                      genome="hg19",
                                      type="l",
                                      name="SV",
                                      window=10,
                                      #group = rep(grouplist),
                                      # strand = "+", "-",
                                      chromosome=chr)
  }

  output$downloadData <- downloadHandler(
    filename = function() {
      "svdataset.csv"
    },
    content = function(file) {
      write.csv(filteredsvs(), file, row.names = FALSE)
    }
  )

  output$testtext <- renderText(descriptivetext()) 
  #output$datatable <- renderDataTable(filteredsvs())
  output$datatable <- renderDataTable(tidyfilteredsvs())
  output$browserimage <- renderImage(getimage(), deleteFile = TRUE)

}

shinyApp(ui, server)

