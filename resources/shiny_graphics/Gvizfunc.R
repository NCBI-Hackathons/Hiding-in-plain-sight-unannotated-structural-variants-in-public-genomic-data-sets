# clear workspace 
rm(list=ls())

# load libraries NRIP1
library("ggplot2")
library("Gviz")
library("GenomicRanges")
library("tidyverse")
sv <- read_delim("test.bed", delim = "\t", col_names = TRUE) %>%
  rename("chr" = "#chr")
# remove the # from header column  
# load data into tidyverse

# practice inputs
instart <- 14367028
instop <- 14373374
inchr="chr21"

# filter data based on inputs 
newsv=sv %>%
  filter(chr == inchr) %>%
  filter(SV_type == "Deletion") %>%
  filter(start >= instart) %>%
  filter(end <= instop) %>%
  mutate(width = (end-start)+1) %>%
  mutate(group = 1:length(chr))
startlist <- newsv$start
widthlist <- newsv$width
grouplist <- newsv$group
# build deletion track 
deletion.track <- AnnotationTrack(data = newsv,
                                  start=startlist, 
                                  width = widthlist, 
                                  genome="hg19", type="l", 
                                  name="deletion", 
                                  window=10,
                                  group = rep(grouplist),
                                 # strand = "+", "-",
                                  chromosome=inchr)

plotTracks(deletion.track, scale = 0.5, labelPos = "below", shape="box", groupAnnotation="group")
# plotting function 
func_chr <- function(chr, start, stop) {
  
  # build chromosome model 
    ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = chr)
  # build gene models 
    axisTrack <- GenomeAxisTrack()
    data(geneModels)
    biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = stop, name = "Gene Model", transcriptAnnotation = "symbol")
    plotTracks(list(ideoTrack,axisTrack,biomTrack, deletion.track), from = start, to = stop, scale = 0.5, labelPos = "below")
}

func_chr(inchr, instart, instop)

