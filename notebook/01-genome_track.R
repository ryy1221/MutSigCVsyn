### This script plot genome track according to gencode v19 annotation file for local mutation heterogeneity plot
rm(list=ls())
setwd('/gpfs/group/epo2/default/yur97/gitlab/pcawg-to-mutsigcv/notebook')

### Load pakcages
library(Gviz)
library(GenomicRanges)

gtrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr8")
plotTracks(ideoTrack, from = 0, to = 1000000
           showId = FALSE, showBandId = TRUE, cex.bands = 0.4)

data(cpgIslands)
cpgIslands
atrack <- AnnotationTrack(cpgIslands, name = "CpG")
plotTracks(atrack)
