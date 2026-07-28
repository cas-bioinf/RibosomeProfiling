#! /usr/bin/env Rscript

################################################################################
## Triplet_periodicity.R                                                      ##
## ---------------------                                                      ##
## A script for visual check of triplet periodicity                           ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2024-08-01                                                    ##
## Released under Apache License 2.0                                          ##
################################################################################


#### Load a script with commonly used functions to prevent duplication of code (except of this one)
(function(){
  paths = unlist(sapply(sys.frames(), function(f) f$ofile))
  if (length(paths)) {
    path = paths[length(paths)]
  } else {
    args = commandArgs()
    positions = which(startsWith(args, "--file="))
    if (length(positions)) {
      path = substring(args[positions], 8)
    } else {
      path = "."
    }
  }
  path = file.path(dirname(path), "common.R")
  if (file.exists(path)) {
    source(path)
  } else {
    stop("File 'common.R' is missing in the expected location, please repair the path or restore the file.")
  }
})()


#### Prints help message and quit the script for possibility to call it from multiple points
help <- function(error = c()) {
  help.common(c("Triplet_periodicity.R -h                                   Prints this message.",
                "Triplet_periodicity.R <annotations> <directory> [options]  The script generates a heatmap that shows where how many",
                "                                                           reads starts/ ends relative to start/ stop codon depending",
                "                                                           on their lengths.",
                "                                                           Log plus one transformation is used to zero counts can",
                "                                                           be plotted in log scale.",
                "                                                           The script requires a filepath annotations file <annotations>",
                "                                                           in a GTF file format and a path to a directory <directory>",
                "                                                           with BAM files in transcript-relative coordinates.",
                "                                                           Unfortunately, the used library is not able to process an",
                "                                                           individual BAM file only.",
                "",
                "Options:",
                "  --img {pdf|svg|png}  What file format (and extension) to use to save the plots.",
                "                       Pdf is used by default.",
                "  --old {TRUE|FALSE}   Whether to use old - original version of rends_heat from riboWaltz library, with an automatic",
                "                       filtering of too short to be discarded; with the new version, no filtering is performed.",
                "  --overlap <length>   How long region beyond the visible region is to be taken into account for filtering.",
                "                       This option is relevant for old mode only.",
                "                       The default overlap is 0.",
                "  --prefix <prefix>    Path prefix for output files.",
                "                       If not specified, files are saved in the current working directory.",
                "  --xin <length>       How long region within CDS should be shown.",
                "                       The default length is 50.",
                "  --xout <length>      How long region within UTR should be shown.",
                "                       The default length is 50.",
                "  --ymin <length>      Minimal read length that should be plotted.",
                "                       No limit is used by default.",
                "  --ymax <length>      Maximal read length that should be plotted.",
                "                       No limit is used by default."),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=T)
# Parse optional arguments
if (length(args) > 2) {
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
      # Output files prefix
      "--prefix"={ prefix=args[i+1] },
      # Output files format
      "--img"={ extension=check.extension(args[i+1]) },
      # Outer boundary (in UTR regions)
      "--xout"={ xout=as.integer(args[i+1]) },
      # Inner boundary (in CDS region)
      "--xin"={  xin=as.integer(args[i+1]) },
      # Minimal read length to be plotted
      "--ymin"={ ymin=as.integer(args[i+1]) },
      # Maximal read length to be plotted
      "--ymax"={ ymax=as.integer(args[i+1]) },
      # Whether original version of RiboWaltz should be used (only sufficiently long transcripts are taken
      #   https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/d6121ce02f5ab5c92bbcc96d6bba2a74284dcbf4/R/read_end_metaheatmap_plot.R#L65 )
      "--old"={ old = parse.boolean(args[i+1], args[i]) },
      # How far beyond the visible region is to be taken into account (relevant in old mode only)
      "--overlap"={ overlap=as.integer(args[i+1]) } )
  }
}


#### External libraries
# Check whether all used libraries are installed
check.installed("ggplot2", "riboWaltz")
# Load intensively used libraries
library(ggplot2)


#### Set defaults if not overwritten
if (! exists("prefix")) {
  prefix=""
}
if (!exists("extension")) {
  extension = "pdf"
}
ylim = c(if(exists("ymin")) ymin-0.5 else NA,
         if(exists("ymax")) ymax+0.5 else NA)
xlim = c(if(exists("xout")) xout+0.5 else 50.5,
         if(exists("xin"))  xin +0.5 else 50.5)
if (! exists("old")) {
  old = FALSE
}
if (! exists("overlap")) {
  overlap = 0
}
overlap = overlap-0.5


#### Initial processing of input files
annotation_dt = riboWaltz::create_annotation(gtfpath = args[1], dataSource = "ensemble", organism = "Homo sapiens")
reads_list = riboWaltz::bamtolist(bamfolder = args[2], annotation = annotation_dt)


#### Processing of individual input files
## Generate a single plot
plot.triplet.periodicity = function(region, end, lower, upper, suffix) {
  ggplot(heatmap[heatmap$region==region,][heatmap$end==end,], aes(dist, length, fill=count)) +
    geom_tile() +
    xlim(-xlim[lower], xlim[upper]) +
    scale_y_continuous(limits=ylim) +
    scale_fill_distiller(palette = "Spectral", trans="log10", limits=rng) +
    coord_fixed() +
    xlab(paste(end, tolower(region))) +
    theme(title = element_text(hjust = 0.5))
  ggsave(paste0(prefix, name, "-", suffix, ".", extension))
}
# Auxiliary variable of keys for plotted counts
keys = c("length", "dist", "region", "end")
# Generate heatmaps for each BAM file
for (name in names(reads_list)) {
  ## Extract read counts for each position and length in given ranges
  heatmap = riboWaltz::rends_heat(reads_list, annotation_dt, utr5l=xlim[1]+overlap, cdsl=xlim[2]+overlap, utr3l=xlim[1]+overlap, old=old, sample=name)[["count_dt"]]

  ## To have continuous and the same y-ranges in all plots (get ranges, make them continuous, expand the original table, and set missing values to default)
  ranges = sapply(keys, function(h) unique(heatmap[[h]]))
  length.min = if(exists("ymin")) as.integer(ymin) else min(ranges$length)
  length.max = if(exists("ymax")) as.integer(ymax) else max(ranges$length)
  ranges$length=seq(length.min, length.max)
  background = expand.grid(ranges)
  heatmap = merge(background, heatmap, by=keys,  all.x=TRUE)
  heatmap[is.na(heatmap)] = 0

  ## To trade with zeros in log-scale
  heatmap$count = heatmap$count+1

  ## To have the same z-range in all plots
  rng=range(heatmap$count, na.rm = T)

  # Plotting
  plot.triplet.periodicity("Distance from start (nt)", "5' end", 1, 2, "begin-5")
  plot.triplet.periodicity("Distance from start (nt)", "3' end", 1, 2, "begin-3")
  plot.triplet.periodicity("Distance from stop (nt)" , "5' end", 2, 1, "end-5")
  plot.triplet.periodicity("Distance from stop (nt)" , "3' end", 2, 1, "end-3")
}
