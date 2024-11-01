#! /usr/bin/env Rscript

################################################################################
## Triplet_periodicity.R                                                      ##
## ---------------------                                                      ##
## A script for visual check of triplet periodicity                           ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2022-10-06                                                    ##
## Released under Apache License 2.0                                          ##
################################################################################

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("The script requires at least two input arguments - a filepath to annotations and a dirpath to translated BAM files")
}
if (length(args) > 2) {
  if (length(args) %% 2 == 1) {
    stop(paste("There is an extra argument:",args[length(args)]))
  }
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
      # Output files prefix
      "--output"={ output=args[i+1] },
      # Title of plots (not currently used)
      "--title"={ title=args[i+1] },
      # Outer boundary (in UTR regions)
      "--xout"={ xout=args[i+1] },
      # Inner boundary (in CDS region)
      "--xin"={  xin=args[i+1] },
      # Minimal read length to be plotted
      "--ymin"={ ymin=args[i+1] },
      # Maximal read length to be plotted
      "--ymax"={ ymax=args[i+1] },
      # Whether original version of RiboWaltz should be used (only sufficiently long transcripts are taken https://github.com/LabTranslationalArchitectomics/riboWaltz/blob/d6121ce02f5ab5c92bbcc96d6bba2a74284dcbf4/R/read_end_metaheatmap_plot.R#L65 )
      "--old"={ arg=tolower(args[i+1])
                if (arg=="true" || arg=="t") { old=TRUE
                } else if (arg=="false" || arg=="f") { old=FALSE
                } else { stop(paste("Unrecognized value of --old parameter:",args[i+1])) }
              },
      # How far beyond the visible region is to be taken into account (relevant in old mode only)
      "--overlap"={ overlap=as.integer(args[i+1]) },
      # Output files format
      "--format"={ extension=args[i+1] } )
  }
}

# Set defaults if not overwritten
if (! exists("old")) {
  old=FALSE
}
if (! exists("output")) {
  output=""
}
ylim = c(if(exists("ymin")) as.integer(ymin)-0.5 else NA,
         if(exists("ymax")) as.integer(ymax)+0.5 else NA)
xlim = c(if(exists("xout")) as.integer(xout)+0.5 else 50.5,
         if(exists("xin"))  as.integer(xin) +0.5 else 50.5)
if (! exists("overlap")) {
  overlap=0
}
overlap=overlap-0.5
if (! exists("extension")) {
  extension="svg"
}

# Load libraries
library(riboWaltz)
library(ggplot2)
library(RColorBrewer)

# Read input files
annotation_dt <- create_annotation(gtfpath = args[1], dataSource = "ensemble", organism = "Homo sapiens")
reads_list <- bamtolist(bamfolder = args[2], annotation = annotation_dt)

# Which columns well be expanded
keys <- c("length", "dist", "region", "end")

# Iterate over each BAM file
for (name in names(reads_list)) {
  # Header is not currently used
  if (exists("title")) {
    header = title;
  } else {
    header = name;
  }

  # Extract read counts for each position and length in given ranges
  heatmap <- rends_heat(reads_list, annotation_dt, utr5l=xlim[1]+overlap, cdsl=xlim[2]+overlap, utr3l=xlim[1]+overlap, old=old, sample=name)[["dt"]]

  # To have continuous and the same y-ranges in all plots
  ranges <- sapply(keys, function(h) unique(heatmap[[h]]))
  length.min = if(exists("ymin")) as.integer(ymin) else min(ranges$length)
  length.max = if(exists("ymax")) as.integer(ymax) else max(ranges$length)
  ranges$length=seq(length.min, length.max)
  background <- expand.grid(ranges)
  heatmap <- merge(background, heatmap, by=keys,  all.x=TRUE)
  heatmap[is.na(heatmap)] = 0

  # To trade with zeros in log-scale
  heatmap$count <- heatmap$count+1

  # To have the same z-range in all plots
  rng=range(heatmap$count, na.rm = T)

  # Plotting
  ggplot(heatmap[heatmap$region=="Distance from start (nt)",][heatmap$end=="5' end",], aes(dist, length, fill=count)) + geom_tile() + xlim(-xlim[1], xlim[2]) + scale_y_continuous(limits=ylim) + scale_fill_distiller(palette = "Spectral", trans="log10", limits=rng) + coord_fixed() + xlab("5' end distance from start codon (nt)") + theme(title = element_text(hjust = 0.5))
  ggsave(paste(output,name,"-begin-5.",extension, sep=""))
  ggplot(heatmap[heatmap$region=="Distance from start (nt)",][heatmap$end=="3' end",], aes(dist, length, fill=count)) + geom_tile() + xlim(-xlim[1], xlim[2]) + scale_y_continuous(limits=ylim) + scale_fill_distiller(palette = "Spectral", trans="log10", limits=rng) + coord_fixed() + xlab("3' end distance from start codon (nt)") + theme(title = element_text(hjust = 0.5))
  ggsave(paste(output,name,"-begin-3.",extension, sep=""))
  ggplot(heatmap[heatmap$region=="Distance from stop (nt)",][heatmap$end=="5' end",], aes(dist, length, fill=count)) + geom_tile() + xlim(-xlim[2], xlim[1]) + scale_y_continuous(limits=ylim) + scale_fill_distiller(palette = "Spectral", trans="log10", limits=rng) + coord_fixed() + xlab("5' end distance from stop codon (nt)")  + theme(title = element_text(hjust = 0.5))
  ggsave(paste(output,name,"-end-5.",extension, sep=""))
  ggplot(heatmap[heatmap$region=="Distance from stop (nt)",][heatmap$end=="3' end",], aes(dist, length, fill=count)) + geom_tile() + xlim(-xlim[2], xlim[1]) + scale_y_continuous(limits=ylim) + scale_fill_distiller(palette = "Spectral", trans="log10", limits=rng) + coord_fixed() + xlab("3' end distance from stop codon (nt)")  + theme(title = element_text(hjust = 0.5))
  ggsave(paste(output,name,"-end-3.",extension, sep=""))
}
