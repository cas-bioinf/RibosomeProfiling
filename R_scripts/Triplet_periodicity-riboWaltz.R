#! /usr/bin/env Rscript

################################################################################
## Triplet_periodicity.R                                                      ##
## ---------------------                                                      ##
## A script for visual check of triplet periodicity                           ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2025-11-13                                                    ##
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
                "  --img {pdf|svg|png}    What file format (and extension) to use to save the plots.",
                "                         Pdf is used by default.",
                "  --names <old>+ <new>+  Rename samples. With i-th <new> name corresponds to the i-th <old> filename",
                "                         without path and extension (see riboWaltz::bamtolist).",
                "                         This option must be the last one. If not specified, original names are used.",
                "  --prefix <prefix>      Path prefix for output files.",
                "                         If not specified, files are saved in the current working directory.",),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
parity=T
# Parse optional arguments
if (length(args) > 2) {
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
      # Output files prefix
      "--prefix"={ prefix=args[i+1] },
      # Output files format
      "--img"={ extension=check.extension(args[i+1]) },
      # Output files format
      "--names"={
        parity=F
        j=(length(args)+i)/2
        sample_names = setNames(args[(j+1):length(args)], args[(i+1):j])
        } )
  }
}
# Check basic validity of command line arguments
check.parity(args, even=parity)


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
if (!exists("sample_names")) {
  sample_names = NULL
}


#### Initial processing of input files
annotation_dt = riboWaltz::create_annotation(gtfpath = args[1], dataSource = "ensemble", organism = "Homo sapiens")
reads_list = riboWaltz::bamtolist(bamfolder = args[2], annotation = annotation_dt, name_samples=sample_names)
psite_offset = riboWaltz::psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list = riboWaltz::psite_info(reads_list, psite_offset)


#### Generate plots for each input file
for (name in names(reads_list)) {
  example_frames_stratified = riboWaltz::frame_psite_length(reads_psite_list, annotation_dt, sample = name, region = "all", cl = 90)
  ggplot2::ggsave(paste0(prefix, name, "-frame_psite_length.", extension),
                  example_frames_stratified[[paste("plot", name, sep="_")]] + ggplot2::scale_fill_gradientn("% P-sites", colours = c("blue", "white", "red", "black"),
                                                                                                                         breaks = c(0, 33, 67, 100),
                                                                                                                         labels = c(0, 33, 67, 100),
                                                                                                                         limits = c(0, 100)))
  example_frames = riboWaltz::frame_psite(reads_psite_list, annotation_dt, sample = name, region = "all")
  ggplot2::ggsave(paste0(prefix, name, "-frame_psite.", extension), example_frames[[paste("plot", name, sep="_")]])
  example_metaprofile = riboWaltz::metaprofile_psite(reads_psite_list, annotation_dt, sample = name, utr5l = 50, cdsl = 50, utr3l = 50)
  ggplot2::ggsave(paste0(prefix, name, "-metaprofile_psite.", extension), example_metaprofile[[paste("plot", name, sep="_")]] + ggplot2::scale_y_continuous(limits = c(0, NA)), width=30,height=10)
}
