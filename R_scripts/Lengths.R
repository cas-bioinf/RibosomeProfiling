#! /usr/bin/env Rscript

################################################################################
## Lengths.R                                                                  ##
## ---------------------                                                      ##
## A script to compute read length frequencies and plot their histograms      ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-08-24                                                    ##
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
  help.common(c("Lengths.R -h                   Prints this message.",
                "Lengths.R <lengths> [options]  The script computes read length frequencies and plots corresponding histograms.",
                "                               The script requires a filepath to a tab-separated table <lengths> with absolute",
                "                               numbers of occurances of individual read lengths.",
                "",
                "Options:",
                "  --img {pdf|svg|png}  What file format (and extension) to use to save the plots.",
                "                       Pdf is used by default.",
                "  --prefix <prefix>    Path prefix for the output files. Empty by default."),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=F)
# Parse optional arguments
if (length(args) > 1) {
  for (i in seq(2, length(args), by=2)) {
    switch(args[i],
      # Prefix path of output files
      "--prefix"={ prefix=args[i+1] },
      # File format of output histograms
      "--img"={ extension=check.extension(args[i+1]) },
      # Default
      { stop(paste("Unrecognized argument:", args[i])) }
    )
  }
}


#### External libraries
# Check whether all used libraries are installed
check.installed("ggplot2", "stringr")
# Load intensively used libraries
library(ggplot2)


#### Set defaults if not overwritten
if (! exists("prefix")) {
  prefix = ""
}
if (!exists("extension")) {
  extension = "pdf"
}


#### Initial processing of input files
lengths = read.csv(args[1], sep='\t', row.names=1)
# Sort columns
lengths = lengths[,order(colnames(lengths))]
lengths.rel = t(t(lengths)/colSums(lengths,na.rm=T))


#### Generate stats and histograms of read length
process = function(joined, pattern, filename, title) {
  ## Identify columns to be examined
  header = grepl(pattern, colnames(lengths.rel))
  
  ## Compute stats and ommit NA rows (to trim lengths unused in the current subset)
  mean = rowSums(lengths.rel[,header]) / sum(header)
  sd = apply(lengths.rel[,header], 1, sd)
  table = na.omit(data.frame(mean, sd))
  
  ## Save the histogram
  ggplot(table) +
    geom_bar(aes(x=as.numeric(rownames(table)), y=mean), stat="identity") +
    geom_errorbar(aes(x=as.numeric(rownames(table)), ymin=mean-sd, ymax=mean+sd)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=paste("Read length distribution of",title), x="Length (nt)", y="Count frequency (%)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(prefix, filename, ".", extension))
  
  ## Add current stats to other
  colnames(table) <- paste(title, colnames(table))
  return(transform(merge(joined,table,by=0,all=T,sort=F), row.names=Row.names, Row.names=NULL))
}
# Initialization of the stats table
joined=data.frame(matrix(ncol=0,nrow=length(rownames(lengths))))
rownames(joined)=rownames(lengths)
# Process each examined group
for (type in c("FP", "mRNA")) {
  for (unit in unique(stringr::str_split(colnames(lengths), "_", n = 3, simplify = T)[,2])) {
    joined = process(joined, paste(type,unit,"",sep="_"), paste(type,unit,sep="-"), paste(type,unit))
  }
  joined = process(joined, paste(type,"",sep="_"), type, type)
}
# Write stats to the output file
write.csv(joined[order(as.numeric(rownames(joined))),], file=paste0(prefix,"read_length_distribution.csv"))
