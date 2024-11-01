#! /usr/bin/env Rscript

################################################################################
## Lengths.R                                                                  ##
## ---------------------                                                      ##
## A script to compute read length frequencies and plot their histograms      ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2022-09-14                                                    ##
## Released under Apache License 2.0                                          ##
################################################################################

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("The script requires at least one input argument - path to a table with read length counts for each sample")
}
if (length(args) > 1) {
  if (length(args) %% 2 == 0) {
    stop(paste("There is an extra argument:",args[-1]))
  }
  for (i in seq(2, length(args), by=2)) {
    switch(args[i],
      # Prefix path of output files
      "--output"={ prefix=args[i+1] },
      # File format of output histograms
      "--img"={ extension=args[i+1] },
      # Default
      { stop(paste("Unrecognized argument:", args[i])) }
    )
  }
}

# Set defaults if not overwritten
if (! exists("prefix")) {
  prefix = ""
}
if (!exists("extension")) {
  extension = "pdf"
}
if (! any(extension == c("pdf","png","svg"))) {
  stop(paste("Extension '",extension,"' is not supported, please use 'pdf' (default), 'svg', or 'png'"))
}

# Load libraries
library(ggplot2)
library(stringr)

# Read and prepare input file
lengths = read.csv(args[1],sep='\t')
rownames(lengths) = lengths[,'Length']
lengths = lengths[,colnames(lengths)!="Length"]
lengths = lengths[,order(colnames(lengths))]
lengths.rel = t(t(lengths)/colSums(lengths,na.rm=T))

# Generate stats of read lengths and histograms
process = function(joined, pattern, filename, title) {
  # Compute stats
  header = grepl(pattern, colnames(lengths.rel))
  mean = rowSums(lengths.rel[,header]) / sum(header)
  sd = apply(lengths.rel[,header], 1, sd)
  table = na.omit(data.frame(mean, sd))
  # Save the histogram
  get(extension)(paste0(prefix, filename, ".", extension))
  print(ggplot(table) +
    geom_bar( aes(x=as.numeric(rownames(table)), y=mean), stat="identity") +
    geom_errorbar( aes(x=as.numeric(rownames(table)), ymin=mean-sd, ymax=mean+sd)) +
    scale_y_continuous(labels = scales::percent) +
    labs(title=paste("Read length distribution of",title), x="Length (nt)", y="Count frequency (%)") +
    theme(plot.title = element_text(hjust = 0.5)))
  dev.off()
  # Add current stats to other
  colnames(table) <- paste(title,colnames(table))
  return(transform(merge(joined,table,by=0,all=T,sort=F), row.names=Row.names, Row.names=NULL))
}
joined=data.frame(matrix(ncol=0,nrow=length(rownames(lengths))))
rownames(joined)=rownames(lengths)
for (type in c("FP", "mRNA")) {
  for (unit in unique(str_split(colnames(lengths), "_", n = 3, simplify = T)[,2])) {
    joined = process(joined, paste(type,unit,"",sep="_"), paste(type,unit,sep="-"), paste(type,unit))
  }
  joined = process(joined, paste(type,"",sep="_"), type, type)
}
# Write summary stats to the output file
write.csv(joined[order(as.numeric(rownames(joined))),], file=paste0(prefix,"read_length_distribution.csv"))
