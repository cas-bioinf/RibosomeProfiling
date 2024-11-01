#! /usr/bin/env Rscript

################################################################################
## Significant_TE_heatmap.R                                                   ##
## ---------------------                                                      ##
## A script to generate a heatmap showing TE of significantly changed genes   ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2022-09-21                                                    ##
## Released under Apache License 2.0                                          ##
################################################################################

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("The script requires at least one input argument - input file(s) (results from DESeq2)")
}
for (i in seq(1, length(args), by=2)) {
  switch(args[i],
    # Output file
    "--output"={ output=args[i+1] },
    # Space-separated headers for input files (in the same order)
    "--header"={ headers <- unlist(strsplit(args[2], split=" ")) },
    # P_adjusted threshold
    "--padj"={ padj=as.double(args[i+1]) },
    # Number of bins to visualize
    "--bins"={ bins=as.integer(args[i+1]) },
    # Maximal (and -minimal) value to be displayed 
    "--max"={ bound=as.integer(args[i+1]) },
    # Default
    {
      if (startsWith(args[i], "--")) {
        stop(paste("Unrecognized argument:", args[i]))
      } else {
        inputs=args[(i:length(args))]
        break
      }
    }
  )
}
if (! exists("inputs")) {
  stop("No input file is provided")
}


# Load libraries
library(dplyr)
library(gplots)
library(tools)
library(RColorBrewer)


# Set defaults if not overwritten
length=length(inputs)
if (exists("output")) {
  extension=tolower(file_ext(output))
  if (nchar(extension) == 0 || (extension != "pdf" && extension != "png" && extension != "svg" && extension != "tif" && extension != "tiff")) {
    save="pdf"
    output=paste(args[1],"pdf",sep=".")
  } else {
    save=extension
  }
} else {
  save="pdf"
  output="Significant_TE_heatmap.pdf"
}
if (exists("headers")) {
  if (length(headers) != length) {
    stop("The number of headers must be the same as the number of input files")
  }
} else {
  headers=inputs
}
if (! exists("padj")) {
  padj=0.05
}
if (! exists("bins")) {
  bins=255
}
if (! exists("bound")) {
  bound=2
}


# Join all data (first column is used for initialization; last column is added separately to easier adding of identifying suffix)
data=read.csv(inputs[1],sep='\t')
for (i in 1:(length-2)) {
  data=inner_join(data,read.csv(inputs[i+1],sep='\t'),by="gene_id",suffix=c(paste0(".",headers[i]),""))
}
data=inner_join(data,read.csv(inputs[length],sep='\t'),by="gene_id",suffix=c(paste0(".",headers[length-1]),paste0(".",headers[length])))

# Filter genes that are signifficant for at least one col
data=data[which(rowSums(data[grepl("^padj.", colnames(data))]<padj, na.rm=TRUE)>0),]

# Select matrix of log_2 fold changes
rownames(data)=data[,"gene_id"]
log2foldchange=data[grepl("^log2FoldChange.", colnames(data))]
colnames(log2foldchange)=headers
log2foldchange=as.matrix(log2foldchange)

# Plot heatmap
colors=colorRampPalette(rev(brewer.pal(11,"RdBu")))(bins)
get(save)(output)
heatmap.2(log2foldchange, trace="none", col=colors, labRow=F, breaks = seq(-bound, bound, length.out = bins+1), key.xlab="log2FoldChange")
dev.off()
