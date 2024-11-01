#! /usr/bin/env Rscript

################################################################################
## Significant_TE_intersection.R                                              ##
## ---------------------                                                      ##
## A script to generate scatterplots of log2(FoldChange) for pairs of         ##
## knockdowns with highlighted signifficant genes                             ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-03-22                                                    ##
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
  help.common(c("Significant_TE_intersection.R -h                               Prints this message.",
                "Significant_TE_intersection.R [options] <deseq>+ [<exclude>+]  The script requires a list of filepaths to DESeq2",
                "                                                               output stats <deseq>+ and plots a scatterplot of",
                "                                                               log_2(FoldChange) for each pair of files with",
                "                                                               highlighted significantly changed genes.",
                "                                                               For <exclude>+ see --exclude option.",
                "",
                "Options:",
                "  --exclude {TRUE|FALSE}               Whether there would be provided an additional list of DESeq2",
                "                                       output statswhose significantly changed genes should not be",
                "                                       considered as significantly changed in <deseq>+. List of these",
                "                                       files should be provided as the last, should have the same",
                "                                       length as <deseq>+ and the i-th item of <exclude>+ corresponds",
                "                                       to the i-th item of <deseq>+.",
                "  --header <headers>                   Space separated headers corresponding to <deseq>+ files.",
                "                                       The headers should have the same order as files in <deseq>+.",
                "                                       <deseq>+ filenames are used by default.",
                "  --img {pdf|svg|png}                  What file format (and extension) to use to save the plots.",
                "                                       Pdf is used by default.",
                "  --padj <threshold>                   Threshold on p_adj to be considered as significantly changed.",
                "                                       The default threshold is 0.05.",
                "  --palette <palette>                  Space separated colors in HEX format to be used for plots.",
                "                                       It must contains either 4 colours (unsignifficant,",
                "                                       intersection, first only, second only), i.e. the same palette",
                "                                       for each plot; or n+2 colours (unsignifficant, intersection,",
                "                                       i-th file only...), i.e. constant color for <deseq> files",
                "                                       accross all relevant plots. By default, GREY, RED, GREEN, BLUE",
                "                                       is used with with 0.125 alpha channel for the GREY colour and",
                "                                       0.5 for the rest.",
                "  --prefix <prefix>                    Path prefix for the output files. Empty by default."),
              error)
}


#### Parse arguments
args = get.arguments()
# Parse functional arguments
for (i in seq(1, length(args), by=2)) {
  switch(args[i],
    # Prefix for output files
    "--prefix"={ prefix=args[i+1] },
    # Plots file format (supported are pdf, svg and png)
    "--img"={ extension=args[i+1] },
    # Space-separated headers for input files (in the same order)
    "--header"={ headers <- unlist(strsplit(args[i+1], split=" ")) },
    # Space-separated palette with either 4 colours (Unsignifficant, intersection, 1st only, 2nd only),
    # or n+2 colours (unsignificant, intersection, xth only)
    "--palette"={ colours <- unlist(strsplit(args[i+1], split=" ")) },
    # P_adjusted threshold
    "--padj"={ padj=as.double(args[i+1]) },
    # Whether an additional set with its significant genes to be considered as unchanged will be provided
    "--exclude"={ exclude = parse.boolean(args[i+1], "exclude") },
    # Default
    {
      if (startsWith(args[i], "--")) {
        help(error = c("Unknown argument:", args[i:(i+1)]))
      } else if (exists("exclude") && exclude) {
        count = length(args)+1-i
        if (length(args) %% 2 != 0) {
          stop(paste("Exclude set does not have the same size as the main set:", paste0(args[(i:length(args))], collapse="; ")))
        }
        count = count/2-1
        inputs = args[i:(i+count)]
        excludes = args[(length(args)-count):length(args)]
      } else {
        exclude = F
        inputs = args[(i:length(args))]
      }
      length = length(inputs)
      break
    }
  )
}
# Check basic validity of command line arguments
if (! exists("inputs")) {
  help(error = c("No input file is provided:", args))
}
if (exists("headers")) {
  if (length(headers) != length) {
    pairwise.error("headers", headers, "input files", inputs)
  }
}
if (exists("colours", mode="character")) {
  if (length(colours) == 4) {
    palette = c(colours[2], colours[3], colours[4], colours[1])
  }
  else if (length(colours)-2 == length) {
    individual_palettes = T
  } else {
    pairwise.error("colours", colours, "input files", inputs,
                   "The number of colours must be either 4, or the {number of input files}+2:")
  }
}


#### External libraries
# Check whether all used libraries are installed
check.installed("dplyr", "ggplot2")
# Load intensively used libraries
library(ggplot2)


#### Set defaults if not overwritten
if (!exists("prefix")) {
  prefix = ""
}
if (!exists("extension")) {
  extension = "pdf"
} else {
  check.extension(extension)
}
if (!exists("headers")) {
  headers=inputs
}
if (!exists("colours", mode="character")) {
  palette = c("#FF000080", "#00FF0080", "#0000FF80", "#80808020")
}
if (!exists("padj")) {
  padj=0.05
}


#### Load all input data
load_data = function(paths) {
  return(lapply(paths, read.csv, sep='\t'))
}
data = load_data(inputs)
if (exclude) {
  data.exclude = load_data(excludes)
}


#### Generate scatter plots for each pair of input files
# Auxiliary functions for easier filtering
construct_header = function(variable, id) { return(paste(variable, id, sep=".")) }
log2fc_header = function(id) { return(construct_header("log2FoldChange", id)) }
padj_header = function(id) { return(construct_header("padj", id)) }
# For each pair of files
for (i in 1:(length-1)) {
  for (j in (i+1):length) {
    # Join the current files and prepare common tables (lables and axis limits)
    subdata = dplyr::inner_join(data[[i]], data[[j]], by="gene_id",
                                suffix=c(paste0(".",headers[i]), paste0(".",headers[j])))
    # Classes
    items = c("Both", headers[i], headers[j], "None")
    # Range of the examined columns
    limits = range(subdata[,grep("^log2FoldChange", colnames(subdata))])
    # Generate a palette if it should be customized for each pair
    if (exists("individual_palettes", mode="logical") && individual_palettes) {
      palette = c(colours[2], colours[i+2], colours[j+2], colours[1])
    }
    # Prepare column with significance of the genes and reorder table by it
    if (exclude) {
      subdata$Signifficant = factor(dplyr::case_when(subdata[padj_header(headers[i])] <= padj & data.exclude[[i]]$padj > padj &
                                                     subdata[padj_header(headers[j])] <= padj & data.exclude[[j]]$padj > padj ~ items[1],
                                                     subdata[padj_header(headers[i])] <= padj & data.exclude[[i]]$padj > padj ~ items[2],
                                                     subdata[padj_header(headers[j])] <= padj & data.exclude[[j]]$padj > padj ~ items[3],
                                                     TRUE ~ items[4]),
                                    items)
    } else {
      subdata$Signifficant = factor(dplyr::case_when(subdata[padj_header(headers[i])] <= padj &
                                                     subdata[padj_header(headers[j])] <= padj ~ items[1],
                                                     subdata[padj_header(headers[i])] <= padj ~ items[2],
                                                     subdata[padj_header(headers[j])] <= padj ~ items[3],
                                                     TRUE ~ items[4]),
                                    items)
    }
    # Add number of occurrences of classes to them to be shown in the legend
    items.count = as.data.frame(table(subdata$Signifficant))
    subdata$Signifficant = factor(subdata$Signifficant,
                                  items,
                                  labels=paste0(items.count$Var1, " (", items.count$Freq, ")"))
    # Reorder table to have significant points on the top of the plot
    subdata = subdata[order(subdata$Signifficant, decreasing=TRUE),]
    # Generate the scatterplot
    ggplot(subdata, aes(x=!!sym(paste("log2FoldChange", headers[i], sep=".")),
                        y=!!sym(paste("log2FoldChange", headers[j], sep=".")),
                        col=Signifficant)) +
      geom_point(size=1) +
      theme(axis.line=element_line(colour="black"), panel.background=element_blank(), panel.grid=element_blank()) +
      scale_color_manual(drop=FALSE, values=palette) +
      guides(colour = guide_legend(override.aes=list(size=2, alpha=1))) +
      coord_equal() +
      xlim(limits) + 
      ylim(limits)
    ggsave(paste0(prefix, "log2FoldChange-", headers[i], "_", headers[j], ".", extension))
  }
}
