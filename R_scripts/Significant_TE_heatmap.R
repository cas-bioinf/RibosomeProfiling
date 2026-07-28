#! /usr/bin/env Rscript

################################################################################
## Significant_TE_heatmap.R                                                   ##
## ---------------------                                                      ##
## A script to generate a heatmap showing TE of significantly changed genes   ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-09-27                                                    ##
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
  help.common(c("Significant_TE_heatmap.R -h                  Prints this message.",
                "Significant_TE_heatmap.R [options] <input>+  The script generates a heatmap that shows TE of significantly changed genes.",
                "                                             The script requires at least one filepath to a tab-separated table <input>",
                "                                             ...",
                "",
                "Options:",
                "  --bins <bins>        A number of bins that should be used to split data. The default number of bins is 255.",
                "  --header <headers>   Space separated headers corresponding to <input>+ files. The headers should have",
                "                       the same order as files in <input>+. <input>+ filenames are used by default.",
                "  --img {pdf|svg|png}  What file format (and extension) to use to save the plots.",
                "                       Pdf is used by default.",
                "  --max <max>          The upper limit of bins. -<max> is used for the lower limit. The default boundary is 2.",
                "  --output <path>      Where the ouptut plot should be saved.",
                "                       If <path> is an existing directory, <path>/Significant_TE_heatmap.[file_format] is used",
                "                       as the filename. If file format extension is missing, it is appended. If path is not",
                "                       specified, the current working directory is used as <path>.",
                "  --padj <threshold>   Space-separated thresholds on p_adj to be considered as significantly changed.",
                "                       The default threshold is 0.05."),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Parse optional arguments
for (i in seq(1, length(args), by=2)) {
  switch(args[i],
    # Output file
    "--output"={ output=args[i+1] },
    # File format of output histograms
    "--img"={ extension=check.extension(args[i+1]) },
    # Space-separated headers for input files (in the same order)
    "--header"={ headers <- unlist(strsplit(args[i+1], split=" ")) },
    # P_adjusted threshold
    "--padj"={ padj=as.double(args[i+1]) },
    # Number of bins to visualize
    "--bins"={ bins=as.integer(args[i+1]) },
    # Maximal (and -minimal) value to be displayed 
    "--max"={ bound=as.double(args[i+1]) },
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


#### External libraries
# Check whether all used libraries are installed
check.installed("dplyr", "gplots", "tools", "RColorBrewer")


#### Set defaults if not overwritten
output = construct.filename(output, extension, "Significant_TE_heatmap")
extension = output$extension
output = output$filepath
length = length(inputs)
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


#### Initial processing of input files
# Read all input files; rename all columns that are not used for join to ensure uniqueness of column names; and finaly join by 'gene_id' column
data = Reduce(function(x, y) dplyr::inner_join(x, y, by="gene_id"),
              sapply(headers, function(x, data) setNames(data[[x]], gsub("^(?!gene_id$)(.*)$", paste("\\1", x, sep="."), names(data[[x]]), perl=T)),
                     data=setNames(lapply(inputs, read_deseq), headers), simplify=F))
# Filter genes that are significant for at least one column
data=data[which(rowSums(data[grepl("^padj.", colnames(data))]<padj, na.rm=TRUE)>0),]
# Select matrix of log_2 fold changes
rownames(data)=data[,"gene_id"]
log2foldchange=as.matrix(data[grepl("^log2FoldChange.", colnames(data))])
colnames(log2foldchange)=headers


#### Plot the heatmap
colors=colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(bins)
get(extension)(output)
gplots::heatmap.2(log2foldchange, trace="none", col=colors, labRow=F, breaks = seq(-bound, bound, length.out = bins+1), key.xlab="log2FoldChange")
dev.off()
