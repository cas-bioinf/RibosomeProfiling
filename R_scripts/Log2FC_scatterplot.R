#! /usr/bin/env Rscript

################################################################################
## Log2FC_scatterplot.R                                                       ##
## ---------------------                                                      ##
## A script to generate scatterplots of log2(FoldChange) for pairs of         ##
## knockdowns with highlighted signifficant genes                             ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-10-19                                                    ##
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
  help.common(c("Log2FC_scatterplot.R -h                             Prints this message.",
                "Log2FC_scatterplot.R [options] <mRNA>+ <FP>+ <TE>+  The script requires a list of filepaths to DESeq2",
                "                                                    output stats <deseq>+ and plots a scatterplot of",
                "                                                    log_2(FoldChange) for each pair of files with",
                "                                                    highlighted significantly changed genes.",
                "                                                    For <exclude>+ see --exclude option.",
                "",
                "Options:",
                "  --class <class>      Highlight genes of class <class>:",
                "                       e  exclusive - signifficant change in TE and not signifficantly changed in mRNA.",
                "                       f  forwarded - signifficant change in mRNA and not signifficantly changed in TE.",
                "                       i  intensified - signifficant changes in both mRNA and TE; and the changes are in the same direction.",
                "                       b  buffered - signifficant changes in both mRNA and TE; but the changes are in the opposite direction.",
                "                       Multiple classes can be selected. If only one class is selected, it is divided",
                "                       into Up and Down subclasses. All classes are highlighted by default.",
                "  --header <headers>   Space separated headers corresponding to <deseq>+ files. The headers should have",
                "                       the same order as files in <deseq>+. <deseq>+ filenames are used by default.",
                "  --img {pdf|svg|png}  What file format (and extension) to use to save the plots.",
                "                       Pdf is used by default.",
                "  --label-x <labels>   Labels of x-axis. 'mRNA log_2(FC)' is used by default.",
                "  --label-y <labels>   Labels of y-axis. 'FP log_2(FC)' is used by default.",
                "  --padj <threshold>   Threshold on p_adj to be considered as significantly changed.",
                "                       The default threshold is 0.05.",
                "  --palette <palette>  Space separated colors in HEX format to be used for plots. It must contains",
                "                       either 2 or 3 colours if one class is selected (up, down and optionally none),",
                "                       or n or n+1 colours if more than one class is selected (in the same order with",
                "                       an optionally none-colour as the last one). By default colours #FFC00080,",
                "                       #00A00080, #0000FF80 and #FF00C080 are used for 4 classes; #FF000080, #00FF0080",
                "                       and #0000FF80 are used for 3 classes; and #FF800080 and #5555FF80 are used for",
                "                       1 or 2 classes. The default none-colour is in all cases #80808020.",
                "  --prefix <prefix>    Path prefix for the output files. Empty by default."),
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
         "--img"={ extension=check.extension(args[i+1]) },
         # Space-separated headers for input files (in the same order)
         "--header"={ headers <- unlist(strsplit(args[i+1], split=" ")) },
         # Space-separated headers for input files (in the same order)
         "--label-x"={ label.x <- args[i+1] },
         # Space-separated headers for input files (in the same order)
         "--label-y"={ label.y <- args[i+1] },
         # Classes that should be highlighted
         "--class"={ classes <- tolower(args[i+1]) },
         # Space-separated palette with either 4 colours (Unsignifficant, intersection, 1st only, 2nd only),
         # or n+2 colours (unsignificant, intersection, xth only)
         "--palette"={ palette <- unlist(strsplit(args[i+1], split=" ")) },
         # P_adjusted threshold
         "--padj"={ padj=as.double(args[i+1]) },
         # Default
         {
           if (startsWith(args[i], "--")) {
             help(error = c("Unknown argument:", args[i:(i+1)]))
           } else {
             length = length(args)-i+1
             if (length %% 3 != 0) {
               help(error = c("Unequal length of file sets:", args[i:length(args)]))
             }
             length = length/3
             inputs = args[(i:length(args))]
           }
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
    pairwise.error("headers", headers, "each set of input files", inputs)
  }
}
if (exists("classes")) {
  if (nchar(classes) != length(unique(unlist(strsplit(classes, ""))))) {
    help(error = c("Some class is selected multiple times:", classes))
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
}
if (!exists("headers")) {
  headers=inputs[(2*length+1):(3*length)]
}
if (!exists("label.x")) {
  label.x = "mRNA"
}
if (!exists("label.y")) {
  label.y = "FP"
}
if (!exists("classes")) {
  classes="eifb"
}
if (!exists("padj")) {
  padj=0.05
}
if (!exists("palette", mode="character")) {
  switch (as.character(nchar(classes)),
          "1" =,
          "2" = { palette = c("#FF800080", "#5555FF80") },
          "3" = { palette = c("#FF000080", "#00FF0080", "#0000FF80") },
          "4" = { palette = c("#FFC00080", "#00A00080", "#0000FF80", "#FF00C080") },
          { help(error = c("Unexpected length of classes:", classes)) }
  )
}
if ((nchar(classes)==1 && length(palette)==2) || (nchar(classes)>1 && nchar(classes)==length(palette))) {
  palette = c(palette, "#80808020")
} else if ((nchar(classes)==1 && length(palette)!=3) || (nchar(classes)>1 && (nchar(classes)+1)!=length(palette))) {
  help(error = c("Number of colours in palette does not correspond to the number of selected classes: ", paste(classes, " vs. ", paste(palette, collapse=" "), sep="'")))
}


#### Identify classes that should be highlighted and attach colours
labels = c("Exclusive", "Intensified", "Forwarded", "Buffered")
labels.shortcuts = tolower(substr(labels,1,1))
items = rep("None", 2*length(labels)+1)
items.order = items[1]
for (i in nchar(classes):1) {
  j = match(substr(classes,i,i), labels.shortcuts)
  if (is.na(j)) {
    help(error = c("Unknown class name:", substr(classes,i,i)))
  }
  if (nchar(classes) == 1) {
    items.order = c(paste(labels[j], c("up", "down")), items.order)
    items[(2*j-1):(2*j)] = items.order[1:2]
  } else {
    items.order = c(labels[j], items.order)
    items[(2*j-1):(2*j)] = labels[j]
  }
}


#### Load all input data
load_data = function(paths) {
  return(lapply(paths, read_deseq))
}
data=load_data(inputs)
data.mRNA = data[1:length]
data.FP = data[(length+1):(2*length)]
data.TE = data[(2*length+1):(3*length)]


#### Generate scatter plots for each pair of input files
# Auxiliary functions for easier filtering
construct_header = function(variable, id) { return(paste(variable, id, sep=".")) }
log2fc_header = function(id) { return(construct_header("log2FoldChange", id)) }
padj_header = function(id) { return(construct_header("padj", id)) }
# Main plot function
plot.sample = function(mRNA, FP, TE, header) { 
  # Join the current files and prepare common tables (labels and axis limits)
  subdata = dplyr::inner_join(dplyr::inner_join(mRNA, FP, by="gene_id", suffix=c(".mRNA", "")), TE, by="gene_id", suffix=c(".FP", ".TE"))
  # Remove genes with some undefined padj
  subdata = subdata[rowSums(is.na(subdata[,grep("^padj", colnames(subdata))]))==0,]
  # Range of the examined columns
  limits = range(subdata[,grep("^log2FoldChange", colnames(subdata))])
  # Prepare column with significance of the genes and reorder table by it
  subdata$Signifficant = factor(dplyr::case_when(subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] >  padj &
                                                   sign(subdata[log2fc_header("TE")]) > 0 ~ items[1],
                                                 subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] >  padj &
                                                   sign(subdata[log2fc_header("TE")]) < 0 ~ items[2],
                                                 subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("TE")]) == sign(subdata[log2fc_header("mRNA")]) &
                                                   sign(subdata[log2fc_header("TE")]) > 0 ~ items[3],
                                                 subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("TE")]) == sign(subdata[log2fc_header("mRNA")]) &
                                                   sign(subdata[log2fc_header("TE")]) < 0 ~ items[4],
                                                 subdata[padj_header("TE")] > padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("mRNA")]) > 0 ~ items[5],
                                                 subdata[padj_header("TE")] > padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("mRNA")]) < 0 ~ items[6],
                                                 subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("TE")]) != sign(subdata[log2fc_header("mRNA")]) &
                                                   sign(subdata[log2fc_header("TE")]) > 0 ~ items[7],
                                                 subdata[padj_header("TE")] <= padj & subdata[padj_header("mRNA")] <= padj &
                                                   sign(subdata[log2fc_header("TE")]) != sign(subdata[log2fc_header("mRNA")]) &
                                                   sign(subdata[log2fc_header("TE")]) < 0 ~ items[8],
                                                 TRUE ~ items[9]),
                                items.order)
  # Add number of occurrences of classes to them to be shown in the legend
  items.count = as.data.frame(table(subdata$Signifficant))
  subdata$Signifficant = factor(subdata$Signifficant,
                                labels=paste0(items.count$Var1, " (", items.count$Freq, ")"),
                                items.order)
  # Reorder table to have significant points on the top of the plot
  subdata = subdata[order(subdata$Signifficant, decreasing=TRUE),]
  
  ##TODO
  # Generate the scatterplot
  ggplot(subdata, aes(x=log2FoldChange.mRNA, y=log2FoldChange.FP, col=Signifficant)) +
    geom_point(size=1) +
    theme(axis.line=element_line(colour="black"), panel.background=element_blank(), panel.grid=element_blank()) +
    scale_color_manual(drop=FALSE, values=palette) +
    guides(colour = guide_legend(override.aes=list(size=2, alpha=1))) +
    coord_equal() +
    xlim(limits) + 
    ylim(limits)
  ggsave(paste0(prefix, "log2FoldChange-", header, ".", extension))
  
}
# For each pair of files
for (i in 1:length) {
  plot.sample(data.mRNA[[i]], data.FP[[i]], data.TE[[i]], headers[i])
}
