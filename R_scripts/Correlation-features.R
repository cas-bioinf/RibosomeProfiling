#!/usr/bin/env Rscript

################################################################################
## Correlation-features.R                                                     ##
## ---------------------                                                      ##
## A script to evaluate correlations between log_2(FoldChange) and RNA        ##
## characteristics, and to generate boxplots, and kernell density plots or    ##
## bar plots to visualize changes of the RNA characteristics in               ##
## up/downregulated genes                                                     ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-04-28                                                    ##
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
  help.common(c("Correlation-features.R -h                         Prints this message.",
                "Correlation-features.R <stats> <deseq> [options]  The script requires a filepath to DESeq2 output stats <deseq> and",
                "                                                  filepath to background statistics <stats> and create correlation plots.",
                "",
                "Options:",
                "  --exclude <path>                     Path to DESeq2 output stats whose significantly changed genes",
                "                                       should not be considered as significantly changed in <stats>.",
                "  --img {pdf|svg|png}                  What file format (and extension) to use to save the plots.",
                "                                       Pdf is used by default.",
                "  --method {pearson|spearman|kendall}  What method to use to compute correlation coefficients.",
                "                                       Spearman's rank correlation coefficient is used by default.",
                "  --outliers {TRUE|FALSE}              Whether to plot outliers in boxplots. They are plotted by default.",
                "  --padj <threshold>                   Threshold on p_adj to be considered as significantly changed.",
                "                                       The default threshold is 0.05.",
                "  --prefix <prefix>                    Path prefix for the generated plots. Empty by default.",
                "  --table <output>                     Do not generate plots, just save the result table to <output>.",
                "                                       Options other than --exclude are not used with this option."),
             error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=T)
# Parse optional arguments
if (length(args) > 2) {
  # All switchers are in format --keyword value, so this is possible
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
           # Path prefix for output files
           "--prefix"={ prefix=args[i+1] },
           # Plots file format (supported are pdf, svg and png)
           "--img"={ extension=args[i+1] },
           # P_adjusted threshold
           "--padj"={ threshold=as.double(args[i+1]) },
           # Correlation coeficient to be evaluated
           "--method"={ method=args[i+1] },
           # What significant genes should be considered as unchanged
           "--exclude"={ excludePath = args[i+1] },
           # Whether outliers should be plotted
           "--outliers"={ outliers = parse.boolean(args[i+1], "outliers") },
           # Just save the final table in this file
           "--table"={ output.table = args[i+1] },
           # Unknown switcher
           { help(error = c("Unknown argument:", args[i:(i+1)])) }
    )
  }
}
# Assign mandatory arguments for easier changes
maneStatsPath = args[1]
deseqPath     = args[2]

#### Set defaults if not overwritten
if (!exists("prefix")) {
  prefix = ""
}
if (!exists("extension")) {
  extension = "pdf"
} else {
  check.extension(extension)
}
if (!exists("threshold")) {
  threshold = 0.05
}
if (!exists("method")) {
  method="spearman"
}
if (!exists("outliers")) {
  outliers=TRUE
}


#### External libraries
# Check whether all used libraries are installed
check.installed("dplyr", "ggplot2", "stringr")
# Load intensively used libraries
library(ggplot2)


#### Generation of plots
# Evaluate p-value
significancy = function(pvalue) {
  if      (is.na(pvalue)) return("")
  else if (pvalue<1e-10)  return("***")
  else if (pvalue<1e-5)   return("**")
  else if (pvalue<0.05)   return("*")
  else                    return("")
}
# Generate kernell density plots, resp. bar plots if there should be only small range of integer x-axis
plot_density = function(data, direction, description.direction, column, description="",
                        logarithmic=F, percentage=F, smallint=F, prefix="") {
  if(length(unique(data[,column])) > 1 && length(which(data[,direction] == "Background")) > 1
                                       && length(which(data[,direction] != "Background")) > 1) {
    # Compute and label p-value
    Pval = t.test(reformulate(direction, if(logarithmic) paste0('log(',column,')') else column), data=data )$p.value
    label = trimws(paste("P =", formatC(Pval,digits=2,format="G"),significancy(Pval)))
    # Plot the plot
    ggplot(data, aes(!!sym(column), fill=!!sym(direction), colour=!!sym(direction)))+
      {if(smallint) geom_bar(aes(y=after_stat(prop)), position=position_dodge(preserve="single"))
        else geom_density(alpha=0.25)}+
      {if(logarithmic)scale_x_log10()}+
      {if(percentage)scale_x_continuous(labels=scales::percent_format(accuracy=1))}+
      {if(smallint)scale_y_continuous(labels=scales::percent_format(accuracy=1))}+
      annotate("text", x=Inf, y=Inf, label=label, size=8, vjust=1, hjust=1) +
      labs(x=description, y={if(smallint) "Frequency" else "Density"}, fill="Legend", colour="Legend")
    ggsave(paste0(prefix, column, "-", direction, "-", {if(smallint)"bar" else "density"}, ".", extension))
  }
}
# Generate boxplot
plot_box = function(data, column, description, logarithmic, percentage, prefix="") {
  # No changes, no plot needed
  if (sum(data$Direction=="Unchanged") == nrow(data)) { return() }
  # Compute p-values; ranked test, so log can be ommited
  if (sum(data$Direction=="Upregulated") == 0) {
    Pval = c(1,
             1,
             wilcox.test(x = data[data$Direction=="Downregulated",column],
                         y = data[data$Direction=="Unchanged"    ,column])$p.value)
  } else if (sum(data$Direction=="Downregulated") == 0) {
    Pval = c(wilcox.test(x = data[data$Direction=="Upregulated",column],
                         y = data[data$Direction=="Unchanged"  ,column])$p.value,
             1,
             1)
  } else {
    Pval=c(wilcox.test(x = data[data$Direction=="Upregulated"  ,column],
                       y = data[data$Direction=="Unchanged"    ,column])$p.value,
           wilcox.test(x = data[data$Direction=="Upregulated"  ,column],
                       y = data[data$Direction=="Downregulated",column])$p.value,
           wilcox.test(x = data[data$Direction=="Downregulated",column],
                       y = data[data$Direction=="Unchanged"    ,column])$p.value)
  }
  # Evaluate p-values
  label = sapply(Pval, significancy)
  # Set colors and positions of labels
  if (sum(data$Direction=="Upregulated") == 0) {
    colors = c("#0000FF", "#FFFFFF")
    left   = c("",        label[3])
    right  = c(label[2],  label[1])
  } else if (sum(data$Direction=="Downregulated") == 0) {
    colors = c("#FFFFFF", "#FF0000")
    left   = c(label[3],  label[2])
    right  = c(label[1],  "")
  } else {
    colors = c("#0000FF", "#FFFFFF", "#FF0000")
    left   = c("",        label[3],  label[2])
    right  = c(label[2],  label[1],  "")
  }
  # Plot the boxplot
  get(extension)(paste0(prefix,column,"-box.",extension))
  boundaries = boxplot(reformulate("Direction", if(percentage) paste("100",column,sep="*") else column), data=data,
                       log=if(logarithmic) "y" else "", xlab=NULL, ylab=description, col=colors, outline=outliers)
  text(x=c(1:3), y=boundaries$stats[nrow(boundaries$stats)-1,], left,  col="#0000FF", adj=c(1,0), cex=3)
  text(x=c(1:3), y=boundaries$stats[nrow(boundaries$stats)-1,], right, col="#FF0000", adj=c(0,0), cex=3)
  dev.off()
}
# Generate plots and evaluate the correlation
process_feature = function(data, column, description, logarithmic, percentage, smallint, prefix="") {
  # Common filter to do not repeat it each times - remove NA, and zeros if x-axis should be in logarithmic scale
  data.na = data[which(!is.na(data[,column]) & (!logarithmic | data[,column]>0)),]
  # Plot densities/ frequencies
  plot_density(data.na, "Up",   "upregulated",   column, description, logarithmic, percentage, smallint, prefix)
  plot_density(data.na, "Down", "downregulated", column, description, logarithmic, percentage, smallint, prefix)
  # Plot boxplots
  plot_box(data.na, column, description, logarithmic, percentage, prefix)
  # Compute correlation coefficents and their p-values, and return them
  cor = data.frame(cor.test(data.na$log2FoldChange, data.na[,column], method=method)[c("estimate","p.value")],
                   stringsAsFactors=FALSE)
  colnames(cor)[1] = rownames(cor)[1]
  rownames(cor) = column
  return(cor)
}

#### Read input files and create the final table
# Read DESeq2 output file and repair its first header if missing
read_deseq = function(path) {
  file = read.csv(path, sep='\t')
  if (colnames(file)[1] == "X") {
    colnames(file)[1] = "gene_id"
  }
  return(file)
}
# Read input files
maneStats = read.csv(maneStatsPath, sep='\t')
deseq = read_deseq(deseqPath)
if (exists("excludePath")) {
  exclude = read_deseq(excludePath)
  if (! all(deseq$gene_id == exclude$gene_id)) {
    stop(paste("Files ", args[2], " and ",args[i+1],
               " does not contain the same genes or they are in a different order.", sep="'"))
  }
}
# Filter out Ids without defined p_adj (e.g. low mean normalized count i.e. its qualified decision cannot be made)
data = deseq[!is.na(deseq$padj),]
if (exists("exclude")) {
  data = subset(data, gene_id %in% exclude[!is.na(exclude$padj),]$gene_id)
}
# Classify up- and down-regulated genes
data$Direction = ifelse(data$padj>threshold, "Unchanged", ifelse(data$log2FoldChange>0, "Upregulated", "Downregulated"))
if (exists("exclude")) {
  data[data$gene_id %in% subset(exclude,padj<=threshold)$gene_id,]$Direction = "Unchanged"
}
data$Up =   ifelse(data$Direction ==   "Upregulated",   "Upregulated", "Background")
data$Down = ifelse(data$Direction == "Downregulated", "Downregulated", "Background")
# Construct a table with all informations about examined genes only
data = merge(data, maneStats, by="gene_id")


#### Just save the the final table to process it elsewhere
if (exists("output.table")) {
  write.table(data, output.table, sep="\t", row.names=FALSE)
  quit()
}


#### Determine settings for each examined feature
# Generate labels from sufixes of uORFs detail
process_suffix = function(detail) {
  return(dplyr::case_when(detail == ""          ~ "",
                          detail == "before"    ~ "and classified as non-overlapping",
                          detail == "across"    ~ "and classified as overlapping",
                          detail == "extension" ~ "and classified as N-terminal extension",
                          TRUE                  ~ gsub('_',' ',detail)))
}
# Set settings for all features
settings = setNames(data.frame(colnames(maneStats)[-1],stringsAsFactors=F), "header")
settings[c("feature", "detail")] = stringr::str_split_fixed(settings$header, '\\.', 2)
settings$logarithmic = ifelse(settings$feature == "length",     T, F)
settings$percentage  = ifelse(settings$feature == "gc_content", T, F)
settings$smallint    = ifelse(settings$feature == "uORF",       T, F)
settings = dplyr::mutate(settings, description=trimws(dplyr::case_when(detail == "exon" ~ "all exons",
                                                                       detail == "UTR5" ~ "5'UTR",
                                                                       detail == "UTR3" ~ "3'UTR",
                                                                       startsWith(detail, "ATG_TRR")  
                                                                       ~ paste("with ATG start codon and TRR stop codon",
                                                                               process_suffix(substring(detail,9))),
                                                                       startsWith(detail, "complete") 
                                                                       ~ paste("with defined both start and stop codons",
                                                                               process_suffix(substring(detail,10))),
                                                                       TRUE ~ gsub('_',' ',detail))))
settings = dplyr::mutate(settings, description=dplyr::case_when(feature == "gc_content" ~ paste("GC content in",
                                                                                                description, "(%)"),
                                                                feature == "length"     ~ paste("Length of",
                                                                                                description, "(bp)"),
                                                                feature == "uORF"       ~ paste("Number of uORFs",
                                                                                                description, "(#)"),
                                                                TRUE                    ~ paste(gsub('_',' ',feature)
                                                                                                , " - ", description)))


#### Process all featurs and print results (correlation coefficients and their p-values)
t(mapply(process_feature, settings$header, settings$description,
         settings$logarithmic, settings$percentage, settings$smallint,
         MoreArgs = list(data=data, prefix=prefix)))
