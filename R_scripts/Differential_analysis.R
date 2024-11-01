#! /usr/bin/env Rscript

################################################################################
## Differential_analysis.R                                                    ##
## ---------------------                                                      ##
## A script for differential analysis of Ribo-Seq data                        ##
##                                                                            ##
## Based on DESeq2 tutorials Beginner's guide to using the DESeq2 package;    ##
## Differential analysis of count data – the DESeq2 package; and Analyzing    ##
## RNA-seq data with DESeq2.                                                  ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2022-09-15                                                    ##
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
  help.common(c("Differential_analysis.R -h                           Prints this message.",
                "Differential_analysis.R <design> <counts> [options]  The script performs a differential expression analysis of",
                "                                                     Ribo-seq data, attach some annotations to genes, and",
                "                                                     generates some common plots.",
                "                                                     The script requires a filepath to a tab-separated table <design>",
                "                                                     with design of the experiment (containing rownames - sample",
                "                                                     identifiers - and columns 'assay' and 'siRNA'; with colnames)",
                "                                                     and a tab-separated table <counts> with gene counts (containing",
                "                                                     rownames - Ensembl gene identifiers - and one column for each",
                "                                                     sample identifier).",
                "",
                "Options:",
                "  --img {pdf|svg|png}          What file format (and extension) to use to save the plots.",
                "                               Pdf is used by default.",
                "  --output <dirname>           Directory where to save result. If not specified, <counts> filename without",
                "                               its extension is taken and appends the lowest number such that a new directory",
                "                               will be created.",
                "  --padj <threshold>           Space-separated thresholds on p_adj to be considered as significantly changed.",
                "                               The default threshold is 0.05.",
                "  --pca_ids {TRUE|FALSE}       Whether to show sample identifiers in PCA plots. They are not shown by default.",
                "  --use_bm_cache {TRUE|FALSE}  Whether to use cache when querying biomaRt. Cache is not use by default (it",
                "                               makes a problem in some older versions of R."),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=T)
# Parse optional arguments
for (i in seq(3, length(args), by=2)) {
  switch(args[i],
    # Directory where to store results (will be created)
    "--output"={ workingDir=args[i+1] },
    # Plots file format (supported are pdf, svg and png
    "--img"={ extension=args[i+1] },
    # Threshold on p_adjusted value to consider a gene as significantly differentially expressed
    "--padj"={ padjs = as.numeric(unlist(strsplit(args[i+1], split=" "))) },
    # Whether PCA plots should contain sample identifiers
    "--pca_ids"={ pcaIds = parse.boolean(args[i+1], "pca_ids") },
    # Whether BioMart cache should be used (it does not works in some combinations of installed packages)
    "--use_bm_cache"={ useBMcache = parse.boolean(args[i+1], "use_bm_cache") }
  )
}
# Assign mandatory arguments for easier changes
designPath = args[1]
countsPath = args[2]


#### Set defaults if not overwritten
if (! exists("workingDir")) {
  # Derive dirname from counts name by adding the lowest unsigned integer to get non-existing dirname
  baseName = tools::file_path_sans_ext(countsPath)
  i = 1
  repeat {
    if (!file.exists(paste(baseName, i, sep="-"))) {
      workingDir = paste(baseName, i, sep="-")
      break
    }
    i = i+1
  }
}
if (!exists("extension")) {
  extension = "pdf"
} else {
  check.extension(extension)
}
if (! exists("padjs")) {
  padjs = c(0.05)
}
if (!exists("pcaIds")) {
  pcaIds = FALSE
}
if (!exists("useBMcache")) {
  useBMcache = FALSE
}


#### External libraries
# Check whether all used libraries are installed
check.installed("biomaRt", "DESeq2", "dplyr", "ggplot2", "gplots", "ggrepel", "RColorBrewer", "tools")
# Load intensively used libraries
library(DESeq2)
library(ggplot2)


#### Initial processing of inputs
# Read input files, select only samples required by designTable and make samples in the same order (required by DESeq2)
designTable = read.delim(designPath, sep="\t", header=TRUE, row.names=1)
countsTable = read.delim(countsPath, sep="\t", header=TRUE, row.names=1)[rownames(designTable)]
# Set up dataframe for DESeq2
data = DESeqDataSetFromMatrix(countData=countsTable, colData=designTable, design=~assay*siRNA)
# Collapse samples that should be collapsed
if (dim(designTable)[1] != dim(unique(dplyr::group_by_all(designTable)))[1]) {
  data = collapseReplicates(data, factor(paste(designTable$assay, designTable$siRNA, designTable$id, sep="_")), colnames(data))
}
# Set reference levels to mRNA and nt
data$assay = relevel(data$assay, "mRNA")
data$siRNA = relevel(data$siRNA, "nt")
# Remove rows with zero counts
data = DESeq(data[rowSums(counts(data)) >0])
# normalize to library size
data.transformed = rlog(data)


#### Create working dir and move into it
dir.create(workingDir, recursive=T)
setwd(workingDir)


#### Basic plots section
# Plot dispersion estimates
get(extension)(paste("dispersion_estimate",extension, sep="."))
plotDispEsts(data, ylim=c(1e-5, 1e1))
dev.off()
# Visualize sample distances
get(extension)(paste("sample_distance-heatmap",extension, sep="."), width=10,height=10)
gplots::heatmap.2(as.matrix(dist(t(assay(data.transformed)))), trace="none", col=colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255))
dev.off()
# Make PCA plots
plot_PCA = function(data, filename) {
  pca = plotPCA(data, intgroup=c("siRNA","assay"), returnData=TRUE)
  percentVar = round(100*attr(pca, "percentVar"))
  plot = ggplot(pca, aes(PC1, PC2, fill=siRNA, shape=assay)) +
         geom_point(size=3) +
         guides(fill=guide_legend(override.aes=list(shape=22))) +
         scale_shape_manual(values=c(21,24)) +
         xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) +
         coord_fixed(ratio=1)
  if (pcaIds) {
    plot = plot + ggrepel::geom_text_repel(aes(label=name), size=2, max.overlaps = Inf)
  }
  ggsave(paste(filename, extension, sep="."), plot)
}
plot_PCA(data.transformed, "PCA")
for (assay in c("FP", "mRNA")) {
  plot_PCA(data.transformed[, data.transformed$assay %in% c(assay)], paste("PCA",assay, sep="-"))
}


#### Create results tables and show summaries
ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
process_results = function(result, type) {
  # Show summaries
  lapply(padjs,summary,object=result)
  # Add other gene IDs and GO IDs
  result$ensembl = rownames(result)
  mapping = biomaRt::getBM(attributes=c("ensembl_gene_id","entrezgene_id","hgnc_symbol","go_id","name_1006"), filters="ensembl_gene_id", values=result$ensembl, mart=ensembl, useCache=useBMcache)
  mapping = mapping[match(result$ensembl, mapping$ensembl_gene_id),]
  result$entrez      = mapping$entrezgene_id
  result$hgnc_symbol = mapping$hgnc_symbol
  result$go_id       = mapping$go_id
  result$name_1006   = mapping$name_1006
  result.df = as.data.frame(result)
  # Save result tables
  write.table(result.df, sep="\t", col.names=NA, file=paste0(type,"-",label,"_vs_nt.tbl"))
  for (padj in padjs) {
    sig = result.df[which(result.df$padj < padj),]
    write.table(sig[order(sig$log2FoldChange),], sep="\t", col.names=NA, file=paste0("sig",padj,"-",type,"-",label,"_vs_nt.tbl"))
  }
  # Create MA plot
  get(extension)(paste0("MA_plot-", type, "-", label, "_vs_nt.", extension))
  plotMA(result, ylim=c(-2,2))
  dev.off()
  # Generate histogram of pvalue distribution
  get(extension)(paste0("pvalue-histogram-", type, "_", label, "_vs_nt.", extension))
  hist(result$pvalue, breaks=20, col="grey", main=paste("Histogram of p-values for",label,"vs nt"), xlab=paste("p-value",type,label,"vs nt"))
  dev.off()
  # Create Volcano plot
  result.df$Expression = dplyr::case_when(result.df$log2FoldChange > 0 & result.df$padj <= 0.05 ~ "Upregulated",
                                          result.df$log2FoldChange < 0 & result.df$padj <= 0.05 ~ "Downregulated",
                                          TRUE ~ "Unchanged")
  ggplot(result.df, aes(x=log2FoldChange, y=-log10(padj), col=Expression)) +
    geom_point(size = 0.5) +
    xlab(expression("log"[2]*"(FC)")) + 
    ylab(expression("-log"[10]*"(p-adjusted)")) +
    scale_color_manual(values=c("blue", "grey", "red")) +
    guides(colour = guide_legend(override.aes=list(size=2), reverse=T)) 
  ggsave(paste0("Volcano-", type, "-", label, "_vs_nt.", extension))
}
# For each non-nt siRNA process mRNA, FP and TE
labels = unique(designTable[,"siRNA"])
for (label in labels[labels!="nt"]) {
  process_results(results(data, name=paste("siRNA",label,"vs_nt", sep="_")), "mRNA")
  process_results(results(data, contrast=list(c(paste("siRNA",label,"vs_nt", sep="_"), paste0("assayFP.siRNA",label)))), "FP")
  process_results(results(data, name=paste0("assayFP.siRNA",label)), "TE")
}
