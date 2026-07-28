#! /usr/bin/env Rscript

################################################################################
## Correlation-samples.R                                                      ##
## ---------------------                                                      ##
## A script to compute correlation coefficients of gene counts between two    ##
## samples and visualize them in a density scatterplot                        ##
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
  help.common(c("Correlation-samples.R -h                               Prints this message.",
                "Correlation-samples.R <sample_1> <sample_2> [options]  The script computes correlation coefficients of gene",
                "                                                       counts between two samples and visualizes them in",
                "                                                       a density scatterplot.",
                "                                                       The script requires filepaths to two tab-separated tables",
                "                                                       <sample_1> and <sample_2> with readcount for each gene/",
                "                                                       transcript in format [id]{TAB}[count].",
                "",
                "Options:",
                "  --img {pdf|svg|png}                  What file format (and extension) to use to save the plots.",
                "                                       Pdf is used by default.",
                "  --main <title>                       Main title of the plot.",
                "                                       The default title is '[x-axis title] vs. [y-axis title]'.",
                "  --method {Pearson|Spearman|Kendall}  What correlation coefficient should be computed. The value is interpreted",
                "                                       case-insensitive and you may enter only the first letter.",
                "                                       Spearman's rank correlation coefficient is used by default.",
                "  --min <threshold>                    Do not plot genes with less than <threshold> counts in both samples.",
                "                                       No threshold is used by default.",
                "  --output <path>                      Where the ouptut plot should be saved. If <path> is an existing directory,",
                "                                       <path>/[basename of <sample_1>]-vs-[basename of <sample_1>].[file_format]",
                "                                       is used as the filename. If file format extension is missing, it is appended.",
                "                                       If not specified, the current working directory is used as <path>.",
                "  --xlab <title>                       A name of the <sample_1>.",
                "                                       Basename of <sample_1> (without extension) is used by default.",
                "  --ylab <title>                       A name of the <sample_2>.",
                "                                       Basename of <sample_2> (without extension) is used by default"),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=T)
# Parse optional arguments
if (length(args) > 2) {
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
      # Where to store output file
      "--output"={ output=args[i+1] },
      # File format of output histograms
      "--img"={ extension=check.extension(args[i+1]) },
      # Title of the plot
      "--main"={ label.main=args[i+1] },
      # Description of the x axis
      "--xlab"={ label.1=args[i+1] },
      # Description of the y axis
      "--ylab"={ label.2=args[i+1] },
      # Genes having in both samples less than <min> counts should be ommited
      "--min"={ threshold=as.integer(args[i+1]) },
      # Genes having in both samples less than <min> counts should be ommited
      "--method"={ method=parse.correlation.method(args[i+1], args[i])  }
    )
  }
}
# Assign mandatory arguments for easier changes
sample.x = args[1]
sample.y = args[2]


#### External libraries
# Check whether all used libraries are installed
check.installed("ggplot2", "ggpointdensity", "rlang", "tools", "viridis")
# Load intensively used libraries
library(ggplot2)


#### Set defaults if not overwritten
output = construct.filename(output, extension, sample.x, sample.y)$filepath
if (! exists("label.1")) {
  label.1 = tools::file_path_sans_ext(basename(sample.x));
}
if (! exists("label.2")) {
  label.2 = tools::file_path_sans_ext(basename(sample.y));
}
if (! exists("label.main")) {
  label.main = paste(label.1, "vs.", label.2)
}
if (! exists("method")) {
  method = "spearman"
}
switch(method,
  "pearson" ={ symbol="r[P]" },
  "spearman"={ symbol="r[s]" },
  "kendall" ={ symbol="tau" }
)


#### Initial processing of input files
# Read input files and preprocess them; special counters ("no feature", "ambiguous", etc.) starting with '__' are removed (HTSeq-count convention)
process = function(filename, id) {
  replicate = read.csv(filename, sep='\t', header=F)
  colnames(replicate) = c("gene", paste("replicate",id,sep="_"))
  replicate = replicate[!grepl("^__", replicate$gene),]
  replicate[2] = replicate[2] + 1
  return(replicate)
}
replicates = merge(process(sample.x,1), process(sample.y,2), by="gene", all = TRUE)
replicates[is.na(replicates)] = 1
# Use threshold if defined; '>' instead of '>=' is because there are +1 for each gene count to have log meaningfull even for 0
if (exists("threshold")) {
  replicates = replicates[(replicates$replicate_1>threshold | replicates$replicate_2>threshold),]
}


#### Compute and print correlation coefficients
cc = cor(log(replicates$replicate_1), log(replicates$replicate_2), method=method)
print(cc)


#### Plot the scatterplot
ggplot(replicates, aes(replicate_1, replicate_2)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  labs(colour = "Local density") + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(label.main) +
  xlab(label.1) +
  ylab(label.2) +
  labs(fill="local density") +
  annotate("text", label=paste0(symbol," == '",sprintf("%0.3f",cc),"'"), parse=T,
           x=max(replicates$replicate_1), y=min(replicates$replicate_2), hjust="inward", vjust="inward")
ggsave(output)
