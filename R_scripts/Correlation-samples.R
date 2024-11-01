#! /usr/bin/env Rscript

################################################################################
## Correlation-samples.R                                                      ##
## ---------------------                                                      ##
## A script to compute correlation coefficients of gene counts between two    ##
## samples and visualize them in a density scatterplot                        ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2022-09-12                                                    ##
## Released under Apache License 2.0                                          ##
################################################################################

# Parse arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("The script requires at least two input arguments - filepaths to replicates")
}
if (length(args) > 2) {
  if (length(args) %% 2 == 1) {
    stop(paste("There is an extra upaired argument:",args[-1]))
  }
  for (i in seq(3, length(args), by=2)) {
    switch(args[i],
      # Where to store output file
      "--output"={ output=args[i+1] },
      # Title of the plot
      "--main"={ label.main=args[i+1] },
      # Description of the x axis
      "--xlab"={ label.1=args[i+1] },
      # Description of the y axis
      "--ylab"={ label.2=args[i+1] },
      # Genes having in both samples less than <min> counts should be ommited
      "--min"={ threshold=args[i+1] },
      # Genes having in both samples less than <min> counts should be ommited
      "--method"={ arg=tolower(args[i+1])
                   if        (arg=="pearson"  || arg=="p") { method="pearson"
                   } else if (arg=="spearman" || arg=="s") { method="spearman"
                   } else if (arg=="kendall"  || arg=="k") { method="kendall"
                   } else { stop(paste("Unrecognized value of --method parameter:",args[i+1])) }}
    )
  }
}

# Set defaults if not overwritten
if (! exists("output")) {
  output = paste(tools::file_path_sans_ext(basename(args[1])), "vs", tools::file_path_sans_ext(basename(args[2])), sep="-")
}
if (! any(endsWith(output, c(".bmp",".eps",".jpeg",".jpg",".pdf",".png",".ps",".svg",".tex",".tiff")))) {
  output = paste(output, "pdf", sep=".")
}
if (! exists("label.1")) {
  label.1 = tools::file_path_sans_ext(basename(args[1]));
}
if (! exists("label.2")) {
  label.2 = tools::file_path_sans_ext(basename(args[2]));
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

# Load libraries
library(ggplot2)
library(ggpointdensity)
library(viridis)

# Read input files and preprocess them; special counters ("no feature", "ambiguous", etc.) starting with '__' are removed (HTSeq-count convention)
process = function(filename, id) {
  replicate = read.csv(filename, sep='\t', header=F)
  colnames(replicate) = c("gene", paste("replicate",id,sep="_"))
  replicate = replicate[!grepl("^__", replicate$gene),]
  replicate[2] = replicate[2] + 1
  return(replicate)
}
replicates = merge(process(args[1],1), process(args[2],2), by="gene")

# Use threshold if defined; '>' instead of '>=' is because there are +1 for each gene count to have log meaningfull even for 0
if (exists("threshold")) {
  replicates = replicates[(replicates$replicate_1>threshold | replicates$replicate_2>threshold),]
}

# Compute and print correlation coefficients
cc = cor(log(replicates$replicate_1), log(replicates$replicate_2), method=method)
print(cc)

# Plot the scatterplot
ggplot(replicates, aes(replicate_1, replicate_2)) +
  geom_pointdensity() +
  scale_color_viridis() +
  labs(colour = "Local density") + 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  coord_fixed() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  ggtitle(label.main) +
  xlab(label.1) +
  ylab(label.2) +
  labs(fill="local density") +
  annotate("text", label=paste0(symbol," == '",sprintf("%0.3f",cc),"'"), parse=T, x=max(replicates$replicate_1), y=min(replicates$replicate_2), hjust="inward", vjust="inward")
ggsave(output)
