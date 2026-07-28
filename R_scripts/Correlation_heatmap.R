#! /usr/bin/env Rscript

################################################################################
## Correlation-heatmap.R                                                      ##
## ---------------------                                                      ##
## Generate heatmap of samples based on their similarity                      ##
##                                                                            ##
## Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)                         ##
## Last update: 2023-09-01                                                    ##
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
  help.common(c("Correlation_heatmap.R -h                 Prints this message.",
                "Correlation_heatmap.R <input> [options]  The script transforms a similarity table from long to wide format",
                "                                         and generates a heatmap visualizing the similarity of samples.",
                "",
                "Options:",
                "  --heatmap <path>     Where the ouptut heatmap should be saved. If <path> is an existing directory,",
                "                       <path>/[basename of <input>].[file_format] is used as the filename. If file",
                "                       format extension is missing, it is appended. If not specified, the current",
                "                       working directory is used as <path>.",
                "  --img {pdf|svg|png}  What file format (and extension) to use to save the heatmap. Pdf is used by default.",
                "  --table <path>       Where the ouptut table in wide format should be saved. If <path> is an existing",
                "                       directory, <path>/[basename of <input>].tbl is used as the filename. If file",
                "                       format extension is missing, 'it'.tbl' is appended. If not specified, the current",
                "                       working directory is used as <path>."),
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
      # Output files format
      "--img"={ extension=check.extension(args[i+1]) },
      # Path to the output heatmap
      "--heatmap"={ heatmap=args[i+1] },
      # Path to the output table in wide format
      "--table"={ widetable=args[i+1] }
    )
  }
}
# Assign mandatory arguments for easier changes
input = args[1]


#### External libraries
# Check whether all used libraries are installed
check.installed("gplots", "reshape2", "tools")
# Load intensively used libraries


#### Set defaults if not overwritten
heatmap = construct.filename(heatmap, extension, input)
if (exists("widetable", mode="character")) {
  if (nchar(tools::file_ext(widetable)) == 0) {
    widetable = paste0(widetable, ".tbl")
  }
} else {
  widetable = paste0(tools::file_path_sans_ext(input), ".tbl")
}


#### Generate outputs and save them in output files
# Read the table and transform it to wide format
coefs = reshape2::acast(read.delim(input,sep="\t",header=F), V1~V2, value.var="V3")
# Write the table
write.table(coefs, widetable, col.names=NA)
# Generate and save the heatmap
get(heatmap$extension)(heatmap$filepath)
gplots::heatmap.2(coefs, trace="none", cexRow=0.5, cexCol=0.5)
dev.off()
