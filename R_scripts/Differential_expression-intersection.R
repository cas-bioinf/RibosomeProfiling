#! /usr/bin/env Rscript

################################################################################
## Differential_analysis-intersection.R                                       ##
## ---------------------                                                      ##
## A script to visualize intersection gene sets using Euler diagram. It is   ##
## expected that gene ids are in the first column of tab-separated files.     ##
## Note that at most 5 gene sets can be visualized.                           ##
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
  help.common(c("Differential_expression-intersection.R -h                          Prints this message.",
                "Differential_expression-intersection.R <output> <input>+ <label>+  The script visualizes intersection of gene",
                "                                                                   sets using Euler diagram.",
                "                                                                   The diagram is saved in  <output> file. If",
                "                                                                   <output> does not contains an extension,",
                "                                                                   PDF fileformat is used and '.pdf' extension",
                "                                                                   is attached to the file name.",
                "                                                                   <input>+ files should be in tab-separated",
                "                                                                   fileformat and only their first columns",
                "                                                                   are used.",
                "                                                                   <label>+ are used as labels for the sets",
                "                                                                   where i-th <label> should correspond to",
                "                                                                   i-th <output> file."),
              error)
}


#### Parse command line arguments
args = get.arguments()
# Check basic validity of command line arguments
check.parity(args, even=F)


#### External libraries
# Check whether all used libraries are installed
check.installed("futile.logger", "ggplot2", "RColorBrewer", "tools", "VennDiagram")


#### Process mandatory arguments
# Save <output> path, will be validated later
filename = args[1]
args = skip(args, 1)
# Read input files
tables = lapply(args[1:(length(args)/2)], function(x) rownames(read.csv(x, sep='\t', header=T, row.names=1)))
args = skip(args, length(args)/2)
# Generate sample names
names = mapply(function(label, data) paste0(label, " (", length(data), ")"), args, tables)
# Validation of <output> path
filename = construct.filename(filename, NA, args)$filepath


#### Plot the diagram
# Choose colors for for the diagram; "[1:length(tables)]" is for the case less than 3 samples are visualized
lines = RColorBrewer::brewer.pal(length(tables), "Set1")[1:length(tables)]
area = ggplot2::alpha(lines, 0.3)
# Generate the diagram (and suppress its log)
futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger")
diagram = VennDiagram::venn.diagram(x=tables, category.names=names, col=lines, fill=area, filename=NULL)
ggplot2::ggsave(filename, diagram)
