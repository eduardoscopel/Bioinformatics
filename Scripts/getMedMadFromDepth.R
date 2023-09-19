#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).depth", call.=FALSE)
}
depthtab <- fread(args[1], col.names = c("Chromosome", "Position", "Depth"))
cat(sub(".depth","",args[1]), median(depthtab$Depth, mad(depthtab$Depth), '\n',sep = '\t')
