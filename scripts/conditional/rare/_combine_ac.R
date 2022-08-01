#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){


  stopifnot(file.exists(args$file_cts)) 
  stopifnot(file.exists(args$file_bin)) 
  d1 <- fread(args$file_cts)
  d2 <- fread(args$file_bin)
  d <- merge(d1, d2, by = c("chr","pos","id","ref","alt")) 
  outfile <- paste0(args$out_prefix, ".txt.gz")
  fwrite(d, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--file_cts", default=NULL, required = TRUE, help = "Path to cts phenotypes")
parser$add_argument("--file_bin", default=NULL, required = TRUE, help = "Path to binary phenotypes")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix path")
args <- parser$parse_args()

main(args)









