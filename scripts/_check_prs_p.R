#!/usr/bin/env Rscript

library(argparse)

main <- function(args){

  stopifnot(file.exists(args$ldsc))
  ldsc <- readRDS(args$ldsc)
  pvalue <- ldsc$coefficients$pvalue[2]
  out <- ifelse(pvalue < as.numeric(args$p_cutoff), 1 , 0)
  write(out, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ldsc", default=NULL, help = "?")
parser$add_argument("--p_cutoff", default=0.05, help = "?")
args <- parser$parse_args()

main(args)

