#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

  print(args)
  stopifnot(dir.exists(args$input_dir))



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

