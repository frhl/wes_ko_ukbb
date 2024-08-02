#!/usr/bin/env Rscript

# read merged saige file and get the ranked top x p.value

library(argparse)
library(data.table)


main <- function(args){
    
    if (!file.exists(args$input_path)) stop(paste(args$input_path, "does not exist!"))
    d <- fread(args$input_path)
    d <- d[grepl(d$MarkerID, pattern="ENSG"),]
    # allow conditional P-values
    if ("p.value_c" %in% colnames(d)){
      p <- as.numeric(d$p.value_c)
    } else {
      p <- as.numeric(d$p.value)
    } 
    sorted_p <- sort(p)
    top_p <- sorted_p[as.numeric(args$select_min_p)]
    top_p <- top_p * as.numeric(args$multiply_p)
    write(top_p, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--conditional", default=NULL, help = "path to the input")
parser$add_argument("--select_min_p", default=1, help = "path to the input")
parser$add_argument("--multiply_p", default=1, help = "path to the input")
args <- parser$parse_args()

main(args)

