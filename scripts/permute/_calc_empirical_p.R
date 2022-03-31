#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){
    
    stopifnot(file.exists(args$input_path))
    d <- fread(args$input_path)
    tstat <- as.numeric(d$Tstat)
    true_t <- as.numeric(args$true_tstat)

    p_value <- sum(tstat >= true_t)/length(tstat)
    write(p_value, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--true_tstat", default=NULL, help = "path to the input")
args <- parser$parse_args()

main(args)

