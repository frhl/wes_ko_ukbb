#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){
    
    input_path1 <- args$input_path1
    input_path2 <- args$input_path2
    operation <- args$operation
    outfile <- args$outfile
    
    d1 <- fread(input_path1, header = FALSE)
    d2 <- fread(input_path2, header = FALSE)
    if (operation == "intersect") {
        d <- data.table(intersect(d1$V1, d2$V1))
    } else {
        stop(paste(operation, "is not allowed"))   
    }

    write(paste("writing to", outfile), stdout())   
    fwrite(d, outfile, col.names = FALSE, quote = FALSE) 
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path1", default=NULL, help = "?")
parser$add_argument("--input_path2", default=NULL, help = "?")
parser$add_argument("--operation", default="intersect", help = "intersect or union")
parser$add_argument("--outfile", default=NULL, help = "")
args <- parser$parse_args()

main(args)

