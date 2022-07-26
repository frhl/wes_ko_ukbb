#!/usr/bin/env Rscript

library(argparse)

main <- function(args){

    x <- args$markers
    stopifnot(x != "")
   
    # combine markers and order
    markers <- unlist(strsplit(x, split = ','))
    d <- data.frame(do.call(rbind, strsplit(markers, split = ':')))
    colnames(d) <- c("chr","pos","a1", "a2")
    d$pos <- as.numeric(d$pos)
    d$chr_int <- as.numeric(gsub("chr","",d$chr))
    d <- d[order(d$chr_int, d$pos),]
    d$chr_int <- NULL

    # generate outout and append
    out <- apply( d , 1 , paste , collapse = ":")
    out <- gsub(" ", "", out)
    out <- paste0(out, collapse = ',') 

    write(out, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--markers", default=NULL, help = "list of comma sperated markers")
args <- parser$parse_args()

main(args)

