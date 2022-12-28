#!/usr/bin/env Rscript

library(argparse)
library(data.table)


is_overlapping <- function(x1,x2,y1,y2){
    return(max(x1,y1) <= min(x2,y2))
    }

# get overlapping intervals
which_overlap <- function(x1, x2, y1, y2){
    nx <- length(x1)
    ny <- length(y1)
    index <- lapply(1:nx, function(i){
        index <- which(unlist(lapply(1:ny, function(j){
            x1i <- x1[i]
            x2i <- x2[i]
            y1j <- y1[j]
            y2j <- y2[j]
            stopifnot(x1i <= x2i)
            stopifnot(y1j <= y2j)
            overlap <- is_overlapping(x1i, x2i, y1j, y2j)
        })))
        return(unlist(index))
    })
    # ensure that missing indexes are NAs
    index[!(index %in% 1:ny)] <- NA
    return(unlist(index))
 }


main <- function(args){



    # (2) append with permuted data
    fwrite(final, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, help = "chromosome")
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--permutations", default=NULL, help = "number of times the gene should be permuted")
parser$add_argument("--seed", default=NULL, help = "seed for randomizer")
parser$add_argument("--vcf_id", default="GENE", help = "Substitute for rsid")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

