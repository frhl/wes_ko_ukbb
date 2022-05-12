#!/usr/bin/env Rscript

library(argparse)


main <- function(args){
    
    a <- as.numeric(as.character(args$a))
    b <- as.numeric(as.character(args$b))
    o <- as.character(args$o)
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    stopifnot(o %in% c("ge","le", "gt", "lt"))

    if (o %in% "ge"){
        result <- as.numeric(a >= b)
    } else if ( o %in% "gt"){
        result <- as.numeric(a > b)
    } else if ( o %in% "le"){
        result <- as.numeric(a <= b)
    } else if ( o %in% "lt"){
        result <- as.numeric(a < b)
    } else {
        stop("operator not yet implemnted")
    }

    write(result, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--a", default=NULL, help = "A value")
parser$add_argument("--b", default=NULL, help = "B value")
parser$add_argument("--o", default=NULL, help = "operator")
args <- parser$parse_args()

main(args)

