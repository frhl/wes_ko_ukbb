#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)

main <- function(args){

    x <- args$markers
    stopifnot(x != "")
    markers <- unlist(strsplit(x, split = ','))
    ordered_markers <- order_markers(markers)
    out <- paste0(ordered_markers, collapse = ',') 
    write(out, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--markers", default=NULL, help = "list of comma sperated markers")
args <- parser$parse_args()

main(args)

