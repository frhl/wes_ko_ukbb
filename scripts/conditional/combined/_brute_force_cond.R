#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

  common_exists <- (file.exists(args$file_common_markers))
  rare_exists <- (file.exists(args$file_rare_markers))
    
  if (common_exists & rare_exists) {
    d1 <- fread(args$file_common_markers, header = FALSE)
    d2 <- fread(args$file_rare_markers, header = FALSE)
    d <- rbind(d1, d2)
  } else if (common_exists) {
    d <- fread(args$file_common_markers, header = FALSE)
  } else if (rare_exists) {
    d <- fread(args$file_rare_markers, header = FALSE)
  } else {
    stop(paste("Expected at least", args$file_rare_marker, "or", args$file_common_markers, "to exist!"))
  }

  # rename and order markers
  markers <- gwastools::order_markers(d$V1)
  fwrite(data.table(x=markers), args$outfile, row.names = FALSE, col.names = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--file_common_markers", default=NULL, required = TRUE, help = "Path to common markers")
parser$add_argument("--file_rare_markers", default=NULL, required = TRUE, help = "Path to rare markers")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









