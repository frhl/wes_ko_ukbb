#!/usr/bin/env Rscript


devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

  print(args)

  stopifnot(file.exists(args$file_common_markers))
  stopifnot(file.exists(args$file_rare_markers))
    
  # common markers may be empty, but there will always be rare markers
  d1 <- fread(args$file_common_markers, header = FALSE)
  d2 <- fread(args$file_rare_markers, header = FALSE)



  # ensure variants follow input formatting required by saige
  #markers <- gwastools::order_markers(markers)
  #fwrite(data.table(x=markers), args$outfile, row.names = FALSE, col.names = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--file_common_markers", default=NULL, required = TRUE, help = "Path to common markers")
parser$add_argument("--file_rare_markers", default=NULL, required = TRUE, help = "Path to rare markers")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









