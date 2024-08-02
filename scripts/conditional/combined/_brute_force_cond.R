#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

  common_exists <- (file.exists(args$file_common_markers))
  rare_exists <- (file.exists(args$file_rare_markers))
  collapsed_exists <- (file.exists(args$file_collapsed_markers))
  exit_without_writing <- FALSE  
  
  if (common_exists & rare_exists & collapsed_exists) {
    write(paste("Using common/rare/collapsed for", args$outfile), stderr())
    d1 <- fread(args$file_common_markers, header = FALSE)
    d2 <- fread(args$file_rare_markers, header = FALSE)
    d3 <- fread(args$file_collapsed_markers, header = FALSE)
    d <- rbind(d1, d2, d3)
  } else if (rare_exists & collapsed_exists) {
    write(paste("Using rare/collapsed for", args$outfile), stderr())
    d1 <- fread(args$file_collapsed_markers, header = FALSE)
    d2 <- fread(args$file_rare_markers, header = FALSE)
    d <- rbind(d1, d2)
  } else if (common_exists & collapsed_exists) {
    write(paste("Using common/collapsed for", args$outfile), stderr())
    d1 <- fread(args$file_common_markers, header = FALSE)
    d2 <- fread(args$file_collapsed_markers, header = FALSE)
    d <- rbind(d1, d2)
  } else if (common_exists & rare_exists) {
    write(paste("Using common/rare for", args$outfile), stderr())
    d1 <- fread(args$file_common_markers, header = FALSE)
    d2 <- fread(args$file_rare_markers, header = FALSE)
    d <- rbind(d1, d2)
  } else if (collapsed_exists) {
    write(paste("Using collapsed for", args$outfile), stderr())
    d <- fread(args$file_collapsed_markers, header = FALSE)
  } else if (common_exists) {
    d <- fread(args$file_common_markers, header = FALSE)
  } else if (rare_exists) {
    write(paste("Using are for", args$outfile), stderr())
    d <- fread(args$file_rare_markers, header = FALSE)
  } else {
    write(paste("Expected at least", args$file_rare_marker, "or", args$file_common_markers, "to exist. No marker file created."), stderr())
    exit_without_writing <- TRUE  
  }

  # rename and order markers
  if (exit_without_writing == FALSE) {
    markers <- gwastools::order_markers(d$V1)
    fwrite(data.table(x=markers), args$outfile, row.names = FALSE, col.names = FALSE)
  }
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--file_common_markers", default=NULL, required = TRUE, help = "Path to common markers")
parser$add_argument("--file_rare_markers", default=NULL, required = TRUE, help = "Path to rare markers")
parser$add_argument("--file_collapsed_markers", default=NULL, required = TRUE, help = "Path to rare markers")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









