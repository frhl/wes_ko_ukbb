#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){
 
    stopifnot(file.exists(args$input_path))

    # read input
    d <- fread(args$input_path)

    # We include the "real" marker in the VCF to extract
    # the reference/true P-value and T-statistic
    true_p_marker_name <- "actual"
    stopifnot(true_p_marker_name %in% d$MarkerID)
    d_actual <- d[d$MarkerID == true_p_marker_name,]
    if ("p.value_c" %in% colnames(d_actual)){
      true_p <- as.numeric(d_actual$p.value_c)
      true_t <- as.numeric(d_actual$Tstat_c)
    } else {
      true_p <- as.numeric(d_actual$p.value)
      true_t <- as.numeric(d_actual$Tstat)
    }
    write(true_p, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "Merged SAIGE file with 'actual' (true) marker included.")
args <- parser$parse_args()

main(args)

