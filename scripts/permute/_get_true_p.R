#!/usr/bin/env Rscript

#python_bin <- findpython::find_python_cmd()
#write(paste("Python path:", python_bin), stderr())
#options('python_cmd'='/well/lindgren/users/mmq446/conda/skylake/envs/rpy/bin/python')

#library(argparse)
library(data.table)

main <- function(args){
 
    stopifnot(file.exists(args$input_path))
    stopifnot(args$what %in% c("p", "t"))
        
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

    # there is one true marker per VCF, we ensure
    # they are all the same and take the unique
    true_p <- unique(true_p)
    true_t <- unique(true_t)
    stopifnot(length(true_p) == 1)
    stopifnot(length(true_t) == 1)
    if (is.na(true_p)) stop(paste("true P/t-statistic is NA for", args$input_path))
    
    if (args$what == "p"){
        write(true_p, stdout())
    } else if (args$what == "t"){
        write(true_t, stdout())
    } else {
        stop("arg 'what' should be either 'p' (p-value) or 't' (t-statistic)")
    }

}

# add arguments
#parser <- ArgumentParser()
#parser$add_argument("--input_path", default=NULL, help = "Merged SAIGE file with 'actual' (true) marker included.")
#args <- parser$parse_args()

my_args <- commandArgs(trailingOnly=TRUE)
#write(paste("args:", my_args), stderr())
args <- list(input_path=my_args[1], what=my_args[2])
main(args)

