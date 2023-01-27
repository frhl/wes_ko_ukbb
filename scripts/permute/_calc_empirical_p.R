#!/usr/bin/env Rscript

library(argparse)
library(data.table)

r_temp = "data/tmp/rtmp"

main <- function(args){
    
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    stopifnot(is.numeric(as.numeric(args$true_tstat)))
    stopifnot(is.numeric(as.numeric(args$true_p)))
    
    # get true values from non-permuted data
    true_t <- as.numeric(args$true_tstat)
    true_p <- as.numeric(args$true_p)
    true_AC <- as.numeric(args$true_AC)
    true_marker <- as.character(args$true_marker)

    # read input
    d <- fread(args$input_path)
   
    # exclude real markers (non-permuted stuff)
    stopifnot("MarkerID" %in% colnames(d))
    if (args$exclude_real_markers) {
        bool_real_marker <- !grepl("ENSG", d$MarkerID)
        n_real_marker <- sum(bool_real_marker)
        if (n_real_marker > 0) {
            write(paste("Excluded", n_real_marker, "real marker(s)."), stderr()) 
            d <- d[!bool_real_marker, ]
        }
    }

    # get permuted values
    marker <- as.character(d$MarkerID)
    tstat <- as.numeric(d$Tstat)
    pvalue <- as.numeric(d$p.value)
    pvalue_orig <- pvalue
    AC <- as.numeric(d$AC_Allele2)

    # Check for conditional
    if ("p.value_c" %in% colnames(d)){ 
      pvalue <- as.numeric(d$p.value_c)
    } else {
      pvalue_orig <- rep(NA, length(pvalue_orig))
    }

       # write a matrix of t-statistics and p-values
    out_t <- c(true_t, tstat)
    out_p <- c(true_p, pvalue)
    out_p_orig <- c(true_p, pvalue_orig)
    is_permuted <- c(rep(0, 1), rep(1, length(tstat)))
    out_marker <- c(true_marker, marker)
    out_ac <- c(true_AC, AC)
    dt <- data.table(Tstat=out_t, p=out_p, p_orig=out_p_orig, out_marker, out_ac, is_permuted)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dt, outfile, sep = "\t")
    
    # count how many times the original P-value could be found
    idx_same_p <- which(true_p %in% dt$p[dt$is_permuted == 1])  
    n_same_p <- length(idx_same_p)
    idx_string <- paste0(idx_same_p, collapse = ",")
    #write(paste("\nThe permuted P-value is equal to the true P-value", n_same_p, "times (idx=",idx_string,"). Saving to", outfile), stdout())

    # calculate empirical P-value (note, need to ensure that
    # the correct side is evaluated). 
    count_p_ge_true <- sum(true_p >= pvalue)/length(pvalue) 
    empirical_p <- count_p_ge_true
    write(empirical_p, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "path to saige (merge) file of markers")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for output")
parser$add_argument("--true_tstat", default=NULL, help = "the true t-statistic from non-permuted analysis")
parser$add_argument("--true_p", default=NULL, help = "the true p-value from non-permuted analysis")
parser$add_argument("--true_AC", default=NULL, help = "the true allele count from non-permuted analysis")
parser$add_argument("--true_marker", default=NULL, help = "the true marker name")
parser$add_argument("--exclude_real_markers", default=FALSE, action="store_true", help = "Exclude real conditioning markers (from conditioning analysis)")
args <- parser$parse_args()

main(args)

