#!/usr/bin/env Rscript

library(argparse)
library(data.table)

r_temp = "data/tmp/rtmp"

main <- function(args){
    
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_prefix)))
    stopifnot(is.numeric(as.numeric(args$true_tstat)))
    stopifnot(is.numeric(as.numeric(args$true_p)))

    # read input 
    d <- fread(args$input_path)
    tstat <- as.numeric(d$Tstat)
    pvalue <- as.numeric(d$p.value)
    pvalue_orig <- pvalue
    true_t <- as.numeric(args$true_tstat)
    true_p <- as.numeric(args$true_p)
    
    # Check for conditional
    if ("p.value_cond" %in% colnames(d)){ 
      pvalue <- as.numeric(d$p.value_cond)
    } else {
      pvalue_orig <- rep(NA, length(pvalue_orig))
    }

    # write a matrix of t-statistics and p-values
    out_t <- c(true_t, tstat)
    out_p <- c(true_p, pvalue)
    out_p_orig <- c(true_p, pvalue_orig)
    is_permuted <- c(rep(0, 1), rep(1, length(tstat)))
    dt <- data.table(Tstat=out_t, p=out_p, p_orig=out_p_orig, is_permuted)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dt, outfile, sep = "\t")

    # calculate empirical P-value (note, need to ensure that
    # the correct side is evaluated). 
    count_p_ge_true <- sum(true_p >= pvalue)/length(pvalue) 
    empirical_p <- 1-count_p_ge_true
    write(empirical_p, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--out_prefix", default=NULL, help = "path to the input")
parser$add_argument("--true_tstat", default=NULL, help = "path to the input")
parser$add_argument("--true_p", default=NULL, help = "path to the input")
args <- parser$parse_args()

main(args)

