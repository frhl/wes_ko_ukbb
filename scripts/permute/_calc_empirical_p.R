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
    true_t <- as.numeric(args$true_tstat)
    true_p <- as.numeric(args$true_p)

    # write a matrix of t-statistics and p-values
    out_t <- c(true_t, tstat)
    out_p <- c(true_p, pvalue)
    is_permuted <- c(rep(0, 1), rep(1, length(tstat)))
    dt <- data.table(Tstat=out_t, p=out_p, is_permuted)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(dt, outfile, sep = "\t")

    # calculate empirical P-value
    empirical_p_gt <- sum(true_p >= pvalue)/length(pvalue) 
    empirical_p_ge <- sum(true_p > pvalue)/length(pvalue) 
    empirical_p <- (empirical_p_gt + empirical_p_ge) / 2
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

