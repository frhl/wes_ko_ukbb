#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# helper, read a hail list
clean_hail_list <- function(x, comma_sub = ', ') {
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}


# Read with and annotate basename
fread_with_basename <- function(f){
    d <- fread(f)
    file_sans_gz <- tools::file_path_sans_ext(basename(f))
    file_sans_txt <- tools::file_path_sans_ext(file_sans_gz)
    d$basename <- file_sans_txt
    return(d)
}

# generate a table of genes in file
tabulate_MarkerID_basename <- function(M){
    stopifnot(all(c("CHR","MarkerID",'basename') %in% colnames(M)))
    # enumerate all MarkerID basenames
    d <- as.data.frame(table(M$MarkerID, M$basename))
    colnames(d) <- c("MarkerID", "basename", "Freq")
    d <- d[d$Freq > 0,]
    d$Freq <- NULL
    # merge with chromosome 
    chroms <- M[,c("CHR","MarkerID")]
    chroms <- chroms[!duplicated(chroms),]
    d <- merge(d, chroms, all.x = TRUE)
    return(d)
}

main <- function(args){

    print(args)
    stopifnot(dir.exists(args$spa_cts_dir))     
    stopifnot(dir.exists(args$spa_bin_dir))

    spa_cts_files <- list.files(args$spa_cts_dir, pattern = ".txt.gz", full.names = TRUE)
    spa_bin_files <- list.files(args$spa_bin_dir, pattern = ".txt.gz", full.names = TRUE)
    
    # read in all the SPA step 2 results
    spa_cts_full <- do.call(rbind, lapply(spa_cts_files, fread_with_basename))
    spa_bin_full <- do.call(rbind, lapply(spa_bin_files, fread_with_basename))

    # write out each gene involved in each comparison (avg of ~1700 genes per SPA)
    spa_full <- rbind(spa_cts_full[,c("CHR","MarkerID","basename", "p.value")], spa_bin_full[,c("CHR","MarkerID","basename", "p.value")])
    out_prefix_true_p <- paste0(args$out_prefix, "_true_p.tsv.gz")
    fwrite(spa_full, out_prefix_true_p, sep = '\t')
    out_prefix_genes <- paste0(args$out_prefix, "_genes.tsv.gz")
    fwrite(tabulate_MarkerID_basename(spa_full), out_prefix_genes, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_cts_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--spa_bin_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--tsv_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--p_cutoff", default=10e-12, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

