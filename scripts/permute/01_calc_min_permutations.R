#!/usr/bin/env Rscript

library(argparse)
library(data.table)



# Read with and annotate basename
fread_with_basename <- function(f){
    d <- fread(f)
    if (! "V1" %in% colnames(d)){
        file_sans_gz <- tools::file_path_sans_ext(basename(f))
        file_sans_txt <- tools::file_path_sans_ext(file_sans_gz)
        d$basename <- file_sans_txt
        return(d)
    } else {
        warning(paste(f, "had NA's in colname! Skipping.."))
        return(NULL)
    }
}

main <- function(args){

    print(args)
    stopifnot(dir.exists(args$spa_cts_dir))
    stopifnot(dir.exists(args$spa_bin_dir))

    spa_cts_files <- list.files(args$spa_cts_dir, pattern = ".txt.gz", full.names = TRUE)
    spa_bin_files <- list.files(args$spa_bin_dir, pattern = ".txt.gz", full.names = TRUE)

    # read in all the SPA step 2 results
    spa_cts_full <- rbindlist(lapply(spa_cts_files, fread_with_basename), fill = TRUE)
    spa_bin_full <- rbindlist(lapply(spa_bin_files, fread_with_basename), fill = TRUE)
    write("Loaded all cts/binary files..", stdout())

    print(head(spa_cts_full, n=2))
    print(head(spa_bin_full, n=2))

    # combine and subset
    keep <- c("CHR","MarkerID","basename", "p.value", "Tstat", "p.value_c", "Tstat_c")
    spa_full <- rbindlist(list(spa_cts_full, spa_bin_full))
    spa_full <- spa_full[, colnames(spa_full) %in% keep, with = FALSE]
    spa_full <- spa_full[,..keep]
    
    n <- nrow(spa_full) 

    if (args$use_cond_p) {
        
       # use p-value from conditional values if applicable 
       bool_c <- !is.na(spa_full$p.value_c)
       spa_full$pvalue <- spa_full$p.value
       spa_full$tstat <- spa_full$Tstat
       spa_full$cond <- bool_c
       spa_full$pvalue[bool_c] <- spa_full$p.value_c[bool_c] 
       spa_full$tstat[bool_c] <- spa_full$Tstat_c[bool_c]
    }    

    # format to avoid scientific notation 
    out_prefix_true_p_detailed <- paste0(args$out_prefix, "_true_p_detailed.tsv.gz")
    fwrite(spa_full, out_prefix_true_p_detailed, sep = '\t', quote = FALSE, na = NA)

    # remove old columns
    spa_full$p.value_c <- NULL
    spa_full$Tstat_c <- NULL
    spa_full$cond <- NULL
 
    # the filtered version
    spa_full$pvalue <- format(spa_full$pvalue, scientific = FALSE)
    spa_full$tstat <- format(spa_full$tstat, scientific = FALSE)
    out_prefix_true_p <- paste0(args$out_prefix, "_true_p.tsv.gz")
    fwrite(spa_full, out_prefix_true_p, sep = '\t', quote = FALSE)
    write(paste("Note: wrote", out_prefix_true_p),stderr())
    
    # write out genes and chromosomes used 
    spa_genes <- spa_full[,c("MarkerID","CHR")]
    spa_genes <- spa_genes[!duplicated(spa_genes),]
    out_prefix_genes <- paste0(args$out_prefix, "_genes.tsv.gz")
    fwrite(spa_genes, out_prefix_genes, sep = '\t')
    write(paste("Note: wrote", out_prefix_genes),stderr())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_cts_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--spa_bin_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--tsv_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--use_cond_p", default=FALSE, action = "store_true", help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

