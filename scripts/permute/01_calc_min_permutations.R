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
    spa_cts_full <- do.call(rbind, lapply(spa_cts_files, fread_with_basename))
    spa_bin_full <- do.call(rbind, lapply(spa_bin_files, fread_with_basename))
    write("Loaded all cts/binary files..", stdout())

    print(head(spa_cts_full, n=2))
    print(head(spa_bin_full, n=2))

    # write out gene-spa p-value pairs
    spa_full <- rbind(
        spa_cts_full[,c("CHR","MarkerID","basename", "p.value", "Tstat")], 
        spa_bin_full[,c("CHR","MarkerID","basename", "p.value", "Tstat")]
    )
    
    # format to avoid scientific notation 
    spa_full$p.value <- format(spa_full$p.value, scientific = FALSE)
    spa_full$Tstat <- format(spa_full$Tstat, scientific = FALSE)
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
args <- parser$parse_args()

main(args)

