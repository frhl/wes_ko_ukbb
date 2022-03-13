#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# helper
clean_hail_list <- function(x, comma_sub = ', ') {
    x <- gsub('(\\[)|(\\])|(\\")','',x)
    x <- gsub('\\,',comma_sub,x)
    return(x)
}

main <- function(args){

    print(args)
    stopifnot(dir.exists(args$spa_cts_dir))     
    stopifnot(dir.exists(args$spa_bin_dir))

    spa_cts_files <- list.files(args$spa_cts_dir, pattern = ".txt.gz", full.names = TRUE)
    spa_bin_files <- list.files(args$spa_bin_dir, pattern = ".txt.gz", full.names = TRUE)
    
    spa_cts <- do.call(rbind, lapply(spa_cts_files, function(f){fread(f)}))
    spa_bin <- do.call(rbind, lapply(spa_bin_files, function(f){fread(f)}))

    spa_cts <- spa_cts[,c('MarkerID','p.value', 'CHR')]
    spa_bin <- spa_bin[,c('MarkerID','p.value', 'CHR')]
    spa <- rbind(spa_cts, spa_bin)
    spa$p.value <- as.numeric(spa$p.value)
    NAs <- suppressWarnings(sum(is.na(spa$p.value)))
    if (NAs > 0) write(paste("[verbose]",NAs, "NAs in SPA P.values"), stderr())

    spa_chr <- spa[,c('MarkerID', 'CHR')]
    spa_chr <- spa_chr[!duplicated(spa_chr),]
    spa_aggr <- aggregate(p.value ~ MarkerID, data = spa, FUN = function(x) min(x, na.rm = TRUE))
    spa_aggr$permut <- ceiling(2/spa_aggr$p.value)

    tsv <- do.call(rbind, lapply(1:22, function(chr){
        path <- gsub('CHR',chr, args$tsv_path)
        fread(path)
    }))

    tsv <- tsv[,c('gene_id','varid')]
    tsv <- tsv[!duplicated(tsv),]
    tsv$varid <- clean_hail_list(tsv$varid)

    aggr <- aggregate(varid ~ gene_id, data = tsv, FUN = function(x) paste(x, collapse = ', '))
    aggr$varid <- unlist(lapply(strsplit(aggr$varid, split = ', '), function(x) paste0(unique(x), collapse = ', ')))
    colnames(aggr) <- c('MarkerID','varid')

    spa_aggr <- merge(spa_aggr, aggr, all.x = TRUE) # combine SPA / tsv
    spa_aggr <- merge(spa_aggr, spa_chr, all.x = TRUE) # add chromsome
    spa_aggr <- spa_aggr[rev(order(spa_aggr$permut)),]

    outfile <- paste0(args$out_prefix, ".tsv.gz")
    write(paste("[verbose] writing",outfile),stderr())
    fwrite(spa_aggr, outfile, sep = '\t')



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--spa_cts_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--spa_bin_dir", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--tsv_path", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

