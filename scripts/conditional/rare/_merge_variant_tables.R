#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    stopifnot(dir.exists(dirname(args$out_prefix)))
    expected_chunks <- as.numeric(args$chunks) 
    stopifnot(expected_chunks > 0)

    # setup input
    files <- list.files(
        dirname(args$in_prefix), 
        pattern = basename(args$in_prefix), 
        full.names = TRUE
    )

    # regex to correct files
    files <- files[grepl("chr", files)]
    files_hash <- files[grepl("hash.txt.gz",files)]
    files_ac <- files[grepl("AC.txt.gz",files)]

    lst_hash <- lapply(files_hash, fread)
    lst_ac <- lapply(files_ac, fread)
    
    # ensure all files are recovereed
    n_hash <- length(lst_hash)
    n_ac <- length(lst_ac)
    stopifnot(n_hash == n_ac)
    all_chunks_found <- n_hash >= expected_chunks
    if (!all_chunks_found) stop(paste("Found", n_hash, "chunks but expected", expected_chunks, "for", args$in_prefix))
    

    # combine them all
    d_hash <- rbindlist(lst_hash)
    d_ac <- rbindlist(lst_ac)

    # remove overlaps
    d_hash <- d_hash[!duplicated(d_hash),]
    d_ac <- d_ac[!duplicated(d_ac),]

    # write to disk
    outfile_hash <- paste0(args$out_prefix, "_hash.txt.gz")
    outfile_ac <- paste0(args$out_prefix, "_AC.txt.gz")
    write(paste("writing", outfile_hash, "and", outfile_ac), stdout())
    fwrite(d_hash, outfile_hash, sep = "\t")
    fwrite(d_ac, outfile_ac, sep = "\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_prefix", default=NULL, required = TRUE, help = "Prefix of path")
parser$add_argument("--chunks", default=NULL, required = TRUE, help = "Expected number of chunks")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix")
args <- parser$parse_args()

main(args)

