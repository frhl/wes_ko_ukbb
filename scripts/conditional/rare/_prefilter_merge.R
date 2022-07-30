#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    stopifnot(dir.exists(dirname(args$out_prefix)))

    # setup input
    files <- list.files(
        dirname(args$in_prefix), 
        pattern = basename(args$in_prefix), 
        full.names = TRUE
    )
    files <- files[grepl("chr", files)]
    files_hash <- files[grepl("hash.txt.gz",files)]
    files_ac <- files[grepl("AC.txt.gz",files)]

    # get columns with IDs
    col_id <- c("chr","pos","id","ref","alt")

    # get starting columns with markers
    markers <- fread(files[[1]])[,col_id,with=FALSE]
    marker_id <- markers$id

    # get hashes
    lst_hash <- lapply(files_hash, function(f){
        d <- fread(f)
        stopifnot(all(marker_id == d$id))
        d <- d[,-col_id,with=FALSE]
        unlink(f)
        return(d)
    })

    # get AC
    lst_ac <- lapply(files_ac, function(f){
        d <- fread(f)
        stopifnot(all(marker_id == d$id))
        d <- d[,-col_id,with=FALSE]
        unlink(f)
        return(d)
    })

    # combine them all
    d_hash <- cbind(markers, do.call(cbind, lst_hash))
    d_ac <- cbind(markers, do.call(cbind, lst_ac))

    # write to disk
    outfile_hash <- paste0(args$out_prefix, "_hash.txt.gz")
    outfile_ac <- paste0(args$out_prefix, "_AC.txt.gz")
    fwrite(d_hash, outfile_hash, sep = "\t")
    fwrite(d_ac, outfile_ac, sep = "\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_prefix", default=NULL, required = TRUE, help = "Prefix of path")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix")
args <- parser$parse_args()

main(args)









