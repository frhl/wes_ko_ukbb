#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    stopifnot(dir.exists(args$in_dir))
    files <- list.files(args$in_dir, pattern = args$regex, full.names = TRUE)
    files_index <- paste0(files, ".index")
    n <- length(files)
    if (n < 22) stop(paste0("Expected more than 21 files but found ", n, ". Merge regex: ", args$in_dir))    
    if (n > 23) stop(paste0("Expected less than24 files but found ", n, ". Merge regex: ", args$in_dir))
    # read files
    lst <- lapply(files, fread)
    d <- rbindlist(lst, fill = TRUE)
    # clean up files
    lapply(files, unlink)
    lapply(files_index, unlink)
    # order by chr and position
    d$chr_int <- as.numeric(gsub("chr","",d$CHR))
    d <- d[order(d$chr_int, d$POS),]
    d$chr_int <- NULL
    if (args$cond){
        


    }
    # setup outfile
    outfile <- args$out
    write(paste0("writing merge to ", outfile), stderr())
    fwrite(d, outfile, sep = "\t", na = NA, quote = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "directory of input file")
parser$add_argument("--cond", default=0, help = "Add columns to the data, so that cond P-values are always taken when available.")
parser$add_argument("--regex", default=NULL, help = "basename of input file")
parser$add_argument("--out", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

