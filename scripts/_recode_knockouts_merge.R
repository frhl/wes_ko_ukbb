#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    stopifnot(dir.exists(args$in_dir))
    files <- list.files(args$in_dir, pattern = args$regex, full.names = TRUE)
    suffix <- args$suffix
    if (!is.null(suffix)) files <- files[grepl(suffix,files)]
    lst <- lapply(files, fread)
    stopifnot(length(lst) == as.numeric(args$n_expt))
    d <- rbindlist(lst, fill = TRUE)
    outfile <- args$out
    write(paste0("writing merge to ", outfile), stderr())
    fwrite(d, outfile, sep = "\t", na = NA, quote = FALSE)

    # cleaning up the chunks
    for (f in files){
        if (file.exists(f)){
            write(paste0("Removing chunk ", basename(f)), stderr())
            unlink(f)
        }
    }

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "directory of input file")
parser$add_argument("--regex", default=NULL, help = "basename of input file")
parser$add_argument("--suffix", default=NULL, help = "?")
parser$add_argument("--out", default=NULL, help = "?")
parser$add_argument("--n_expt", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

