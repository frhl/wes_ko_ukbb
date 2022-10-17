#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(stringr)

main <- function(args){

    # parser
    print(args)
    stopifnot(!file.exists(args$out_path))
    AUTOSOMES <- paste0("chr",1:22)

    # read all autosomes
    d <- do.call(rbind, lapply(AUTOSOMES, function(chrom){
      fpath <- gsub("chrCHR",chrom, args$in_prefix)
      if (!file.exists(fpath)) stop(paste(fpath, "does not exist!"))
      d <- fread(fpath) 
      d <- d[d$pKO > 0 & d$pKO < 1]
      return(d)
    }))

    # clean up 
    d$varid <- gsub('(\\[)|(\\])|(\\")', '', d$varid)
    d$gts <- gsub('(\\[)|(\\])|(\\")', '', d$gts)
    d$chromosome <- str_extract(d$varid,"chr[0-9]+")
    d$knockout <- NULL

    # write combined SNPs 
    outfile = args$out_path
    write(paste("writing to", outfile), stdout())
    fwrite(d, outfile, sep = "\t") 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_prefix", default=NULL, required = TRUE, help = "Input with chrCHR format in file name.")
parser$add_argument("--out_path", default=NULL, required = TRUE, help = "Outout file destination?")
args <- parser$parse_args()

main(args)

