library(argparse)
library(stringr)

main <- function(args){
    
    files <- list.files(dirname(args$in_prefix), pattern = basename(args$in_prefix), full.names = TRUE)
    files <- files[grepl(".vcf.bgz$", files)]
    stopifnot(length(files) > 0)

    chunks <- stringr::str_extract(files,"[0-9]+of[0-9]+")
    chunks <- data.frame(do.call(rbind, strsplit(chunks, "of")))
    chunks$X1 <- as.numeric(chunks$X1)
    chunks$X2 <- as.numeric(chunks$X2)

    stopifnot(length(unique(chunks$X2)) == 1)
    stopifnot(length(unique(chunks$X1)) == nrow(chunks))
        
    out <- files[order(chunks$X1)]
    write(paste(out, collapse = " "), stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_prefix", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
args <- parser$parse_args()

main(args)


