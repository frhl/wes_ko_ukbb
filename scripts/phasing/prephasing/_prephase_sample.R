# extract sample ID from chunk list

library(argparse)
library(data.table)

main <- function(args){

    chunk_idx = args$chunk_idx
    sample_idx = args$sample_idx
    interval_path = args$interval_path
    stopifnot(file.exists(interval_path))
    d <- readLines(interval_path)
    chunk_idx <- as.numeric(chunk_idx) 
    sample_idx <- as.numeric(sample_idx)
    samples <- unlist(strsplit(d[chunk_idx], split = ','))
    stopifnot(length(samples) >= 1)
    s <- samples[sample_idx] 
    out <- paste(s, collapse = " ")
    write(out, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chunk_idx", default=NULL, required = TRUE, help = "")
parser$add_argument("--sample_idx", default=NULL, required = TRUE, help = "")
parser$add_argument("--interval_path", default=NULL, help = "")
args <- parser$parse_args()

main(args)


