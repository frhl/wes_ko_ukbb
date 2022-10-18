
library(data.table)

main <- function(args){

    interval_idx = args$interval_idx
    interval_path = args$interval_path
    read_placeholder = args$read_placeholder
    stopifnot(file.exists(interval_path))
    d <- readLines(interval_path)
    idx <- as.numeric(interval_idx) 
    samples <- unlist(strsplit(d[idx], split = ','))
    stopifnot(length(samples) > 1)
    # change all occurences of SAMPLE to actual sampleID
    sample_read_path <- unlist(lapply(samples, function(s) {
        readpath <- gsub("SAMPLE", s, read_placeholder)
        stopifnot(file.exists(readpath))
        return(readpath)
    }))
    # collapse and write out
    out <- paste(sample_read_path, collapse = " ")
    write(out, stdout())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--interval_idx", default=NULL, required = TRUE, help = "")
parser$add_argument("--interval_path", default=NULL, help = "")
parser$add_argument("--read_placeholder", default=NULL, help = "")
args <- parser$parse_args()

main(args)


