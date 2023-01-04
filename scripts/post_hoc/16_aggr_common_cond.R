
library(data.table)
library(argparse)

main <- function(args){

    stopifnot(dir.exists(dirname(args$markers_dir)))
    files  <- list.files(args$markers_dir, full.names = TRUE, pattern = ".markers$")
    genes <- stringr::str_extract(files, "ENSG[0-9]+")
    traits <- stringr::str_extract(files, "200k_([0-9]|[A-Z]|[a-z]|\\_)+_pLoF")
    traits <- gsub("(200k_)|(_pLoF)", "", traits)
    
    # get a few stats
    stats <- data.table(
        genes = genes,
        traits = traits
    )

    # open files and extract conditional markers
    open_files <- do.call(rbind, lapply(files, function(f){
        d <- fread(f)
        d$converged <- (d$V1<30)
        d <- d[nrow(d),c(1,2,4,5,8,7)]
        colnames(d) <- c("iteration", "current_marker", "last_p","cutoff_p", "converged","conditioning_markers")
        return(d)
    }))

    # combine the files
    final <- cbind(stats, open_files)
    final$files <- files

    # out
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(final, outfile, sep = "\t", na="NA", quote=FALSE)
        

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--markers_dir", default=NULL, required = TRUE, help = "Directory for markers")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


