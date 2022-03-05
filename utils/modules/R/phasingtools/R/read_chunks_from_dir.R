#' @title read and append a directory of chunks
#' @param directory string to directory
#' @export

read_chunks_from_dir <- function(directory, regex = NULL){
    files <- list.files(directory, full.names = TRUE, pattern = '.txt')
    if (!is.null(regex)) files <- files[grepl(regex, files)]
    if (length(files) == 0) stop(paste("no files in", directory))
    lst <- lapply(files, function(f) {
        d <- fread(f)
        chunk <- extract_chunk_id(f)
        d$chunk <- chunk
        d$cumsum <- cumsum(d$switches)
        d$dataset <- basename(f) 
        return(d)
    })
    df <- do.call(rbind, lst)
    return(df)
}


