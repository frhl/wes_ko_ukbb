#' @title read and append a directory of chunks
#' @param fname filename
#' @exporti

read_chunks_combined <- function(fname){
      d <- fread(fname)
      d$cumsum <- cumsum(d$switches)
      d$dataset <- basename(fname)
      return(d)
}

