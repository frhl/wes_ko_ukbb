#' @title extract the current chunk number based on filename
#' @param fname string corresponding to filename)
#' @export


extract_chunk_id <- function(fname){
      return(as.numeric(unlist(strsplit(str_extract(fname,'[0-9]+of[0-9]+'), split = 'of'))[1]))
}
