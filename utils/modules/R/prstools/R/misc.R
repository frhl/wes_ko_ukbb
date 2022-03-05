#' @title file exists with extension
#' @param fname string. file name
#' @param ext string. extension
#' @return boolean 
#' @export

file.exists.ext <- function(fname, ext){
    fname <- tools::file_path_sans_ext(fname)
    fname <- paste0(fname, ext)
    return(file.exists(fname))
}


