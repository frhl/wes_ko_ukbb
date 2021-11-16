#' @title fread with zcat
#' @param path path to read in
#' @param ... commands passend to fread
#' @export

zcat <- function(path,...){
    cmd = paste("zcat", path)
    return(fread(cmd = cmd,...))
}





