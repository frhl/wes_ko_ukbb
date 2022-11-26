#' @title check if inptu is SNV based on refernce and alternate
#' @param ref vector of refence alleles
#' @param alt vector of alternate alleles
#' @export


is_snv <- function(ref, alt){
    stopifnot(!is.null(alt))
    stopifnot(!is.null(ref))
    stopifnot(length(ref) == length(alt))
    n <- length(ref)
    out <- unlist(lapply(1:n, function(idx){
        a1 <- ref[idx]
        a2 <- alt[idx]
        a1_nchar <- nchar(a1)
        a2_nchar <- nchar(a2)
        if ((a1_nchar == 1) & (a2_nchar == 1)){
            return(TRUE)
        } else {
            return(FALSE)
        }
    }))
    return(out)
}


