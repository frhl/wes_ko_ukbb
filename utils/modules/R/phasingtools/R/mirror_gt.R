#' @title mirror genotypes by flipping alternate allele
#' @param x string or vector of genotype string, e.g. "1|0".
#' @export

# mirror genotype, e.g. mirror_gt(c("1|0","0|1"))
mirror_gt <- function(x){
    stopifnot(grepl("[0-1]\\|[0-1]",x))
    gt1 <- substring(x, 1,1)
    gt2 <- substring(x, 3,3)
    return(paste0(gt2,"|",gt1))
}




