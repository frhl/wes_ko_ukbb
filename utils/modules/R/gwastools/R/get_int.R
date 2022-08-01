#' @title inverse normalization transformation (INT)
#' @param x numerics, typically phenotype values.
#' @return a vector of inverse rank noramlized phenotype values.

get_int <- function(x) {
    return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}


