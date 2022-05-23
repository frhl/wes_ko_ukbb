#' @title inverse rank normalization
#' @param x numerics, typically phenotype values.
#' @return a vector of inverse rank noramlized phenotype values.

inv_rank_norm <- function(x) {
    return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
}


