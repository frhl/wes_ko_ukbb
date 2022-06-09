#' @title inverse normalization transformation (RINT)
#' @param x numerics, typically phenotype values.
#' @param k adjustable offset" (default: Blom offset = 3/8). 
#' @references \url{https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html#inverse-normal-transformation}
#' @return a vector of inverse rank noramlized phenotype values.


get_rint <- function(x, k = 3/8){

    stopifnot((k > 0) & (k<1/2))

    # keep original vector
    x_orig <- df$x
    defined <- !is.na(x_orig)
    x <- x_orig[defined]

    # samples and ties
    n <- length(x)
    xrank <- rank(x, ties.method = "min")

    # Transformation
    # Based on equation in this R function documentation:
    # https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html#inverse-normal-transformation
    
    x_rint <- qnorm( (xrank-k) / (n-2*k + 1))
    x_orig[defined] <- x_rint
    return(x_orig)
    
}
