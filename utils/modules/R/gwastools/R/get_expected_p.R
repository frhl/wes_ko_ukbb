#' @title Get expected P-value based on uniform quantiles
#' @param observed_p float for current p-values
#' @param na.rm Should NAs be discarded?
#' @exporti

get_expected_p <- function(observed_p, na.rm = FALSE){
    stopifnot(is.numeric(observed_p))
    stopifnot(any(observed_p >= 0))
    sum_na <- sum(is.na(observed_p))
    if (sum_na > 0){
           if (na.rm){
                observed_p <- observed_p[!is.na(observed_p)]
           } else {
                warning("NAs detected in P-values. Set na.rm = TRUE to remove them.")
           }
    }
    n <- length(observed_p)
    observed_rank <- rank(observed_p)
    uniform <- (1:n)/(n+1)
    uniform <- uniform[observed_rank]
    return(uniform)
}


