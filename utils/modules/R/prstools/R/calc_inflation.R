#' @title calculate lambda squared inflation statistic
#' @param P P-values from GWAS
#' @return Inflation statistics

calc_inflation <- function(P){
    return(median(qchisq(1-P,1)) / qchisq(0.5, 1))
}




