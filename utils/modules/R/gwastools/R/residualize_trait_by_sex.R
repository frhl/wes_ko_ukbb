#' @title residualize trait seperately for sexes
#' @description fit fixed effects and residualize the trait 
#' @param trait string, trait of interest.
#' @param covars vector of fixed effects to be included
#' @param data phenotype data
#' @param sex string, column in data that contains sex information (expect 1 or 0)
#' @references From workflow in \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6488973/}
#' @return a vector of residuals

residualize_trait_by_sex <- function(trait, data, sex = "sex", covars = c("age", "age2", "sex", "sex_age")){
    stopifnot(sex %in% colnames(data))
    stopifnot(data[[sex]] %in% c(0,1))
    bool_sex0 <- data[[sex]] == 0
    bool_sex1 <- data[[sex]] == 1
    resids <- rep(NA, nrow(data))
    resids[bool_sex0] <- residualize_trait(trait, data[bool_sex0], covars)
    resids[bool_sex1] <- residualize_trait(trait, data[bool_sex1], covars)
    return(resids)
}


