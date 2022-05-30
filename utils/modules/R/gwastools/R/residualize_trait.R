#' @title residualize trait
#' @description fit fixed effects and residualize the trait 
#' @param trait string, trait of interest.
#' @param covars vector of fixed effects to be included
#' @param data phenotype data
#' @return a vector of residuals

residualize_trait <- function(trait, data, covars = c("age", "age2", "ukbb.centre")){
    stopifnot(!is.null(data))
    stopifnot(trait %in% colnames(data))
    stopifnot(all(covars %in% colnames(data)))
    covars <- paste0(covars, collapse = "+")
    #na_covar_rows <- rowSums(is.na(data[,covars, with = FALSE])) > 0
    #na_pheno_rows <- 
    f <- as.formula(paste0(trait, "~", covars))
    fit <- lm(f, data = data, na.action=na.exclude)
    return(resid(fit))
}

