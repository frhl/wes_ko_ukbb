
#' @title get path to ldsc validation file
#' @return path to file containing phenotypes teste
get_ldsc_path <- function(){
    path <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/validation/ldsc_summary.txt.gz"
    stopifnot(file.exists(path))
    return(path)
}

#' @title get ldsc heritability estimate
#' @param traits a vector of traits
#' @param get either "h2", "p", "std_error" or "n_eff"
#' @return estimate for each phenotype
get_ldsc_param <- function(traits, get = "h2"){
    d <- fread(get_ldsc_path())
    d <- d[d$coef == "h2",]
    if (get == "h2"){
        mapping <- d$estimate
        names(mapping) <- d$phenotype
    } else if (get == "p"){
        mapping <- d$pvalue
        names(mapping) <- d$phenotype
    } else if (get == "std_error"){
        mapping <- d$std_error
        names(mapping) <- d$phenotype
    } else if (get == "n_eff"){
        mapping <- d$n_eff
        names(mapping) <- d$phenotype
    }
    return(mapping[traits])
}




