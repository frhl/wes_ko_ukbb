
#' @title get path to pheno header
#' @return path to file containing phenotypes tested
get_phenos_header_path <- function(){
    path <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/phenotypes/dec22_phenotypes_binary_200k_header.tsv"
    stopifnot(file.exists(path))
    return(path)
}

#' @title get path to phenos run with PRS
#' @param use_bonf_corrected boolean indicating whether only
#' phenotypes passing bonferroni correction should be returned
#' @return path to file containg PRS that have been run
get_phenos_prs_path <- function(use_bonf_corrected = TRUE){
    dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/validation/"
    if (use_bonf_corrected) path <- file.path(dir,"data/prs/validation/ldsc_summary_keep_phenos.txt")
    stopifnot(file.exists(path))
    return(path)
}


#' @title get phenotype from path
#' @param string or vector of paths
#' @return string or vector of phenotypes
gsub_phenotype_from_path <- function(paths){
    sub <- gsub("(ukb_eur_wes_200k_)|(_pLoF_dam.*)", "", basename(paths))
    stopifnot(length(sub) > 0)
    invalid_phenos <- sum(!sub %in% get_phenos_tested())
    if (invalid_phenos > 0) warning(paste(invalid_phenos, "invalid phenotypes."))
    return(sub)
}

#' @title get phenotypes
#' @param prs string. Either "with" or "without". Default (NULL) 
#' @param use_bonf_corrected bool. Only use bonferroni corrected phenotypes
#' will include both with and without PRS
#' @return vector of phenotypes tested
get_phenos_tested <- function(prs=NULL, use_bonf_corrected=TRUE){
    phenos <- readLines(get_phenos_header_path())
    if (!is.null(prs)){
        phenos_w_prs <- readLines(get_phenos_prs_path(use_bonf_corrected))[-1]
        if (prs == "with"){
            phenos <- phenos_w_prs
        } else if (prs == "without"){
            phenos <- phenos[!phenos %in% phenos_w_prs]
        } else { 
            stop("Argument 'prs' should be 'with' or 'without' or 'NULL'.")
        }
    }
    return(phenos)
}

#' @title get vector of booleans in which we include PRS
#' @param phenotypes vector of phenotypes 
#' @param use_bonf_corrected bool. Only use bonferroni corrected phenotypes
#' @return boolean indicating whether we condition on PRS
grepl_cond_prs <- function(phenotypes, use_bonf_corrected=TRUE){
    phenos_all <- get_phenos_tested(prs=NULL)
    phenos_prs <- get_phenos_tested(prs="with", use_bonf_corrected)
    n <- sum(!phenotypes %in% phenos_all)
    if (n>0) warning(paste(n, "invalid phenotype(s) included."))
    return(phenotypes %in% phenos_prs)    
}





