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
    path <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/prs/validation/ldsc_summary_nom_sig_phenos.txt"
    if (use_bonf_corrected) path <- "data/prs/validation/ldsc_summary_bonf_sig_phenos.txt"
    stopifnot(file.exists(path))
    return(path)
}

#' @title get phenotypes
#' @param prs string. Either "with" or "without". Default (NULL) 
#' will include both with and without PRS
#' @return vector of phenotypes tested
get_phenos_tested <- function(prs=NULL){
    phenos <- readLines(get_phenos_header_path())
    if (!is.null(prs)){
        phenos_w_prs <- readLines(get_phenos_prs_path(use_bonf_corrected = TRUE))[-1]
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





