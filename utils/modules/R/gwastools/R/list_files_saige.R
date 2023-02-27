#' @title list files/results after saige analysis
#' @param cond either "none", "common", "rare", "combined"
#' @param prs should PRS be included, either "all", "exclude", "only" or "prefer" which resturns PRS when available.
#' @param regex regex for saige files to grab
#' @return a vector of file paths to results

list_files_saige <- function(dir = "", prs = "all", regex = "\\.txt.\\.gz"){

    # get phenotypes we are running
    pheno_w_prs <- get_phenos_tested(prs="with")
    pheno_wo_prs <- get_phenos_tested(prs="without")
    phenos <- c(pheno_w_prs, pheno_wo_prs)
    
    # get files that we have created
    files <- sort(list.files(get_saige_dir(dir), full.names = TRUE, pattern = regex))
    files_w_prs <- files[grepl("locoprs.txt.gz", files)]
    files_wo_prs <- files[!grepl("locoprs.txt.gz", files)]
    
    # find name interesection between tested phenos and created files
    files_tested <- unique(unlist(lapply(phenos, function(p) extract_path_by_phenotype(p, files))))
    files_tested_only_prs <- unique(unlist(lapply(pheno_w_prs, function(p) extract_path_by_phenotype(p, files_w_prs))))
    files_tested_only_wo_prs <- unique(unlist(lapply(pheno_wo_prs, function(p) extract_path_by_phenotype(p, files_wo_prs))))
                         
    # list all files with PRS
    if (prs %in% "all") {
        warning("Listing all files including files with PRS that does not pass cutoffs!")
        return(files_tested)
    }
    else if (prs %in% "exclude") {
        return(files_tested_only_wo_prs)
    }
    else if (prs %in% "only") {
        return(files_tested_only_prs)
    }
    else if (prs %in% "prefer") {
        return(c(files_tested_only_prs, files_tested_only_wo_prs))
    }
    else {
        stop(paste(prs, "is not valid. Must be either 'include','exclude','prefer' or 'only'."))
    }
}

# helper to find subsets of file paths and phenotypes that have
# been included in the downstream analysis
extract_path_by_phenotype <- function(phenotype, paths) {
    stopifnot(length(paths) > 1)
    regex <- paste0("ukb_eur_wes_200k_", phenotype, "_pLoF_damaging_missense")
    paths <- paths[grepl(pattern = regex, paths)]
    return(paths)
}
                                      


