#' @title list files/results after saige analysis
#' @param cond either "none", "common", "rare", "combined"
#' @param prs should PRS be included, either "all", "exclude", "only" or "prefer" which resturns PRS when available.
#' @param regex regex for saige files to grab
#' @return a vector of file paths to results

list_files_saige <- function(cond = NULL, prs = "include", regex = "\\.txt\\.gz"){

    # deal with old version
    if (cond == "none") cond <- ""
    if (is.null(cond)) cond <- ""

    # get phenotypes we are running
    trait_allow_prs <- get_phenos_tested(prs="with", use_bonf_corrected = TRUE)
    trait_disallow_prs <- get_phenos_tested(prs="without")
    all_traits <- c(trait_allow_prs, trait_disallow_prs)

    # get files that we have created
    files <- sort(list.files(get_saige_dir(cond), full.names = TRUE, pattern = ".txt.gz"))
    files_pheno <- gsub_phenotype_from_path(files)
    files_prs <- grepl("locoprs.txt.gz", files)
    allow_prs <- files_pheno %in% trait_allow_prs

    # combine into data.frame
    df <- data.frame(
      phenotype = files_pheno,
      prs_available = files_prs,
      prs_allowed = allow_prs,
      path = files
    )

    # subset phenotype
    df <- df[df$phenotype %in% all_traits,]

    # list all files with PRS
    if (prs %in% "include") {
        df <- df[(df$prs_allowed) | (!df$prs_available),]
    } else if (prs %in% "exclude") {
        df <- df[!df$prs_available,]
    } else if (prs %in% "only") {
        df <- df[(df$prs_available) & (df$prs_allowed),]
    } else if (prs %in% "prefer") {
        df1 <- df[(df$prs_available) & (df$prs_allowed),]
        df2 <- df[(!df$prs_available) & (!df$prs_allowed),]
        df <- rbind(df1, df2)
    } else {
        stop(paste(prs, "is not valid. Must be either 'include','exclude','prefer' or 'only'."))
    }
    return(df$path)

}


