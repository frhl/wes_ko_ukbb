#' @title list files/results after saige analysis
#' @param cond either "none", "common", "rare", "combined"
#' @param prs should PRS be included, either "include", "exclude", "only" or "prefer" which resturns PRS when available.
#' @param regex regex for saige files to grab
#' @return a vector of file paths to results

list_files_saige <- function(
        cond="none", prs="include", regex = "\\.txt\\.gz",
        step2_dir = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2/min_mac4",
        step2_common_dir = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_common/min_mac4",
        step2_rare_dir = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_rare_cond/min_mac4",
        step2_combined_dir = "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_rare_cond/min_mac4"
    ){

    # subset paths based on condions
    if (cond %in% "none"){
        files <- list.files(step2_dir, full.names = TRUE, pattern = regex)
        if (length(files) == 0) stop(paste("No files found in", step2_dir))
    } else if (cond %in% "common"){
        files <- list.files(step2_common_dir, full.names = TRUE, pattern = regex)
        if (length(files) == 0) stop(paste("No files found in", step2_common_dir))
    } else if (cond %in% "rare"){
        files <- list.files(step2_rare_dir, full.names = TRUE, pattern = regex)
        if (length(files) == 0) stop(paste("No files found in", step2_rare_dir))
    } else if (cond %in% "combined"){
        files <- list.files(step2_combined_dir, full.names = TRUE, pattern = regex)
        if (length(files) == 0) stop(paste("No files found in", step2_combined_dir))
    } else {
        stop(paste(cond, "is not a valid. Try 'none','common','rare' or 'combined'"))
    }
    
    # perform final subset by PRS
    files <- sort(files)
    is_prs <- grepl("locoprs.txt.gz", files)
    if (prs %in% "include"){
        return(files)
    } else if (prs %in% "exclude"){
        # exclude any PRS files
        return(files[!is_prs])
    } else if (prs %in% "only") {
        # only include PRS
        return(files[is_prs])
    } else if (prs %in% "prefer") {
        # if two files exists, one with PRS and without
        # this option will preferentially select PRS
        directory <- unique(dirname(files))
        stopifnot(length(directory) == 1)
        file_df <- data.frame(bname=basename(files))
        file_df$sans_ext <- gsub(".txt.gz", "", file_df$bname)
        file_df$sans_locoprs <- gsub("_locoprs", "", file_df$sans_ext)
        # count how many occours twice
        counts <- data.frame(table(file_df$sans_locoprs))
        colnames(counts) <- c("sans_locoprs", "n")
        file_df <- merge(file_df, counts)
        # subset to files with locoprs if there are matches
        file_df_n1 <- file_df[file_df$n == 1,]
        file_df_n2 <- file_df[file_df$n >= 2,]
        file_df_n2 <- file_df[grepl("_locoprs", file_df$bname), ]
        file_df <- rbind(file_df_n1, file_df_n2)
        files <- file.path(directory, file_df$bname)
        return(files)
    } else {
        stop(paste(prs, "is not valid. Must be either 'include','exclude','prefer' or 'only'."))
    }

}


