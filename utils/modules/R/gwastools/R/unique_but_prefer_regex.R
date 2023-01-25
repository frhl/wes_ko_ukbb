#' @title Get unique files but prefer files with regex
#' @param files list of files
#' @param regex if there are two files differing only by the string
#' in the regex, then the files containing the regex will be selected
#' and the other discarded.
#' @param ext extension
#' @export

unique_but_prefer_regex <- function(files, regex = "_locoprs", ext = ".txt.gz"){
    
    # do some house keeping on filenames
    directory <- unique(dirname(files))
    stopifnot(length(directory) == 1)
    file_df <- data.frame(bname=basename(files))
    file_df$sans_ext <- gsub(ext, "", file_df$bname)
    file_df$sans_str <- gsub(regex, "", file_df$sans_ext)
    
    # count how many occours twice
    counts <- data.frame(table(file_df$sans_str))
    colnames(counts) <- c("sans_str", "n")
    
    # if there is both a prs and non-prs file, take the PRS file
    file_df <- merge(file_df, counts)
    if (any(file_df$n > 2)) stop("Some files with regex appear more than twice!")
    file_df_n1 <- file_df[file_df$n == 1,]
    file_df_n2 <- file_df[file_df$n >= 2,]
    file_df_n2 <- file_df[grepl(regex, file_df$bname), ]
    file_df <- rbind(file_df_n1, file_df_n2)
    files <- file.path(directory, file_df$bname)  
    
}






