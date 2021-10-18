# extract hail field
extract_hail_field <- function(data){
    
    n = length(data)
    # create a sub matrix for each row
    lst <- lapply(1:n, function(i){
        row <- strsplit(data[i], split = '\\,')[1]
        unlisted_row <- unlist(row)
        count_skip <- lengths(regmatches(unlisted_row, gregexpr(":", unlisted_row)))
        unlisted_row <- unlisted_row[count_skip == 1]
        matrix <- matrix(unlist(lapply(list(unlisted_row), function(x) strsplit(x, split = '\\:'))), nrow = 2)
        entries <- data.table(t(matrix[2,]))
        colnames(entries) <- matrix[1,]
        return(entries)
                                   
    })
                                
    # note that sub-matrices have different column names
    combined <- suppressMessages(dplyr::bind_rows(lst))
    combined <- combined[,!grepl('(name\\.\\.)|(db\\.\\.)',colnames(combined)), with = FALSE]
    return(combined)
} 



