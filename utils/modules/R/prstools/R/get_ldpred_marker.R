#' @title get ldpred marker column from other columns
#' @param d data.frame/table with columns chr, pos, a0, a1
#' @export

get_ldpred_marker <- function(d){
    return(
        paste(
        d$chr, 
        d$pos, 
        d$a0, 
        d$a1, 
        sep = ':')
    )
}


