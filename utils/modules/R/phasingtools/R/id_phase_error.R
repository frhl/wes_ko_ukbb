#' @title id phase error
#' @param switches vector with switches encoded as 0 or 1
#' @note ID any single or double switch error
#' @export

id_phase_errors <- function(switches){
    n <- length(switches)
    stopifnot(switches %in% c(1,0))
    stopifnot(n > 1)
    switches <- unlist(
        lapply(1:n, function(idx){
            switch <- as.logical(switches[idx])
            if (switch) {
                return(TRUE)
            } else {
                return(FALSE)
            }
            })
        )
    return(switches)
    
}


