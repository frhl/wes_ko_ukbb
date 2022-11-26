#' @title id single switch error
#' @param switches vector with switches encoded as 0 or 1
#' @export

id_single_switch_errors <- function(switches){
    stopifnot(all(switches) %in% c(1,0))
    n <- length(switches)
    stopifnot(n > 1)
    switches <- unlist(
        lapply(1:n, function(idx){
            idx_lower <- max(1, idx - 1)
            idx_upper <- min(n, idx + 1)
            switch_lower <- as.logical(switches[idx_lower])
            switch <- as.logical(switches[idx])
            switch_upper <- as.logical(switches[idx_upper])
            # first entry of switches
            if (idx_lower == idx){
                if (switch & !switch_upper){
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            # last entry of switches
            } else if (idx_upper == idx){
                if (switch & !switch_lower){
                    return(TRUE)
                } else {
                    return(FALSE)
                }
            # not last or first entry of switches
            } else {
                if (switch & !switch_lower & !switch_upper){
                    return(TRUE)
                } else {
                    return(FALSE)
                }

            }
        })
    )
    return(switches)
}



