#' @title Get distance between single switches
#' @param switches vector of switches
#' @param pos vector of positions
#' @export

bp_between_single_switches <- function(switches, pos){
    stopifnot(!is.null(pos))
    stopifnot(length(switches) == length(pos))
    stopifnot(switches %in% c(1,0))
    # get position for single switches
    single_switch_idx <- which(id_single_switch_errors(switches))
    single_switch_pos <- pos[single_switch_idx]
    # get distance between each point
    n <- length(single_switch_pos)
    ranges <- unlist(lapply(1:(n-1), function(idx){
        x1 <- single_switch_pos[idx]
        x2 <- single_switch_pos[idx+1]
        distance <- x2-x1
        return(distance)
    }))
    # get position for all switch errors
    #metric <- ifelse(unit == "Mb", 1e+6, unit)
    return(ranges)
    
}


