#' @title get interval overlaps
#' @param x1 range 1 start
#' @param x2 range 1 end
#' @param y1 range 1 start
#' @param y2 range 2 end
#' @return a vector indexes for x in y

which_overlap <- function(x1, x2, y1, y2){
    # check overlap 
    is_overlapping <- function(x1,x2,y1,y2){
        return(max(x1,y1) <= min(x2,y2))
    }
    nx <- length(x1)
    ny <- length(y1)
    index <- lapply(1:nx, function(i){
        index <- which(unlist(lapply(1:ny, function(j){
            x1i <- x1[i]
            x2i <- x2[i]
            y1j <- y1[j]
            y2j <- y2[j]
            stopifnot(x1i <= x2i)
            stopifnot(y1j <= y2j)
            overlap <- is_overlapping(x1i, x2i, y1j, y2j)
        })))
        return(unlist(index))
    })
    # ensure that missing indexes are NAs
    index[!(index %in% 1:ny)] <- NA
    return(unlist(index))
    }




