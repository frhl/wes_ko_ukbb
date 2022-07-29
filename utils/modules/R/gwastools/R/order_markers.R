#' @title Order markers by position
#' @param x vector of markers to be ordered
#' @return a vector of oredered markers.
#' @export

order_markers <- function(x, rm.dup = TRUE) {
    
    # check integrity
    regex <- "chr[0-9]+\\:[0-9]+\\:[a-zA-Z]+\\:[a-zA-Z]"
    ok <- grepl(regex, x)
    if (any(!ok)){
      example <- x[!ok][1]
      n <- sum(!ok)
      stop(paste0(n, " marker(s) were not formatted correctly, e.g. '", example,"'."))
    }

    # format/order by position
    d <- data.frame(do.call(rbind, strsplit(x, split = ':')))
    colnames(d) <- c("chr","pos","a1", "a2")
    d$pos <- as.numeric(d$pos)
    d$chr_int <- as.numeric(gsub("chr","",d$chr))
    d <- d[order(d$chr_int, d$pos),]
    d$chr_int <- NULL
    
    # deal with duplicates
    if (rm.dup){
      dups <- duplicated(d)
      d <- d[!dups,]
      if (sum(dups) > 0) warning("Some markers are duplicated. These have been discarded. Change with 'rm.dup'")
    }
 
    # format out vector
    out <- apply( d , 1 , paste , collapse = ":")
    out <- as.character(gsub(" ", "", out))
    return(out)
}


