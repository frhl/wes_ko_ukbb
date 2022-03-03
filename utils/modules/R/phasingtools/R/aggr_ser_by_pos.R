#' @title aggregate switch errors by site
#' @param data.frame (usually created with hacked trio-stat plugin)
#' @export

aggr_ser_by_site <- function(d){
    # Aggregate switch errors by position 
    
    stopifnot('switches' %in% colnames(d))    
    stopifnot('POS' %in% colnames(d))    
    stopifnot('MAF' %in% colnames(d))    
    stopifnot('HWE' %in% colnames(d))    

    d1 <- aggregate(switches ~ POS, data = d, FUN = sum)
    d2 <- aggregate(switch ~ POS, data = d, FUN = length)
    d3 <- d[,c('POS','MAF','HWE')]
    d3 <- d3[!duplicated(d3),]
    
    colnames(d1) <- c('POS','n_switch')
    colnames(d2) <- c('POS','n_tested')
    mrg <- merge(merge(d1, d2), d3)
    mrg$n_switch_div_n_tested <- mrg$n_switch/mrg$n_tested
    return(mrg)
    
}


