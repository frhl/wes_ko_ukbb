#' @title unpack co-occuring variants
#' @param d a data.table
#' @return a list of list for co-occurecne
#' @export

cooccur_unpack_lib <- function(lib){
    genes <- names(lib)
    unpacked <- rbindlist(lapply(genes, function(g){
        variants <- names(lib[[g]])
        rbindlist(lapply(variants, function(v){
            # extract variants
            splitted_variants <- unlist(strsplit(v, split = "-"))
            if (length(splitted_variants) == 1){
                v1 <- splitted_variants[1]
                v2 <- splitted_variants[1]
            } else {
                v1 <- splitted_variants[1]
                v2 <- splitted_variants[2]
            }

            # combine
            config <- lib[[g]][[v]]
            occ <- data.table(do.call(cbind, config))
            occ <- cbind(g, v1, v2, occ)
            return(occ)
        }))
    }))
    return(unpacked)
}


