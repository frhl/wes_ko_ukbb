#' @title pack co-occuring variants
#' @description a function to zoom over a data.table and collect co-occuring variants.
#' note, currently a special case when >2 damaging variants occour in the same gene
#' then if on opposite haplotypes, those two are sampled, if on the same haplotype, 
#' two are selected at ranom
#' @param d a data.table
#' @return a list of list for co-occurecne
#' @export

cooccur_pack_lib <- function(d){
    
    stopifnot("gene_id" %in% colnames(d))
    stopifnot("varid" %in% colnames(d))
    stopifnot("gts" %in% colnames(d))
    stopifnot("is_cis" %in% colnames(d))
    stopifnot("is_chet" %in% colnames(d))
    stopifnot("is_hom" %in% colnames(d))
    
    lib <- list()
    for (idx in 1:nrow(d)){

        row <- d[idx,]
        gene_key <- row$gene_id
        row_variants <- row$varid
        haplotypes <- row$gts

        # add gene to library
        if (!gene_key %in% names(lib)){
            lib[[gene_key]] <- list(
            )
        }

        # combine splitted variants
        splitted_variants <- unlist(strsplit(row_variants, split = ";"))

        # deal with special case in which >2 variants per gene
        if (length(splitted_variants) > 2){

            # which variants are on opposite haplotypes?
            if (row$is_chet){
                splitted_haplotypes <- unlist(strsplit(haplotypes, split = ";"))
                idx_h1 <- which(splitted_haplotypes == "1|0")[1]
                idx_h2 <- which(splitted_haplotypes == "0|1")[1]
                var_h1 <- splitted_variants[idx_h1]
                var_h2 <- splitted_variants[idx_h2]
                splitted_variants <- c(var_h1, var_h2)
            # same haplotype
            } else {
                splitted_variants <- splitted_variants[1:2]
            }
        }

        # sort variants to make unique key regardless of haplotype
        sorted_variants <- sort(splitted_variants)
        variant_key <- paste0(sorted_variants, collapse = "-")

        # add variant key to library
        if (! variant_key %in% names(lib[[gene_key]])){
            lib[[gene_key]][[variant_key]] <- list(
                cis = 0, chet = 0, hom = 0
            )
        }

        # append to current. library
        if (row$is_cis) {
            lib[[gene_key]][[variant_key]]$cis <- lib[[gene_key]][[variant_key]]$cis + 1
        } else if (row$is_chet) {
            lib[[gene_key]][[variant_key]]$chet <- lib[[gene_key]][[variant_key]]$chet + 1
        } else if (row$is_hom) {
            lib[[gene_key]][[variant_key]]$hom <- lib[[gene_key]][[variant_key]]$hom + 1
        }
    }  
    return(lib)
}



