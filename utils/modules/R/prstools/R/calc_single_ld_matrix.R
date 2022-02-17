#' @title estimate LD matrix 
#' @param G genotypes from obj.bigSNP$genotypes
#' @param POS2 vector of cM distances
#' @param info_snp result from runningn bigsnpr::snp_match between
#' a reference panel and the summary statistics.
#' @param chrs what chromosomes should be evaluated? 
#' @note Chromosome 1 must always be included. 
#' @export


calc_single_ld_matrix <- function(G, POS2, info_snp, chr, ncores = 1, tmp){
    
    write(paste0("Fitting LD for", chr, "using", ncores,"core(s).."), stderr()) 
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    stopifnot(length(ind.chr) > 0)
    corr0 <- snp_cor(
            G,
            ind.col = ind.chr2,
            ncores = ncores,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
    return(list(ld = ld, corr = corr))
}



