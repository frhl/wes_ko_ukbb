#' @title estimate LD matrix 
#' @param G genotypes from obj.bigSNP$genotypes
#' @param POS2 vector of cM distances
#' @param info_snp result from runningn bigsnpr::snp_match between
#' a reference panel and the summary statistics.
#' @param chrs what chromosomes should be evaluated? 
#' @note Chromosome 1 must always be included. 
#' @export


calc_ld_matrix <- function(G, POS2, info_snp, chrs = 1:22, ncores = 1){
    
    chrs <- paste0("chr",chrs)
    for (chr in chrs) {
        write(paste('Evaluating',chr,'..'),stdout())
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
        if (chr == "chr1") {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp)
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
        }
    }
    return(list(ld = ld, corr = corr))
}



