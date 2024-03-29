#' @title match genotype rows with gwas rows
#' @param obj a bigsnp obj (see library bigsnpr)
#' @param gwas a gwas summary statistics data.frame
#' @param ind_row sequence of numbers for samples to select.
#' @param bfile backing file (warning, this file will take up a lot of space! so delete)
#' @param ncores how many cores should be used.
#' @export

match_bigsnp_with_gwas <- function(obj, gwas, ind_row = NULL, bfile = NULL, ncores = 1){

    # add indexes used for later subsetting
    bed_map <- obj$map
    bed_map$marker <- get_ldpred_marker(bed_map)
    bed_map$index <- 1:nrow(bed_map)
    gwas$marker <- get_ldpred_marker(gwas)
    gwas$index <- 1:nrow(gwas)

    # check that all markers are contained in one another
    stopifnot(sum(gwas$marker %in% bed_map$marker) == sum(bed_map$marker %in% gwas$marker))

    # remove any bed markers not in gwas (assuming
    # that the GWAS has already been qced before)
    bed_map <- bed_map[bed_map$marker %in% gwas$marker,]

    # match gwas betas with bed map
    indicies <- match(bed_map$marker, gwas$marker)
    gwas_matched <- gwas[indicies,]
    gwas_matched <- gwas_matched[!(rowSums(is.na(gwas_matched)) == ncol(gwas_matched)),]

    # check that all positions match
    stopifnot(all(gwas_matched$pos == bed_map$pos))

    # subset bigsnp to get the right genotypes
    samples <- obj$bigsnp$fam$sample.ID
    if (is.null(ind_row)) ind_row <- 1:length(samples)
    subsetted_bigsnp <- snp_subset(obj$bigsnp, 
                                   ind.col = bed_map$index, 
                                   ind.row = ind_row,
                                   backingfile = bfile
                                   )
    subsetted_bigsnp <- snp_attach(subsetted_bigsnp)
    genotypes <- subsetted_bigsnp$genotypes
    stopifnot(dim(genotypes)[2] == nrow(gwas_matched))

    return(list(genotypes = genotypes, map = gwas_matched, sid =  samples[ind_row], gwas_indicies = indicies))
}


