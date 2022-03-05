#' @title predict polygenic risk score
#' @param gwas a data.table/frame containing QCed GWAS summary stats
#' @param effects corresponding variant effects after using snp_ldpred2_inf
#' @param obj a bigsnp read in with \code{load_bigsnp_from_bed}.
#' @param ind_row indicies for samples that should be extracted. Default (NULL) means all samples.
#' @param ncores integer. cores to be used.
#' @export

predict_prs <- function(obj, gwas, effects, ind_row = NULL, ncores = 1){
    
    # assume that effect index correspond to gwas index
    stopifnot(length(effects) == nrow(gwas))
    
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
    
    # get subset of effects corresponding to bed_map
    ld_pred_effects <- effects[indicies]
    
    # subset bigsnp to get the right genotypes
    samples <- obj$bigsnp$fam$sample.ID
    if (is.null(ind_row)) ind_row <- 1:length(samples)
    subsetted_bigsnp <- snp_subset(obj$bigsnp, ind.col = bed_map$index, ind.row = ind_row)
    subsetted_bigsnp <- snp_attach(subsetted_bigsnp)
    genotypes <- subsetted_bigsnp$genotypes     
    stopifnot(dim(genotypes)[2] == length(ld_pred_effects))

    # Get PRS
    dot_product <- big_prodVec(genotypes, ld_pred_effects, ncores = ncores)
    return(list(prs = dot_product, sid =  samples[ind_row]))
}



