#' @title get LD matrix for single chromosome 
#' @param info_snp result from runningn bigsnpr::snp_match between
#' a reference panel and the summary statistics.
#' @param chr what chromosomes should be evaluated? 
#' @param sfbm_file where should the LD-matrix backing file be stored?
#' @export

get_single_ld_matrix <- function(info_snp, chr, sfbm_file = tempfile(), ld_dir = "data/prs/hapmap/ld/matrix"){

    # paths to ld-matrix
    path_rds <- file.path(ld_dir, paste0("ld_matrix_",chr,".rds"))
    path_map <- file.path(ld_dir, paste0("ld_matrix_",chr,".txt.gz"))
    
    # retrieve correlation matrix
    corr0 <- readRDS(path_rds)
    map <- fread(path_map)
    map$index <- 1:nrow(map)
    map$marker <- get_ldpred_marker(map)
    
    # subset ld-matrix
    i_index <- na.omit(match(info_snp$marker, map$marker)) 
    matched_map <- map[i_index,]
    j_index <- matched_map$index
    corr0 <- corr0[j_index, j_index]
    
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, sfbm_file)
    outmap <- matched_map

    return(list(ld = ld, corr = corr, map = outmap))
        
}


