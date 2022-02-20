#' @title get LD matrix 
#' @param info_snp result from runningn bigsnpr::snp_match between
#' a reference panel and the summary statistics.
#' @param chrs what chromosomes should be evaluated? 
#' @param sfbm_file where should the LD-matrix backing file be stored?
#' @export

calc_ld_matrix2 <- function(info_snp, chrs = 1:22, sfbm_file = tempfile(), ld_dir = "data/prs/hapmap/ld/matrix"){

    chrs <- paste0("chr",chrs)
    first_chr <- chrs[1]
    for (chr in chrs) {
        write(paste0("Retriving LD for chrom ",chr),stderr())
        
        # paths to ld-matrix
        path_rds <- file.path(ld_dir, paste0("ld_matrix_",chr,".rda"))
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
        
        if (chr == first_chr) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, sfbm_file)
            outmap <- matched_map
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
            outmap <- rbind(outmap, matched_map)
        }
    }
    return(list(ld = ld, corr = corr, map = outmap))
}


