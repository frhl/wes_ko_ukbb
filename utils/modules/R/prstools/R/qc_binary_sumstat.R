#' @title quality control of binary summary statistics
#' @param G genotypes from bigsnpr
#' @param info_snp info from bigsnpr
#' @export

qc_binary_sumstat <- function(G, info_snp){
    ss_trait <- ifelse(trait %in% binary, 0.5, 2)
    maf <- snp_MAF(G, ind.col = info_snp$`_NUM_ID_`, ncores = NCORES)
    sd_val <- sqrt(2 * maf * (1 - maf))
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    is_bad <-
      sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
    return(list(is_bad = is_bad, sd_val = sd_val, sd_ss = sd_ss))
}

