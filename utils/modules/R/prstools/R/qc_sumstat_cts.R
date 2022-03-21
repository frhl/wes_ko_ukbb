#' @title quality control of binary summary statistics
#' @param G genotypes from bigsnpr
#' @param info_snp info from bigsnpr
#' @param trait either "binary" or"cts"
#' @param sd_y standard deviation of cts phenotype (assuming inverse rank normalisation)
#' @export

qc_sumstat_cts <- function(G, info_snp, sd_y, ncores = 1)
    maf <- snp_MAF(G, ind.col = info_snp$`_NUM_ID_`, ncores = ncores)
    sd_val <- sqrt(2 * maf * (1 - maf))
    sd_ss <- with(info_snp, sd_y / sqrt(n_eff * beta_se^2 + beta^2))
    is_bad_sd <-
      sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
    is_bad_maf <- is.na(maf)
    is_bad <- (is_bad_sd | is_bad_maf)
    sum_is_bad <- sum(is_bad, na.rm = TRUE)
    sum_n <- length(is_bad)
    write(paste(sum_is_bad, 'of', sum_n, 'SNPs fail QC'), stderr())
    return(list(is_bad = is_bad, sd_val = sd_val, sd_ss = sd_ss))
}

