#' @title quality control of binary summary statistics
#' @param G genotypes from bigsnpr
#' @param info_snp info from bigsnpr
#' @export

qc_sumstat_binary <- function(G, info_snp, ncores = 1) {
    maf <- snp_MAF(G, ind.col = info_snp$`_NUM_ID_`, ncores = ncores)
    sd_val <- sqrt(2 * maf * (1 - maf))
    # see section 3.4 of https://academic.oup.com/bioinformatics/article/36/22-23/5424/6039173?login=true#233620521
    sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    is_bad_sd <-
      sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
    is_bad_maf <- is.na(maf)
    is_bad <- (is_bad_sd | is_bad_maf)
    sum_is_bad <- sum(is_bad, na.rm = TRUE)
    sum_n <- length(is_bad)
    write(paste(sum_is_bad, 'of', sum_n, 'SNPs fail QC'), stderr())
    return(list(is_bad = is_bad, sd_val = sd_val, sd_ss = sd_ss))
}

