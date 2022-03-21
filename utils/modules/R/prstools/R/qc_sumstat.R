#' @title quality control of binary summary statistics
#' @param G genotypes from bigsnpr
#' @param info_snp info from bigsnpr
#' @param trait either "binary" or"cts"
#' @param sd_y standard deviation of cts phenotype (assuming inverse rank normalisation)
#' @export

qc_sumstat <- function(G, info_snp, n_eff, trait, ncores = 1, sd_y = NULL)
    stopifnot(trait %in% c("binary", "cts"))
    maf <- snp_MAF(G, ind.col = info_snp$`_NUM_ID_`, ncores = ncores)
    sd_val <- sqrt(2 * maf * (1 - maf))
    # see section 3.4 of https://academic.oup.com/bioinformatics/article/36/22-23/5424/6039173?login=true#233620521
    if (trait == "binary") {
      sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2 + beta^2))
    } else {
      # estimate sd(y) from summary statistics
      #sd_y <- min(with(info_snp, sqrt(0.5) * beta_se^2 * n_eff)) 
      stopifnot(!is.null(sd_y))
      sd_ss <- with(info_snp, sd_y / sqrt(n_eff * beta_se^2 + beta^2))
    }
    is_bad_sd <-
      sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
    is_bad_maf <- is.na(maf)
    is_bad <- (is_bad_sd | is_bad_maf)
    sum_is_bad <- sum(is_bad, na.rm = TRUE)
    sum_n <- length(is_bad)
    write(paste(sum_is_bad, 'of', sum_n, 'SNPs fail QC'), stderr())
    return(list(is_bad = is_bad, sd_val = sd_val, sd_ss = sd_ss))
}

