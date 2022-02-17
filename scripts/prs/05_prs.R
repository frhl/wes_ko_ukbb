
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')


main <- function(args){

  stopifnot(file.exists(args$path_bed_pred))
  stopifnot(file.exists(args$path_bed_ld))
  stopifnot(file.exists(args$path_sumstat))

  # setup parallel environment
  NCORES <- nb_cores()
  tmp <- tempfile(tmpdir = "data/tmp/tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  # load data required for setting up LD-matrix 
  ld_data <- load_bigsnp_from_bed(args$path_ld_bed)

  # load summary statistics
  sumstats <- read_hail_sumstat(args$path_sumstat)
  
  # match summary stats and LD data
  info_snp <- snp_match(sumstats, ld_data$map, join_by_pos = TRUE, strand_flip = FALSE)
 
  # QC summary statistics based on LD reference
  qc <- qc_binary_sumstat(ld_data$G, info_snp)
  beta_cols <- c("beta", "beta_se", "n_eff", "_NUM_ID_")
  well_behaved_snps <- (!qc$is_bad)
  df_beta <- info_snp[well_behaved_snps, beta_cols]

  # get ld matrix. Note, that we need 60gb of memory to keep 
  # all the hapmap variants in memory
  snp <- calc_ld_matrix(ld_data$G, ld_data$POS2, info_snp, chrs = 1:22, 
                        ncores = NCORES)

  # perform ld regression
  ldsc <- with(df_beta, snp_ldsc(snp$ld, length(snp$ld), chi2 = (beta / beta_se)^2,
                                 sample_size = df_beta$n_eff, blocks = NULL))  
  ldsc_h2_est <- ldsc[["h2"]]
  
  # Load new data for prediction
  pred <- load_bigsnp_from_bed(args$path_bed_pred)

  # ldpred2 inf model
  beta_inf <- snp_ldpred2_inf(snp$corr, df_beta, ldsc_h2_est)
  final_pred_inf <- big_prodVec(pred$G, beta_inf, ind.col = df_beta$`_NUM_ID_`)  
  
  # ldpred2 auto model
  multi_auto <- snp_ldpred2_auto(
    corr, df_beta, h2_init = ldsc_h2_est,
    vec_p_init = seq_log(1e-4, 0.5, 30),
    ncores = NCORES)  

  # quality controls on chains
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  pred_auto <- big_prodMat(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]],
                           ncores = NCORES)
  sc <- apply(pred_auto, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto <- rowMeans(beta_auto[, keep])

  final_pred_auto <- big_prodVec(G, final_beta_auto,
                               ind.col = df_beta[["_NUM_ID_"]],
                               ncores = NCORES)
  # save results
  results <- data.table(
    sid = NA, 
    fid = NA,
    pred_auto = final_pred_auto,
    pred_inf = final_red_inf
    )

  #out_prefix <- tools::fil
  fwrite(resuts, file = paste0(out_prefix,".txt.gz"), sep = '\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bed_pred", default=NULL, required = TRUE, help = "Path for plink file (bed)")
parser$add_argument("--path_sumstat", default=NULL, required = TRUE, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--path_bed_ld", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









