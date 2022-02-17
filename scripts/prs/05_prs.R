
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  print(args)
  stopifnot(file.exists(args$path_bed_pred))
  stopifnot(file.exists(args$path_ld_matrix))
  stopifnot(file.exists(args$path_sumstat))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)

  # load data required for setting up LD-matrix 
  pred_data <- load_bigsnp_from_bed(args$path_bed_bed)

  # load summary statistics
  sumstats <- read_hail_sumstat(args$path_sumstat)
  
  # match summary stats with input data
  info_snp <- snp_match(sumstats, pred_data$map, join_by_pos = TRUE, strand_flip = FALSE)
 
  # QC summary statistics based on input reference
  qc <- qc_binary_sumstat(ld_data$G, info_snp, NCORES)
  beta_cols <- c("beta", "beta_se", "n_eff", "_NUM_ID_")
  well_behaved_snps <- (!qc$is_bad)
  df_beta <- info_snp[well_behaved_snps, beta_cols]

  # get ld matrix. Note, that we need 60gb of memory to keep 
  # all the hapmap variants in memory
  snp <- readRDS(args$path_ld_matrix)

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
    sid = pred_data$fam$sample.ID, 
    fid = pred_data$fam$sample.ID,
    pred_auto = final_pred_auto,
    pred_inf = final_red_inf
    )

  fwrite(resuts, file = paste0(out_prefix,".txt.gz"), sep = '\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bed_pred", default=NULL, required = TRUE, help = "Path for plink file (bed)")
parser$add_argument("--path_sumstat", default=NULL, required = TRUE, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--path_ld_matrix", default=NULL, required = TRUE, help = "Path to LD matrix (.rda / .rda file)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









