
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  print(args)
  stopifnot(file.exists(args$pred))
  stopifnot(file.exists(args$ld_matrix))
  stopifnot(file.exists(args$gwas)) 

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  
  # Load LD matrix and qced betas
  snp <- readRDS(args$ld_matrix)
  gwas <- fread(args$gwas)
  
  # check that LD and GWAS data was generated simultanously
  stopifnot(length(snp$ld) == nrow(gwas))
  
  # LD-score gives an estimate of h2 for a SNP which we can multiply 
  # by the total number of SNPs. Thus we are better off amalgamating all
  # the SNPs and running that through ldsc and thus avoiding extratpolating
  # that the average h2 explained by snps on chrN extends to all chromosomes 

  # perform LD regressio
  chi2 <- (gwas$beta / gwas$beta_se)^2
  ldsc <- with(gwas, 
               snp_ldsc(
                   snp$ld, 
                   length(snp$ld), 
                   chi2 = chi2,
                   sample_size = gwas$n_eff, 
                   blocks = NULL)
               )  
  
  ldsc_h2_est <- ldsc[["h2"]]
  ldsc_int_est <- ldsc[['int']] 

  # Load new data for prediction
  pred <- load_bigsnp_from_bed(args$pred)
  genotypes <- pred$G

  # ldpred2 inf model
  beta_inf <- snp_ldpred2_inf(snp$corr, gwas, ldsc_h2_est)
  final_pred_inf <- predict_prs(
       obj = pred,
       gwas = gwas,
       effects = beta_inf,
       ncores = NCORES)

  #final_pred_inf <- big_prodVec(
  #   genotypes, 
  #   beta_inf, 
  #   ind.col = gwas$`_NUM_ID_`)  
  
  # ldpred2 auto model
  #multi_auto <- snp_ldpred2_auto(
  #   snp$corr, 
  #   gwas, 
  #   h2_init = ldsc_h2_est,
  #   vec_p_init = seq_log(1e-4, 0.5, 30),
  #   ncores = NCORES)  
  
  # run first estimates
  #beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  #pred_auto <- big_prodMat(
  #    genotypes, 
  #    beta_auto, 
  #    ind.col = gwas[["_NUM_ID_"]],
  #    ncores = NCORES)
  
  # quality controls on chains
  #sc <- apply(pred_auto, 2, sd)
  #keep <- abs(sc - median(sc)) < 3 * mad(sc)
  #final_beta_auto <- rowMeans(beta_auto[, keep])

  # run final estimates
    
  #final_pred_auto <- big_prodVec(
  #    genotypes, final_beta_auto,
  #    ind.col = gwas[["_NUM_ID_"]],
  #    ncores = NCORES)

  # save parameters 
  model <- data.table(
     n_samples = unique(gwas$n),
     n_snps = nrow(gwas),
     ldsc_h2_est = ldsc_h2_est,
     ldsc_int_est = ldsc_int_est,
     inflation = calc_inflation(gwas$P)
     )
  
  # save results
  PGS <- data.table(
    sid = final_pred_inf$sid, 
    prs = final_pred_inf$prs
    )

  fwrite(PGS, file = paste0(args$out_prefix,".txt.gz"), sep = '\t')
  fwrite(model, file = paste0(args$out_prefix,".model"), sep = '\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--gwas", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--pred", default=NULL, required = TRUE, help = "Path to plink (bed) for PGS prediction")
parser$add_argument("--ld_matrix", default=NULL, required = TRUE, help = "Path to LD matrix (.rda / .rda file)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









