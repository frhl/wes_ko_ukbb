
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  stopifnot(file.exists(args$pred))
  stopifnot(file.exists(args$gwas)) 
  stopifnot(file.exists(args$ld_bed))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$method %in% c('inf', 'auto'))
  stopifnot(args$trait %in% c('binary', 'cts'))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  
  # Load LD matrix and summary statistics
  gwas <- read_hail_sumstat(args$gwas, trait = args$trait)
  ld_data <- load_bigsnp_from_bed(args$ld_bed)
  pred <- load_bigsnp_from_bed(args$pred)
  
  # match summary stats and LD data
  info_snp <- snp_match(gwas, ld_data$map, join_by_pos = TRUE, strand_flip = FALSE)

  # qc summary statistics
  qc <- qc_binary_sumstat(ld_data$G, info_snp, NCORES)
  well_behaved_snps <- (!qc$is_bad)
  gwas <- info_snp[well_behaved_snps, ]
  gwas$marker <- get_ldpred_marker(gwas)

  # Get LD matrix for final SNPs
  snp <- get_single_ld_matrix(gwas, chr = args$chrom, ld_dir = args$ld_dir)
 
  # match GWAS with snp-map
  indicies <- na.omit(match(snp$map$marker, gwas$marker))
  gwas <- gwas[indicies,]

  print(nrow(gwas))

  # check that LD-matrix markers and gwas markers have overlap
  # Check that ordering of markers are actually matching
  stopifnot(all(gwas$marker %in% snp$map$marker))
  stopifnot(all(snp$map$marker %in% gwas$marker)) 
  stopifnot(sum(gwas$marker == snp$map$marker) / nrow(gwas) == 1)

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

  if (args$method %in% "inf"){
  
    beta_inf <- snp_ldpred2_inf(snp$corr, gwas, ldsc_h2_est)
    final_pred_inf <- predict_prs(
       obj = pred,
       gwas = gwas,
       effects = beta_inf,
       ncores = NCORES)

    final_pred <- final_pred_inf

  } else if (args$method %in% "auto") {

   write("Running multi_auto", stderr())
   multi_auto <- snp_ldpred2_auto(
      snp$corr,
      gwas,
      h2_init = ldsc_h2_est,
      vec_p_init = seq_log(1e-4, 0.5, 30),
      ncores = NCORES)
    
   print(str(multi_auto))

   # match prediction genotypes with gwas
   write("matching pred & gwas..", stderr())
   pred_match <- match_bigsnp_with_gwas(pred, gwas)
   
   # get estimates with indicies corresponding to pred genotypes
   write("Get estimates for beta_auto..", stderr())
   beta_auto <- sapply(multi_auto, function(auto){
        auto$beta_est[pred_match$gwas_indicies]})
   
   # perform matrix multiplication
   write("Performing mat mul..", stderr())
   pred_auto <- big_prodMat(
      pred_match$genotypes,
      beta_auto,
      ncores = NCORES)

   # quality controls on chains
   write("Quality control om chains..", stderr())
   sc <- apply(pred_auto, 2, sd)
   keep <- abs(sc - median(sc)) < 3 * mad(sc)
   final_beta_auto <- rowMeans(beta_auto[, keep]) 
   
   # get final predicton
   write("Getting final prediction..", stderr())
   final_pred_auto <- big_prodVec(
     pred_match$genotypes, 
     final_beta_auto,
     ncores = NCORES)

   final_pred <- final_pred_auto
   final_pred$sid <- pred_match$sid
  }

  # save parameters 
  model <- data.table(
     method = args$method,
     n_samples = unique(gwas$n),
     n_eff = unique(gwas$n_eff),
     n_snps = nrow(gwas),
     ldsc_h2_est = ldsc_h2_est,
     ldsc_int_est = ldsc_int_est,
     inflation = calc_inflation(gwas$P)
     )
  
  # save results
  PGS <- data.table(
    sid = final_pred$sid, 
    prs = final_pred$prs
    )

  write(paste0(args$pred, ".. done! Writing to ", args$out_prefix, ".txt.gz"), stdout())
  fwrite(PGS, file = paste0(args$out_prefix,".txt.gz"), sep = '\t')
  fwrite(model, file = paste0(args$out_prefix,".model"), sep = '\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--method", default=NULL, required = TRUE, help = "either 'inf' or 'auto'")
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "what chromosome")
parser$add_argument("--trait", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--gwas", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--pred", default=NULL, required = TRUE, help = "Path to plink (bed) for PGS prediction")
parser$add_argument("--ld_bed", default=NULL, required = TRUE, help = "Path to plink file (bed) used to design LD-matrix")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









