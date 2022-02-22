
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  stopifnot(file.exists(args$pred))
  stopifnot(file.exists(args$ldsc))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$method %in% c('inf', 'auto'))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  msg <- paste0("Running ldpred2 with ",NCORES," cores [",Sys.time(),"]") 
  write(msg,stdout())

  # laod ldsc and qced GWAS
  ldsc <- readRDS(args$ldsc)
  gwas <- ldsc$gwas
  h2 <- ldsc$h2_est

  # Estimate h2 chromosome-wide
  N_total <- nrow(gwas)
  N_chr <- sum(gwas$chr == args$chr)
  h2_init <- h2 * (N_chr / N_total)
  
  # check estimates
  stopifnot(h2_init > 0)
  stopifnot(!is.null(gwas)) 

  # get SNP correlations and LD
  snp <- get_single_ld_matrix(gwas, chr = args$chrom, ld_dir = args$ld_dir)
 
  # match GWAS with snp-correlation map
  gwas <- gwas[na.omit(match(snp$map$marker, gwas$marker)), ]

  # check that LD-matrix markers and gwas markers have overlap
  # Check that ordering of markers are actually matching
  stopifnot(all(gwas$marker %in% snp$map$marker))
  stopifnot(all(snp$map$marker %in% gwas$marker)) 
  stopifnot(sum(gwas$marker == snp$map$marker) / nrow(gwas) == 1)
  
   # load data to be used for prediction 
  pred <- load_bigsnp_from_bed(args$pred)
  pred <- match_bigsnp_with_gwas(pred, gwas)
  genotypes <- pred$genotypes  
  indicies <- pred$gwas_indicies

  stopifnot(!is.null(genotypes))
  stopifnot(!is.null(indicies))

  if (args$method %in% "inf"){

    beta_inf <- snp_ldpred2_inf(snp$corr, gwas, h2_init)
    beta_inf <- beta_inf[indicies] 
    final_pred_inf <- big_prodVec(
          genotypes, 
          beta_inf,
          ncores = NCORES)  
    final_pred <- final_pred_inf

  } else if (args$method %in% "auto") {

     write("Running multi_auto", stderr())
     multi_auto <- snp_ldpred2_auto(
        corr = snp$corr,
        df_beta = gwas,
        h2_init = h2_init,
        vec_p_init = seq_log(1e-4, 0.5, 30),
        ncores = NCORES)

     # save data chains
     multi_auto_path <- paste0(args$out_prefix,'_chains.rda')
     saveRDS(multi_auto, multi_auto_path, compress = 'xz')

     # get estimates with indicies corresponding to pred genotypes
     beta_auto <- sapply(multi_auto, function(auto){
          auto$beta_est})
     
     # perform matrix multiplication
     pred_auto <- big_prodMat(
        genotypes,
        beta_auto,
        ncores = NCORES)

     # quality controls on chains
     na_rows <- rowSums(is.na(pred_auto)) > 0 
     if (any(na_rows)) {
        na_rows_pct <- round(100*(na_rows / nrow(pred_auto)), 2)
        write(paste0(na_rows_pct,"% of chain rows contains NAs"), stderr())
     }
     
     sc <- apply(pred_auto, 2, sd, na.rm = TRUE)
     keep <- abs(sc - median(sc)) < 3 * mad(sc)
     final_beta_auto <- rowMeans(beta_auto[, keep], na.rm = TRUE) 
     
     # get final predicton
     final_pred_auto <- big_prodVec(
       genotypes, 
       final_beta_auto,
       ncores = NCORES)

     final_pred <- final_pred_auto
    }

  
  # save parameters 
  model <- data.table(
     method = args$method,
     n_samples = unique(gwas$n),
     n_eff = unique(gwas$n_eff),
     n_snps = nrow(gwas),
     ldsc_h2_genome_wide_est = h2,
     h2_init_est = h2_init,
     inflation = calc_inflation(gwas$P)
     )
  
  # save results
  PGS <- data.table(
    sid = pred$sid, 
    prs = final_pred
    )

  write(paste0(args$pred, ".. done! Writing to ", args$out_prefix, ".txt.gz"), stdout())
  fwrite(PGS, file = paste0(args$out_prefix,".txt.gz"), sep = '\t')
  fwrite(model, file = paste0(args$out_prefix,".model"), sep = '\t')
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "chromosome input")
parser$add_argument("--method", default=NULL, required = TRUE, help = "either 'inf' or 'auto'")
parser$add_argument("--pred", default=NULL, required = TRUE, help = "Path to plink (bed) for PGS prediction")
parser$add_argument("--ldsc", default=NULL, required = TRUE, help = ".rds object containing QCed GWAS and ldsc heritability estimates")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









