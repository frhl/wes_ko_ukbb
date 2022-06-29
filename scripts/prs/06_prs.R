
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html


writelog <- function(msg, log) {
    time <- Sys.time()
    msg <- paste0(time,": ", msg, "\n")
    write(msg, log)
}


# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)

main <- function(args){

  stopifnot(file.exists(args$pred))
  stopifnot(file.exists(args$ldsc))
  stopifnot(!file.exists(args$tmp_bfile))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$chrom %in% paste0("chr",1:22))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  stopifnot(args$method %in% c('inf', 'auto'))

  log <- paste0(args$out_prefix,".log")

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)

  # laod ldsc and qced GWAS
  ldsc <- readRDS(args$ldsc)
  gwas <- ldsc$gwas
  h2 <- ldsc$coefficients$estimate[2]
  pvalue <- ldsc$coefficients$pvalue[2]

  # Estimate h2 chromosome-wide
  N_total <- nrow(gwas)
  N_chr <- sum(gwas$chr == args$chrom)
  h2_init <- h2 * (N_chr / N_total)
  
  # check estimates and ensure that
  # the trait is actually heritable 
  stopifnot(pvalue < 1e-5)
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
  bfile <- tempfile(tmpdir = dirname(args$out_prefix))
  pred <- load_bigsnp_from_bed(args$pred)
  pred <- match_bigsnp_with_gwas(obj=pred, gwas=gwas, bfile=args$tmp_bfile)
  genotypes <- pred$genotypes  
  indicies <- pred$gwas_indicies

  # check that we have genotypes
  stopifnot(!is.null(genotypes))
  stopifnot(!is.null(indicies))
  
  # dimensions  
  cols <- genotypes$`.->ncol` # variants
  rows <- genotypes$`.->nrow` # samples
  missing_gt <- NA

  # need to impute missing SNPs
  if (!is.null(args$impute)){
    sum_cols <- lapply(1:cols, function(i) return(sum(is.na(genotypes[,i])))) 
    missing_gt <- sum(unlist(sum_cols))
    genotypes <- snp_fastImputeSimple(genotypes, method = args$impute) 
  }

  # standardize genotypes 
  means <- NULL; sds <- NULL
  if (args$standardized_gt){
     means <- as.numeric(unlist(lapply(1:cols, function(i) mean(genotypes[,i], na.rm = TRUE))))
     sds <- as.numeric(unlist(lapply(1:cols, function(i) sd(genotypes[,i], na.rm = TRUE))))
  }

  if (args$method %in% "inf"){

    beta_inf <- snp_ldpred2_inf(snp$corr, gwas, h2_init)
    beta_inf <- beta_inf[indicies] 
    final_pred_inf <- big_prodVec(
          genotypes, 
          beta_inf,
          center = means,
          scale = sds,
          ncores = NCORES)  
    final_pred <- final_pred_inf

  } else if (args$method %in% "auto") {

     write("Running multi_auto", stderr())
     multi_auto <- snp_ldpred2_auto(
        corr = snp$corr,
        df_beta = gwas,
        h2_init = h2_init,
        vec_p_init = seq_log(1e-4, 0.9, NCORES), # using cores instead 30
        ncores = NCORES)

     # save data chains
     multi_auto_path <- paste0(args$out_prefix,'_chains.rda')
     saveRDS(multi_auto, multi_auto_path, compress = 'xz')

     # plot chains
     auto <- multi_auto[[1]]
     p <- plot_grid(
        qplot(y = auto$path_p_est) + 
          theme_bigstatsr() + 
          geom_hline(yintercept = auto$p_est, col = "blue") +
          scale_y_log10() +
          labs(y = "p"),
        qplot(y = auto$path_h2_est) + 
          theme_bigstatsr() + 
          geom_hline(yintercept = auto$h2_est, col = "blue") +
          labs(y = "h2"),
        ncol = 1, align = "hv"
      )
     # save plot
     outplot <-
     ggsave(plot=p,filename=paste0(args$outfile, "_chains.png"))

     # get estimates with indicies corresponding to pred genotypes
     beta_auto <- sapply(multi_auto, function(auto){
          auto$beta_est})

     # perform matrix multiplication
     pred_auto <- big_prodMat(
        genotypes,
        beta_auto,
        center = means,
        scale = sds,
        ncores = NCORES)

     # quality controls on chains
     sc <- apply(pred_auto, 2, sd, na.rm = TRUE)
     keep <- abs(sc - median(sc)) < 3 * mad(sc)
     final_beta_auto <- rowMeans(beta_auto[, keep], na.rm = TRUE) 
     
     # get final predicton
     final_pred_auto <- big_prodVec(
       genotypes, 
       final_beta_auto,
       center = means,
       scale = sds,
       ncores = NCORES)

     # keep 
     final_pred <- final_pred_auto
     beta_out <- cbind(gwas, final_beta_auto)
     fwrite(beta_out, file = paste0(args$out_prefix,"_betas.txt.gz"), sep = '\t')

    }

  
  # save parameters 
  model <- data.table(
     method = args$method,
     n_samples = unique(gwas$n),
     n_eff = unique(gwas$n_eff),
     n_snps = nrow(gwas),
     ldsc_h2_genome_wide_est = h2,
     h2_init_est = h2_init,
     inflation = calc_inflation(gwas$P),
     standardized_gt = args$standardized_gt,
     gt_rows = rows,
     gt_cols = cols,
     missing_gt = missing_gt
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
parser$add_argument("--standardized_gt", default=1, required = FALSE, help = "Should genotypes be standardized?")
parser$add_argument("--tmp_bfile", default=NULL, required = TRUE, help = "File path to temporary backing files")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--impute", default=NULL, required = TRUE, help = "Should missing genotypes be imputed? (See https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









