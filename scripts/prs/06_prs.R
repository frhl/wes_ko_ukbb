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

  # cutoff prs if z-score estimate is not good enough
  if (!is.null(args$ldsc_pvalue_cutoff)){
    pvalue_cutoff <- as.numeric(args$ldsc_pvalue_cutoff)
    if (pvalue > pvalue_cutoff) stop("phenotype does not pass threshold")
  }

  # Estimate h2 chromosome-wide
  N_total <- nrow(gwas)
  N_chr <- sum(gwas$chr == args$chrom)
  h2_init <- h2 * (N_chr / N_total)
  
  # check if estimates are ok
  if (N_total < 1) stop(paste0("No rows in GWAS: ", args$ldsc))
  if (N_chr < 1) stop(paste0("No matching chromsomes in gwas rows: ", args$ldsc))
  if (N_chr < 10000) stop(paste0("Less than 10K variants in GWAS rows: ", args$ldsc)) 
  if (h2_init <= 0) stop(paste0("h2_init (",h2_init ,") is zero or negative:",  args$ldsc))

  # check estimates and ensure that
  # the trait is actually heritable 
  #stopifnot(pvalue < 1e-5)
  stopifnot(!is.null(gwas)) 

  # load prediction file (note, that this needs to be done before loading
  # LD-matrix, so that it can be subsetted accordinly).
  pred <- load_bigsnp_from_bed(args$pred)
  gwas <- gwas[gwas$rsid %in% pred$map$rsid,]

  # need to keep track of matching SNPs.
  n_pred_in_gwas <- sum(pred$map$rsid %in% gwas$rsid)
  total <- nrow(pred$map)
  pct <- paste0("(",round((n_pred_in_gwas/total)*100, 2),"%)")
  msg <- paste("Only",n_pred_in_gwas,"of",total,"SNPs", pct,"GWAS and prediction SNVs:", args$out_prefix)
  if (n_pred_in_gwas < 1000) stop(msg)
      
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
  pred <- match_bigsnp_with_gwas(obj=pred, gwas=gwas, bfile=args$tmp_bfile)
  genotypes <- pred$genotypes  
  indicies <- pred$gwas_indicies

  # check that we have genotypes
  stopifnot(!is.null(genotypes))
  stopifnot(!is.null(indicies))
    
  # check that things are mathced
  stopifnot(pred$map$rsid == gwas$rsid)

  # dimensions  
  cols <- genotypes$`.->ncol` # variants
  rows <- genotypes$`.->nrow` # samples
  missing_gt <- NA

  # need to impute missing SNPs
  if (!is.null(args$impute)){
    write(paste0("Imputing genotypes in ", args$out_prefix), stderr())
    sum_cols <- lapply(1:cols, function(i) return(sum(is.na(genotypes[,i])))) 
    missing_gt <- sum(unlist(sum_cols))
    genotypes <- snp_fastImputeSimple(genotypes, method = args$impute) 
  }

  # standardize genotypes 
  means <- NULL; sds <- NULL
  if (args$standardized_gt){
     write(paste0("Standardizing genotypes in ", args$out_prefix), stderr())
     means <- as.numeric(unlist(lapply(1:cols, function(i) mean(genotypes[,i], na.rm = TRUE))))
     sds <- as.numeric(unlist(lapply(1:cols, function(i) sd(genotypes[,i], na.rm = TRUE))))
     stopifnot(sum(is.na(means))==0)
     stopifnot(sum(is.na(sds))==0)
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

     # use proper vec p init
     vec_p_ranges <- max(NCORES, as.numeric(args$vec_p_init_n))
     if (vec_p_ranges < 10) stop("You need at least 10 chains! set vec_p_init_n>=10.")
     
     num_iter = 1000
     burn_in = 1000


     #num_iter=round(runif(1, 100, 200))
     #burn_in=round(runif(1, 100, 200))

     vec_p_init = seq(0.001, 0.60, length.out=vec_p_ranges)
     #vec_p_init = seq(0.001, 0.9, length.out=100)

     write(paste0("Starting LDPred2-auto for ",args$out_prefix,".."), stderr())
     #write(paste0("h2=",h2, "\nh2_chrom=",h2_init,"\nn_variants=",N_chr), stderr())
     #write(paste0("iter=",num_iter,"\nburn_in=",burn_in,"\nvec_p_ranges=",vec_p_ranges), stderr())
     multi_auto <- snp_ldpred2_auto(
        corr = snp$corr,
        df_beta = gwas,
        num_iter = num_iter,
        burn_in = burn_in,
        h2_init = h2_init,
        vec_p_init = vec_p_init,
        ncores = NCORES)
     
     beta_auto <- sapply(multi_auto, function(auto){auto$beta_est})
     converged <- which(colSums(is.na(beta_auto))==0)

     #vec_p_init = seq_log(p_min, 0.6, length.out=vec_p_ranges), # using cores instead 30
     #vec_p_init = seq_log(1e-4, 0.3, length.out=vec_p_ranges), # using cores instead 30

     # ensure that no missing values are present,
     msg <- paste0("Error: Stopping since all chains are NA: ", args$out_prefix)
     if (length(converged) == 0) stop(msg)
     write(paste0("Success! ",length(converged)," chain(s) passed for", args$out_prefix), stderr())

     # save chains
     multi_auto_path <- paste0(args$out_prefix,'_chains.rda')
     saveRDS(multi_auto, multi_auto_path, compress = 'xz')

     # plot example chain of one that is not NA
     converged_index <- converged[1]
     auto <- multi_auto[[converged_index]]
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
     ggsave(plot=p,filename=paste0(args$out_prefix, "_chains.png"))

     # check if initial estimate caused any problems.
     stopifnot(all(rowSums(!is.na(beta_auto)) > 0))
     na_cols <- colSums(is.na(beta_auto)) == nrow(beta_auto)
     n_na_cols <- sum(sum(na_cols))
     n_cols <- ncol(beta_auto)
     if (n_na_cols > 0){
        write(paste0("Note: ",n_na_cols," of ",n_cols," beta_auto NA cols will be discarded (", basename(args$out_prefix), ")"), stdout()) 
        beta_auto <- beta_auto[,!na_cols]
     } 

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
     stopifnot(!any(is.na(keep)))
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


  # count PRS that could not be generated
  missing_prs <- round(sum(is.na(final_pred))/length(final_pred)*100, 2)

  # save parameters 
  model <- data.table(
     ldsc = args$ldsc,
     pred = args$pred,
     chrom = args$chrom,
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
     missing_gt = missing_gt,
     missing_prs_pct = missing_prs
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
parser$add_argument("--ldsc_pvalue_cutoff", default=NULL, help = "cancel the run if the ldsc heritability p-value is not below the given treshold.")
parser$add_argument("--standardized_gt", default=1, required = FALSE, help = "Should genotypes be standardized?")
parser$add_argument("--vec_p_init_n", default=200, required = FALSE, help = "number of intial estimates to sample form (should be at least 5)")
parser$add_argument("--tmp_bfile", default=NULL, required = TRUE, help = "File path to temporary backing files")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--impute", default=NULL, required = TRUE, help = "Should missing genotypes be imputed? (See https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









