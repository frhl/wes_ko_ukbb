# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)
library(ggplot2)
library(tictoc)

main <- function(args){

  stopifnot(file.exists(args$ldsc))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  stopifnot(args$method %in% c('auto'))

  tic(paste0("LDPred2 ", basename(args$out_prefix)))
  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)

  # laod ldsc and qced GWAS
  ldsc <- readRDS(args$ldsc)
  h2_init <- ldsc$coefficients$estimate[2]
  pvalue <- ldsc$coefficients$pvalue[2]
  n_eff <- min(unique(ldsc$gwas$n_eff))
  gwas <- ldsc$gwas

  # cutoff prs if z-score estimate is not good enough
  if (!is.null(args$ldsc_n_eff_cutoff)){
    ldsc_n_eff_cutoff <- as.numeric(args$ldsc_n_eff_cutoff)
    if (n_eff < ldsc_n_eff_cutoff) stop("phenotype does not pass N_eff cutoff")
  }

  # cutoff prs if z-score estimate is not good enough
  if (!is.null(args$ldsc_pvalue_cutoff)){
    pvalue_cutoff <- as.numeric(args$ldsc_pvalue_cutoff)
    if (pvalue > pvalue_cutoff) stop("phenotype does not pass zstat cutoff")
  }

  # get SNP correlations and LD
  stopifnot(!is.null(gwas)) 
  snp <- get_ld_matrix(gwas)
  gwas <- gwas[na.omit(match(snp$map$marker, gwas$marker)), ]

  # check that LD-matrix markers and gwas markers have overlap
  # Check that ordering of markers are actually matching
  stopifnot(all(gwas$marker %in% snp$map$marker))
  stopifnot(all(snp$map$marker %in% gwas$marker)) 
  stopifnot(sum(gwas$marker == snp$map$marker) / nrow(gwas) == 1)
 
  # write SNP order
  marker_order <- gwas[,c("chr","pos","a0","a1","marker")]
  marker_path <- paste0(args$out_prefix,'.txt.gz')
  fwrite(marker_order, marker_path, sep="\t")

  ldpred_with_params <- function(num_iter, burn_in, h2_init, vec_p_init, seed=NULL) { 
     if (!is.null(seed)) set.seed(seed)
     multi_auto <- snp_ldpred2_auto(
         corr = snp$corr,
         df_beta = gwas,
         num_iter = num_iter,
         burn_in = burn_in,
         h2_init = h2_init,
         vec_p_init = vec_p_init,
         ncores = NCORES)
      beta_auto <- sapply(multi_auto, function(auto){auto$beta_est})
      chains_converged <- which(colSums(is.na(beta_auto))==0)
      converged <- length(chains_converged) >= 10 #  want at least 10 chains
      return(list(beta_auto=beta_auto, multi_auto=multi_auto, chains_converged=chains_converged, converged=converged))       
  }

  # use proper vec p init
  vec_p_ranges <- max(NCORES, as.numeric(args$vec_p_init_n))
  if (vec_p_ranges < 10) stop("You need at least 10 chains! set vec_p_init_n>=10.")
 
  # try the following combination of paramters in case of instability
  grid_params <- list(
     list(iter=500, burn_in=200, h2_init=h2_init, vec_p_init=seq_log(1e-4, 0.50, length.out=vec_p_ranges), seed=1),
     list(iter=500, burn_in=200, h2_init=h2_init, vec_p_init=seq_log(1e-4, 0.20, length.out=vec_p_ranges), seed=2),
     list(iter=250, burn_in=100, h2_init=h2_init, vec_p_init=seq_log(1e-4, 0.20, length.out=vec_p_ranges), seed=2)
   )
 
  step <- 0
  for (par in grid_params) {
     step <- step + 1
     write(paste("Beginning step", step, "for", args$out_prefix), stderr())
     ldpred_result <- ldpred_with_params(
         num_iter=par$iter,
         burn_in=par$burn_in,
         h2_init=par$h2_init,
         vec_p_init=par$vec_p_init,
         seed=par$seed
     )
     converged <- ldpred_result$converged
     multi_auto <- ldpred_result$multi_auto
     beta_auto <- ldpred_result$beta_auto
     if (converged) break 
  }
    
  # ensure that no missing values are present,
  msg <- paste0("Error: Stopping since all chains are NA: ", args$out_prefix)
  if (length(converged) == 0) stop(msg)
  write(paste0("Success +10 chain(s) passed for", args$out_prefix), stderr())

  # Chain QC
  (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
  (keep <- (range > (0.95 * quantile(range, 0.95, na.rm=TRUE))))
  keep[is.na(keep)] <- FALSE
  stopifnot(sum(keep) > 0)

  # save weights
  multi_auto_path <- paste0(args$out_prefix,'.rda')
  saveRDS(multi_auto, multi_auto_path, compress = 'xz')
  toc()
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--method", default=NULL, required = TRUE, help = "either 'inf' or 'auto'")
parser$add_argument("--ldsc", default=NULL, required = TRUE, help = ".rds object containing QCed GWAS and ldsc heritability estimates")
parser$add_argument("--ldsc_pvalue_cutoff", default=NULL, help = "cancel the run if the ldsc heritability p-value is not below the given treshold.")
parser$add_argument("--ldsc_n_eff_cutoff", default=NULL, help = "cancel the run if the ldsc N_eff is not below the given treshold.")
parser$add_argument("--vec_p_init_n", default=40, required = FALSE, help = "number of intial estimates to sample form (should be at least 5)")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)



