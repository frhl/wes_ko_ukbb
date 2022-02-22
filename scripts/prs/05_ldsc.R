
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)

main <- function(args){

  stopifnot(file.exists(args$gwas)) 
  stopifnot(file.exists(args$ld_bed))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$trait %in% c('binary', 'cts'))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  
  # Load LD matrix and summary statistics
  gwas <- read_hail_sumstat(args$gwas, trait = args$trait)
  ld_data <- load_bigsnp_from_bed(args$ld_bed)
  
  # match summary stats and LD data
  info_snp <- snp_match(gwas, ld_data$map, join_by_pos = TRUE, strand_flip = FALSE)

  # qc summary statistics
  qc <- qc_binary_sumstat(ld_data$G, info_snp, NCORES)
  well_behaved_snps <- (!qc$is_bad)
  gwas <- info_snp[well_behaved_snps, ]
  gwas$marker <- get_ldpred_marker(gwas)

  # Get LD matrix for final SNPs
  snp <- get_ld_matrix(gwas, ld_dir = args$ld_dir, verbose = TRUE)
 
  # match GWAS with snp-map
  indicies <- na.omit(match(snp$map$marker, gwas$marker))
  gwas <- gwas[indicies,]

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
  
  # what SNPS are used
  ldsc_out <- list(  
    h2_est = ldsc[["h2"]],
    int_est = ldsc[['int']],
    ld_map = snp$map,
    gwas = gwas
  ) 

  #write(paste0(args$pred, ".. done! Writing to ", args$out_prefix, ".rds"), stdout())
  saveRDS(ldsc_out, paste0(args$out_prefix,".rds"))
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--trait", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--gwas", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--ld_bed", default=NULL, required = TRUE, help = "Path to plink file (bed) used to design LD-matrix")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









