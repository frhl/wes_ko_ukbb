
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)


get_sd_y <- function(path_cts_phenotypes, phenotype){
    stopifnot(file.exists(path_cts_phenotypes))
    phenos <- fread(path_cts_phenotypes)
    stopifnot(phenotype %in% colnames(phenos))
    y <- phenos[[phenotype]]
    return(sd(y, na.rm = TRUE))
}

main <- function(args){

  print(args)
  stopifnot(file.exists(args$gwas)) 
  stopifnot(file.exists(args$ld_bed))
  stopifnot(dir.exists(args$ld_dir))
  stopifnot(args$out_prefix != "")
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
  if (args$trait %in% "binary"){
    qc <- qc_sumstat_binary(ld_data$G, info_snp, ncores = NCORES)
  } else {
    # Note, the cts traits require the standard deviation of the phenotype
    sd_y <- get_sd_y(args$path_cts_phenotypes, args$phenotype)
    qc <- qc_sumstat_cts(ld_data$G, info_snp, sd_y, ncores = NCORES)
  }  
  
  # what SNPs should be kept?
  well_behaved_snps <- (!qc$is_bad)
  keep_snps <- well_behaved_snps
  if (args$disable_qc) keep_snps <- TRUE
  gwas <- info_snp[keep_snps, ]
  gwas$marker <- get_ldpred_marker(gwas)

  # get qc data.frame
  d_qc <- data.frame(
    well_behaved_snps = sum(well_behaved_snps), 
    total_snps = length(well_behaved_snps),
    disable_qc = args$disable_qc
  )

  # Get LD matrix for final SNPs
  snp <- get_ld_matrix(gwas, chrs = 1:22, ld_dir = args$ld_dir, verbose = TRUE)
 
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
                   blocks = 200, # if NULL, then SE's are not estimated
                   ncores = NCORES,
                   )
               )  

  # extract coeffecients
  int_est <- ldsc[["int"]]
  h2_est <- ldsc[["h2"]]
  h2_se <- ldsc[["h2_se"]]
  int_se <- ldsc[["int_se"]]

  # measures the proportion of the inflation in the mean chi^2 that 
  # the LD Score regression intercept ascribes to causes other than 
  # polygenic heritability
  ratio <-  (int_est-1)/(mean(chi2, na.rm=TRUE)-1)

  # calculate P-values for linear fit
  h2_z <- h2_est / h2_se
  int_z <- int_est / int_se
  h2_pval <- 2 * pnorm(abs(h2_z), lower.tail = FALSE)
  int_pval <- 2 * pnorm(abs(int_z), lower.tail = FALSE)

  # organize in table
  coefficients <- data.frame(
    estimate = c(int_est, h2_est),
    std_error = c(int_se, h2_se),
    zstat = c(int_z, h2_z),
    pvalue = c(int_pval, h2_pval),
    ratio = ratio
  )

  # rename table
  colnames(coefficients) <- c("estimate", "std_error", "zstat", "pvalue", "ratio")
  rownames(coefficients) <- c("intercept", "h2")

  # what SNPS are used
  ldsc_out <- list(  
    coefficients = coefficients,
    qc = d_qc,
    ldsc = ldsc,
    ld_map = snp$map,
    gwas = gwas
  ) 

  #write(paste0(args$pred, ".. done! Writing to ", args$out_prefix, ".rds"), stdout())
  print(args)
  outfile = paste0(args$out_prefix, ".rds")
  write(paste("Done! writing to", outfile), stdout())
  saveRDS(ldsc_out, outfile)
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--trait", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--gwas", default=NULL, required = TRUE, help = "Path to QCed SNPs")
parser$add_argument("--ld_bed", default=NULL, required = TRUE, help = "Path to plink file (bed) used to design LD-matrix")
parser$add_argument("--ld_dir", default=NULL, required = TRUE, help = "Path to directory with pre-calcualted SNP correlations and LD (.rds files)")
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "String. Current phenotype")
parser$add_argument("--disable_qc", default = FALSE, action='store_true', required = FALSE, help = "String. Current phenotype")
parser$add_argument("--path_cts_phenotypes", default=NULL, required = TRUE, help = "Path to cts phenotypes")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









