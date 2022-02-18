
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
# issues with parallel: https://github.com/privefl/bigstatsr/issues/90

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)


main <- function(args){

  print(args)
  stopifnot(file.exists(args$bed))
  stopifnot(file.exists(args$gwas))

  # setup parallel environment
  NCORES <- max(1, nb_cores())
  bigparallelr::assert_cores(NCORES)
  
  # load data required for setting up LD-matrix 
  ld_data <- load_bigsnp_from_bed(args$bed)

  # load summary statistics
  sumstats <- read_hail_sumstat(args$gwas)
  
  # match summary stats and LD data
  info_snp <- snp_match(sumstats, ld_data$map, join_by_pos = TRUE, strand_flip = FALSE)
   
  #QC summary statistics based on LD reference
  qc <- qc_binary_sumstat(ld_data$G, info_snp, NCORES)
  well_behaved_snps <- (!qc$is_bad)
  df_beta <- info_snp[well_behaved_snps, ] 
  outfile_beta <- paste0(args$out_prefix, "_betas.txt.gz")
  fwrite(df_beta, outfile_beta, sep = '\t')

  # get ld matrix. Note, that we need 60gb of memory to keep 
  # all the hapmap variants in memory
  write("Fitting ld matrix..", stdout())
  snp_corr <- calc_ld_matrix(
                  ld_data$G, 
                  ld_data$POS2, 
                  df_beta, 
                  chrs = 1:22, 
                  ncores = NCORES, 
                  sfbm_file = args$out_prefix)

  # save with link to sfbm file
  outfile = paste0(args$out_prefix, ".rda")
  saveRDS(snp_corr, file = outfile, compress = "xz")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--gwas", default=NULL, required = TRUE, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--bed", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


