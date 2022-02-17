
# details: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html

# wrote a package that contain all dependencies for runnign LDpred2
# and combining various functions into easy to use pipelines. Note,
# that this will also load libraries, e.g. bigsnpr, bigassert
devtools::load_all('utils/modules/R/prstools')
library(argparse)


main <- function(args){

  print(args)
  stopifnot(file.exists(args$path_bed_ld))
  stopifnot(file.exists(args$path_sumstat))

  # setup parallel environment
  NCORES <- nb_cores()
  tmp <- tempfile(tmpdir = "data/tmp/tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  # load data required for setting up LD-matrix 
  ld_data <- load_bigsnp_from_bed(args$path_bed_ld)

  # load summary statistics
  sumstats <- read_hail_sumstat(args$path_sumstat)
  
  # match summary stats and LD data
  info_snp <- snp_match(sumstats, ld_data$map, join_by_pos = TRUE, strand_flip = FALSE)
 
  # QC summary statistics based on LD reference
  qc <- qc_binary_sumstat(ld_data$G, info_snp)
  beta_cols <- c("beta", "beta_se", "n_eff", "_NUM_ID_")
  well_behaved_snps <- (!qc$is_bad)
  df_beta <- info_snp[well_behaved_snps, beta_cols]

  # get ld matrix. Note, that we need 60gb of memory to keep 
  # all the hapmap variants in memory
  print("fitting ld matrix..")
  snp <- calc_ld_matrix(ld_data$G, ld_data$POS2, info_snp, chrs = 1:22, 
                        ncores = NCORES)

  outfile = paste0(args$out_prefix, ".rda")
  saveRDS(snp, file = outfile, compress = "xz")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_sumstat", default=NULL, required = TRUE, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--path_bed_ld", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









