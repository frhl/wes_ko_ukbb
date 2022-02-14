

library(bigsnpr)
library(data.table)
library(argparse)
library(fmsb)
#options(bigstatsr.check.parallel.blas = FALSE)
#options(default.nproc.blas = NULL)

# load helpers (some functions required an active internet connection, 
# these have been re-written to rely on locally available data)
source('utils/modules/R/bigsnpr_helpers.R')


main <- function(x){

  stopifnot(file.exists(args$path_bed))
  stopifnot(file.exists(args$path_snp_match))
  stopifnot(file.exists(args$path_snp_match))

  ## copied from https://choishingwan.github.io/PRS-Tutorial/ldpred/
  NCORES <- nb_cores()
  tmp <- tempfile(tmpdir = "data/tmp/tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

  # Initialize variables for storing the LD score and LD matrix
  corr <- NULL
  ld <- NULL

  # We want to know the ordering of samples in the bed file 
  fam.order <- NULL

  # preprocess the bed file (only need to do once for each data set) and attach 
  snp_readBed(path_bed)
  basename <- tools::file_path_sans_ext(path_bed)
  rds <- paste0(basename,'.rds')
  obj.bigSNP <- snp_attach(rds)

  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")

  # read summary statistics file for SNP matching
  sumstats <- read_hail_sumtat(path_bed)

  # perform SNP matching
  info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)

  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes

  # Rename the data structures
  CHR <- as.numeric(gsub('chr','',map$chr))
  POS <- map$pos
  invalid_CHR <- unique(map$chr[is.na(CHR)])
  if (length(invalid_CHR) > 0){
    write(paste("Chromosome", invalid_CHR[1], "is not among HAPMAP SNPs. Exiting.."), stderr())
  }

  # get the CM information from hapmap SNPs
  POS2 <- snp_asGeneticPosLocal(CHR, POS, mapdir = "data/prs/1000-genomes-genetic-maps",genetic_map = 'hapmap')

  # check for sufficient overlap
  matches <- sum(POS2==0)/length(POS2) # hapmap has many less missing variants than omni
  write(paste0(1-round(matches*100,2),'% of LD panel variants are in genetic map (hapmap).'),stdout())
   
  # calculate LD panel
  chrs <- paste0("chr",1:22)
  for (chr in chrs) {
      # Extract SNPs that are included in the chromosome
      ind.chr <- which(info_snp$chr == chr)
      ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
      stopifnot(length(ind.chr) > 0)
      # Calculate the LD
      corr0 <- snp_cor(
              genotype,
              ind.col = ind.chr2,
              ncores = NCORES,
              infos.pos = POS2[ind.chr2],
              size = 3 / 1000
          )
      if (chr == 1) {
          ld <- Matrix::colSums(corr0^2)
          corr <- as_SFBM(corr0, tmp)
      } else {
          ld <- c(ld, Matrix::colSums(corr0^2))
          corr$add_columns(corr0, nrow(corr))
      }
  }
  # We assume the fam order is the same across different chromosomes
  fam.order <- as.data.table(obj.bigSNP$fam)
  # Rename fam order
  setnames(fam.order,
          c("family.ID", "sample.ID"),
          c("FID", "IID"))

  # perform LD-score regression
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  ldsc <- snp_ldsc(   ld, 
                      length(ld), 
                      chi2 = (df_beta$beta / df_beta$beta_se)^2,
                      sample_size = df_beta$n_eff, 
                      blocks = NULL)
  h2_est <- ldsc[["h2"]]## 4. Perform LD score regression

  # calculate NULL R2

  # Obtain model PRS

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bed", default=NULL, help = "Path for plink file (bed)")
parser$add_argument("--path_snp_match", default=NULL, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--path_ld_bed", default=NULL, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
args <- parser$parse_args()

main()









