

library(bigsnpr)
library(data.table)
library(argparse)
library(fmsb)
#options(bigstatsr.check.parallel.blas = FALSE)
#options(default.nproc.blas = NULL)

# load helpers (some functions required an active internet connection, 
# these have been re-written to rely on locally available data)
source('utils/modules/R/bigsnpr_helpers.R')


main <- function(args){

  stopifnot(file.exists(args$path_bed))
  stopifnot(file.exists(args$path_ld_bed))
  stopifnot(file.exists(args$path_sumstat))
  stopifnot(file.exists(args$path_pheno))

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
  if (!file.exists.ext(path_ld_bed, '.bk')) snp_readBed(path_ld_bed)
  basename <- tools::file_path_sans_ext(path_bed)
  rds <- paste0(basename,'.rds')
  obj.bigSNP <- snp_attach(rds)

  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")

  # remove non-autosomes (e.g. chr8_KI270821v1_alt)
  autosomes <- paste0('chr',1:22)
  obj.bigSNP$map <- obj.bigSNP$map[obj.bigSNP$map$chr %in% autosomes,]

  # read summary statistics file for SNP matching
  sumstats <- read_hail_sumtat(path_sumstat)

  # perform SNP matching
  info_snp <- bigsnpr::snp_match(sumstats, map, strand_flip = FALSE)

  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes

  # Rename the data structures
  CHR <- as.numeric(gsub('chr','',map$chr))
  POS <- map$pos

  # get the CM information from hapmap SNPs
  POS2 <- snp_asGeneticPosLocal(CHR, POS, mapdir = "data/prs/1000-genomes-genetic-maps",genetic_map = 'hapmap')

  # check for sufficient overlap
  matches <- sum(POS2==0)/length(POS2) # hapmap has many less missing variants than omni
  write(paste0(100*(1-round(matches,5)),'% of LD panel variants are in genetic map (hapmap).'),stdout())
   
  # calculate correlation between variants (LD panel)
  # Note: this function MUST start at chr1 (see first IF STATEMENT)
  # Note: Currently this only works for GRCh38. Changes "Chr" downstream.
  chrs <- paste0("chr",1:22)
  stopifnot(all(unique(info_snp$chr) %in% chrs))
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
      if (chr == "chr1") {
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
  h2_est <- ldsc[["h2"]]
  print(h2_est)

  # calculate NULL R2
  # Reformat the phenotype file such that y is of the same order as the 
  
  # sample ordering in the genotype file
  y <- pheno[fam.order, on = c("FID", "IID")]

  # Calculate the null R2
  null_model <- as.formula(paste0(args$response, "~", gsub(',','+', args$covars)))
  if (args$response_type == 'cts'){
    summary <- lm(null_model, data = y)
    null.r2 <- null.model$r.squared
  } else if (args$response_type == 'binary'){
    summary <- glm(null_model, data = y, family=binomial)
    null.r2 <- fmsb::NagelkerkeR2(null.model)
  } else {
    stop('response type must be either "binary" or "cts"')
  }

  # Read in full dataset (all samples) for PRS calculation
  pred_obj.bigSNP <- snp_attach(args$path_bed)
  pred_obj.bigSNP$map <- pred_obj.bigSNP$map[pred_obj.bigSNP$map$chr %in% autosomes,]
  genotype <- pred_obj.bigSNP$genotypes
  ind.test <- 1:nrow(genotype)

  # obtain (genome wide) model PRS
  if (args$model in 'auto'){
    
    # Get adjusted beta from the auto model
    multi_auto <- snp_ldpred2_auto(
        corr,
        df_beta,
        h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
        ncores = NCORES
    )
    beta_auto <- sapply(multi_auto, function(auto)  auto$beta_est)

    # calculate PRS for all samples
    pred_auto <-
        big_prodMat(genotype,
                    beta_auto,
                    ind.row = ind.test,
                    ind.col = info_snp$`_NUM_ID_`)
    # scale the PRS generated from AUTO
    pred_scaled <- apply(pred_auto, 2, sd)
    final_beta_auto <-
        rowMeans(beta_auto[,
                    abs(pred_scaled -
                        median(pred_scaled)) <
                        3 * mad(pred_scaled)])
    pred_auto <-
        big_prodVec(genotype,
            final_beta_auto,
            ind.row = ind.test,
            ind.col = info_snp$`_NUM_ID_`)
    
    pred_final <- pred_auto

    } else if (args$model in 'grid'){
      
      # Prepare data for grid model
      p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
      h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
      grid.param <-
          expand.grid(p = p_seq,
                  h2 = h2_seq,
                  sparse = c(FALSE, TRUE))
      
      # Get adjusted beta from grid model
      beta_grid <- snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)

      pred_grid <- big_prodMat(genotype, 
                            beta_grid, 
                            ind.col = info_snp$`_NUM_ID_`)
      
      pred_final <- pred_grid

    } else if (args$model in 'inf'){
      beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
      pred_inf <- big_prodVec(genotype,
                            beta_inf,
                            ind.row = ind.test,
                            ind.col = info_snp$`_NUM_ID_`)

      pred_final <- pred_inf

    } else {
     stop('current "model" must be either "auto", "grid" or "inf".')
   }
  
   print(str(pred_final))
   print(head(pred_final))       
 

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_bed", default=NULL, required = TRUE, help = "Path for plink file (bed)")
parser$add_argument("--path_sumstat", default=NULL, required = TRUE, help = "Path to a summary statistics file with matching SNPs")
parser$add_argument("--path_ld_bed", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--path_pheno", default=NULL, required = TRUE, help = "Path to files with phenotypes")
parser$add_argument("--covars", default=NULL, required = TRUE, help = "covariates as a string seperated by comma(s).")
parser$add_argument("--response", default=NULL, required = TRUE, help = "Name of response (which is also in phenotype file)")
parser$add_argument("--response_type", default=NULL, required = TRUE, help = "Either 'binary' or 'cts'")
parser$add_argument("--model", default=NULL, required = TRUE, help = "Must be 'auto', 'grid' or 'inf'.")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









