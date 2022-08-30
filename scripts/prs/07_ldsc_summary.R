
library(argparse)
library(data.table)

main <- function(args){
 
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  files <- list.files(args$in_dir, full.names = TRUE, pattern = 'ldsc.+\\.rds')

  d <- do.call(rbind, lapply(files, function(rds){
    phenotype <- gsub("ldsc_","",tools::file_path_sans_ext(basename(rds)))
    d <- readRDS(rds)
    d$log_pvalue <- NULL
    qc <- d$qc
    ldsc <- d$coefficients
    gwas <- d$gwas
    ldsc$n_snps <- ifelse(qc$disable_qc, qc$total_snps, qc$well_behaved_snps)
    ldsc$phenotype <- phenotype
    ldsc$coef <- rownames(ldsc)
    ldsc$n <- gwas$n[1]
    ldsc$n_eff <- gwas$n_eff[1]
    return(ldsc)
  }))

  # write outfile
  outfile1 <- paste0(args$out_prefix, ".txt.gz")
  fwrite(d, outfilei1, sep = "\t") 

  # save list of phenotypes that pass our sig thresholds
  pval_threshold <- 1e-5
  phenotypes <- d$phenotype[d$coef == "h2" & d$pvalue > pval_threshold]
  outfile2 <- paste0(args$out_prefix, "_qc_passed.txt")
  fwrite(data.frame(phenotypes), outfile2, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









