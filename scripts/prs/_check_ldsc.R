library(argparse)

main <- function(args){

  stopifnot(file.exists(args$ldsc))
    
  # laod ldsc and qced GWAS
  ldsc <- readRDS(args$ldsc)
  h2 <- ldsc$coefficients$estimate[2]
  pvalue <- ldsc$coefficients$pvalue[2]
  pvalue_cutoff <- as.numeric(args$ldsc_pvalue_cutoff)
  if (pvalue > pvalue_cutoff) stop(paste(args$ldsc, "does not pass P-value threshold."))
  if (h2 <= 0) stop(paste0("h2 (",h2 ,") is zero or negative:",  args$ldsc))
  write(paste0(ldsc, ".. ok!"), stderr())
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ldsc", default=NULL, required = TRUE, help = ".rds object containing QCed GWAS and ldsc heritability estimates")
parser$add_argument("--ldsc_pvalue_cutoff", default=NULL, help = "cancel the run if the ldsc heritability p-value is not below the given treshold.")
args <- parser$parse_args()

main(args)









