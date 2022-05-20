
library(argparse)
library(data.table)

main <- function(args){
 
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  files <- list.files(args$in_dir, full.names = TRUE, pattern = 'ldsc.+\\.rds')

  d <- do.call(rbind, lapply(files[1:3], function(rds){
    phenotype <- gsub("ldsc_","",tools::file_path_sans_ext(basename(rds)))
    d <- readRDS(rds)
    qc <- d$qc
    ldsc <- d$coefficients
    ldsc$n_snps <- qc$well_behaved_snps
    ldsc$phenotype <- phenotype
    return(ldsc)
  }))
 
  outfile <- paste0(args$out_prefix, ".txt")
  fwrite(d, outfile, sep = "\t") 
 
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









