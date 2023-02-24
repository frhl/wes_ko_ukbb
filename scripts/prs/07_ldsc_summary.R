
library(argparse)
library(data.table)

main <- function(args){
 
  stopifnot(dir.exists(args$in_dir))
  stopifnot(dir.exists(dirname(args$out_prefix)))
  stopifnot(file.exists(args$phenotypes))
  files <- list.files(args$in_dir, full.names = TRUE, pattern = 'ldsc.+\\.rds')
  phenotypes <- readLines(args$phenotypes)

  d <- do.call(rbind, lapply(files, function(rds){
    phenotype <- gsub("ldsc_","",tools::file_path_sans_ext(basename(rds)))
    d <- readRDS(rds)
    qc <- d$qc
    ldsc <- d$coefficients
    gwas <- d$gwas
    ldsc$n_snps <- ifelse(qc$disable_qc, qc$total_snps, qc$well_behaved_snps)
    ldsc$phenotype <- phenotype
    ldsc$coef <- rownames(ldsc)
    ldsc$n <- gwas$n[1]
    ldsc$n_eff <- gwas$n_eff[1]
    ldsc$log_pvalue <- NULL
    return(ldsc)
  }))

  # subset to overlapping phenotypes
  prs_phenotypes <- unique(d$phenotype)
  final_phenotypes <- prs_phenotypes[prs_phenotypes %in% phenotypes]
  d <- d[d$phenotype %in% final_phenotypes,]

  # write outfile
  outfile1 <- paste0(args$out_prefix, ".txt.gz")
  fwrite(d, outfile1, sep = "\t") 

  # bonferroni correction
  n_phenos <- length(final_phenotypes)
  p_cutoff <- 0.05 / n_phenos

  # keep these phenotypes
  phenos_with_nom_sig_prs <- d$phenotype[(d$coef == "h2") & (d$pvalue < 0.05)]
  phenos_with_bonf_sig_prs <- d$phenotype[(d$coef == "h2") & (d$pvalue < p_cutoff)]

  # write bonf and nom sig prs phenos 
  outfile0 <- paste0(args$out_prefix, "_tested_phenos.txt")
  fwrite(data.frame(d$phenotype), outfile0, sep = '\t', col.names = FALSE)
  outfile1 <- paste0(args$out_prefix, "_nom_sig_phenos.txt")
  fwrite(data.frame(phenos_with_nom_sig_prs), outfile1, sep = '\t', col.names = FALSE)
  outfile2 <- paste0(args$out_prefix, "_bonf_sig_phenos.txt")
  fwrite(data.frame(phenos_with_bonf_sig_prs), outfile2, sep = '\t', col.names = FALSE)

 }

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "What phenotypes are being considered?")
args <- parser$parse_args()

main(args)









