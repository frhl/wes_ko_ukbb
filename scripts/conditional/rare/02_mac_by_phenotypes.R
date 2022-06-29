#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(parallel)
library(doParallel)


main <- function(args){


  stopifnot(file.exists(args$in_vcf)) 
  stopifnot(file.exists(args$phenotypes))
  stopifnot(dir.exists(dirname(args$out_prefix)))

  # read in VCF
  cmd <- paste0("zcat ", args$in_vcf, " | grep -v '##'")
  d <- fread(cmd = cmd)

  # get columns with genotypes / metaid
  genotype_cols <- suppressWarnings(!is.na(as.numeric(colnames(d))))
  id_cols <- suppressWarnings(is.na(as.numeric(colnames(d))))

  # get data with genotypes and ID
  G <- d[,genotype_cols, with = FALSE]
  id <- d[,id_cols, with = FALSE]
  id_simple <- data.table(chr = id$`#CHROM`, pos = id$POS, id = id$ID, ref = id$REF, alt = id$alt)

  # read in phenotypes
  pheno_df <- fread(args$phenotypes)
  phenotypes <- colnames(pheno_df)
  if (!is.null(args$subset_phenotypes)){
    new_phenotypes <- strsplit(args$subset_phenotypes, split = ",")
    phenotypes <- intersect(unlist(new_phenotypes), phenotypes)
  }

  cores <- detectCores()
  registerDoParallel(cores-1)
  write(paste("running parallel with",cores,"cores."), stdout())
    
  lst <- (foreach (i=1:length(phenotypes)) %dopar% {
    
    # get defined phenotypes
    pheno <- phenotypes[[i]]
    defined_phenos <- !is.na(pheno_df[[pheno]])
    eid_with_defined_phenos <- pheno_df$eid[defined_phenos]
    
    # get subset of G which contains defiend phenotypes
    G_with_defined_phenos <- colnames(G) %in% eid_with_defined_phenos
    G_subset <- G[,G_with_defined_phenos, with = FALSE]
    sums <- rowSums(G_subset)
    return(sums)
  })

  M <- data.table(do.call(cbind, lst))
  stopifnot(nrow(M) == nrow(G))
  colnames(M) <- phenotypes
  M <- cbind(id_simple, M)


  outfile = paste0(args$out_prefix, ".txt.gz")
  write(paste("Done! writing to", outfile), stdout())
  fwrite(M, outfile)
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "String. Current phenotype")
parser$add_argument("--subset_phenotypes", default=NULL, required = TRUE, help = "String. Current phenotype")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









