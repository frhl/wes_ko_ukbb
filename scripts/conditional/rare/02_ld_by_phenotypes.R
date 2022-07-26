#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(parallel)
library(doParallel)
library(digest)

geno_dosage_hash <- function(G, pheno_df, pheno){
    #' @param G matrix of genotypes (numerics)
    #' @param pheno_df data.table of phenotypes
    #' @param pheno current phenotype (string)
    
    # get defined phenotypes
    defined_phenos <- !is.na(pheno_df[[pheno]])
    eid_with_defined_phenos <- pheno_df$eid[defined_phenos]
    
    # get subset of G which contains defiend phenotypes
    G_with_defined_phenos <- colnames(G) %in% eid_with_defined_phenos
    G_subset <- G[,G_with_defined_phenos, with = FALSE]

    # combine all dosages into a single variant
    dt <- data.table(
        dosage_string = unlist(apply(G_subset, 1, function(x) as.character(paste(x, collapse = '-'))))
    )
    # create hash of the string
    dt$dosage_hash <- unlist(lapply(dt$dosage_string, function(x) digest(x, algo="md5")))
    # clean up
    dt$gt_string <- NULL
    dt$dosage_string <- NULL
    dt$dosage_hash_dup <- as.numeric(duplicated(dt))
    return(dt$dosage_hash)

}



main <- function(args){


  stopifnot(file.exists(args$in_vcf)) 
  stopifnot(file.exists(args$phenotypes))
  stopifnot(file.exists(args$covariates))
  stopifnot(dir.exists(dirname(args$out_prefix)))

  # read in VCF
  cmd <- paste0("zcat ", args$in_vcf, " | grep -v '##'")
  d <- fread(cmd = cmd)

  # get columns with genotypes / metaid
  genotype_cols <- suppressWarnings(!is.na(as.numeric(colnames(d))))
  id_cols <- suppressWarnings(is.na(as.numeric(colnames(d))))

  # Need to ensure that G is all numerics
  G <- d[,genotype_cols, with = FALSE]
  G[, names(G) := lapply(.SD, as.numeric)]
  
  # get ID columns
  id <- d[,id_cols, with = FALSE]
  id_simple <- data.table(chr = id$`#CHROM`, pos = id$POS, id = id$ID, ref = id$REF, alt = id$ALT)

  # read in phenotypes
  pheno_df <- fread(args$phenotypes)
  phenotypes <- colnames(pheno_df)
  if (!is.null(args$subset_phenotypes)){
    new_phenotypes <- strsplit(args$subset_phenotypes, split = ",")
    phenotypes <- intersect(unlist(new_phenotypes), phenotypes)
  }

  # are there samples with missing covariates?
  col_cov <- unlist(strsplit(readLines(args$covariates), split = ","))
  lst <- lapply(col_cov, function(col){row_ok <- is.na(pheno_df[[col]])})
  missing_cov <- rowSums(do.call(cbind, lst)) > 0
  pheno_df <- pheno_df[!missing_cov, ]
  msg <- paste("Removed", sum(missing_cov),"samples with missing covariates.")
  write(msg, stdout())

  # get hash of dosages to count duplicates (find those genotypes
  # that are in perfect LD)
  lst_ds_hash <- lapply(1:length(phenotypes), function(i){
     geno_dosage_hash(G, pheno_df, phenotypes[i])
  }) 

    # combine files
  dt_ds_hash <- data.table(do.call(cbind, lst_ds_hash))
  
  # ensure that they are same lenght
  stopifnot(nrow(dt_ds_hash) == nrow(G))
  
  # clean up header
  colnames(dt_ds_hash) <- phenotypes
  dt_ds_hash <- cbind(id_simple, dt_ds_hash)

  
  outfile = paste0(args$out_prefix, ".txt.gz")
  write(paste("Done! writing to", outfile), stdout())
  fwrite(dt_ds_hash, outfile)
  
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, required = TRUE, help = "either 'binary' or 'cts'")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "String. Current phenotype")
parser$add_argument("--covariates", default=NULL, required = TRUE, help = "Covariates seperated by comma")
parser$add_argument("--subset_phenotypes", default=NULL, required = TRUE, help = "String. Current phenotype")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)









