#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(digest)

main <- function(args){

  stopifnot(file.exists(args$in_vcf))

  # read in vcf
  cmd <- paste0("zcat ", args$in_vcf, "| grep -v '##' ")
  d <- fread(cmd = cmd)

  # setup columns that contains genotypes or id variables
  genotype_cols <- suppressWarnings(!is.na(as.numeric(colnames(d))))
  id_cols <- suppressWarnings(is.na(as.numeric(colnames(d))))

  ## if columns can't be found use this approach to identify ID 
  #total_cols <- 1:ncol(d)
  #id_cols <- 1:which(as.logical(head(d, 1) %in% "GT"))
  #id_cols <- total_cols %in% id_cols
  #genotype_cols <-  (!id_cols)

  # Check that columns are found  
  stopifnot(any(genotype_cols))
  stopifnot(any(id_cols))

  # genotype matrix
  G <- d[,genotype_cols, with = FALSE]

  # setup data.table
  id <- d[,id_cols, with = FALSE]
  if (any("DS" %in% id$FORMAT)) stop("This script only accepts genotypes!")
  dt <- data.table(chr = id$`#CHROM`, pos = id$POS, id = id$ID, ref = id$REF, alt = id$ALT)

  # combine all gentotypes into a single string per variant
  dt$gt_string <- unlist(apply(G, 1, function(x) as.character(paste(x, collapse = '-'))))
  # convert genotype to dosage
  dt$dosage_string <- gsub_to_dosage(dt$gt_string)
  # create hash of the string
  dt$dosage_hash <- unlist(lapply(dt$dosage_string, function(x) digest(x, algo="md5")))
  # clean up
  dt$gt_string <- NULL
  dt$dosage_string <- NULL
  dt$dosage_hash_dup <- as.numeric(duplicated(dt))

  # setup outfile
  outfile <- paste0(args$out_prefix, ".txt.gz")
  fwrite(dt, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix path")
args <- parser$parse_args()

main(args)









