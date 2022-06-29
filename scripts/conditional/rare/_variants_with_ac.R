#!/usr/bin/env Rscript

# subset to variants with at least one defined allele. This is required
# to avoid condiitoning on monomorphic markers in SAIGE

library(data.table)

main <- function(args){


  stopifnot(nchar(args$current_list) > 1) 
  stopifnot(file.exists(args$allele_count_file)) 
  d <- fread(args$allele_count_file)
  stopifnot(nrow(d) > 1)
  stopifnot(args$phenotype %in% colnames(d))
  stopifnot("id" %in% colnames(d))
  
  # read current list of variants
  current_variants <- unlist(strsplit(args$current_variants, split = ','))

  # read variants and allele count by phenotype
  bool_ac_gt_zero <- d[[args$phenotype]] > 0
  variants_with_ac_gt_zero <- d$id[bool_ac_gt_zero]
  
  # subset to variants with AC > 0 and keep string 
  current_variant <- current_variants[current_variants %in% variants_with_ac_gt_zero]
  outstring <- paste0(variants_with_ac_gt_zero, sep = ',')

  # write to stdout so that it can be stored in bash variable
  write(outstring, stdout())
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--current_list", default=NULL, required = TRUE, help = "list of currently used variants seperated by comma")
parser$add_argument("--allele_count_file", default=NULL, required = TRUE, help = "file to allele count by phenotypes")
args <- parser$parse_args()

main(args)









