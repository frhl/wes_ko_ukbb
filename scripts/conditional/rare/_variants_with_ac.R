#!/usr/bin/env Rscript

# subset to variants with at least one defined allele. This is required
# to avoid condiitoning on monomorphic markers in SAIGE

library(argparse)
library(data.table)

main <- function(args){

  print(args)

  stopifnot(file.exists(args$current_marker_file)) 
  stopifnot(file.exists(args$allele_count_file)) 
  stopifnot(as.numeric(args$allele_count_threshold) >= 0) 
  d <- fread(args$allele_count_file)
  stopifnot(nrow(d) > 1)
  stopifnot(args$phenotype %in% colnames(d))
  stopifnot("id" %in% colnames(d))
  
  # read current list of variants
  #current_variants <- unlist(strsplit(args$current_variants, split = ','))
  current_variants <- fread(args$current_marker_file, header = FALSE)
  current_variants <- data.table(t(current_variants))
  rownames(current_variants) <- NULL
  colnames(current_variants) <- "id"

  # read variants and allele count by phenotype
  cutoff <- as.numeric(args$allele_count_threshold) 
  bool_ok <- floor(d[[args$phenotype]]) > cutoff
  variants_ok <- d$id[bool_ok]

  # subset to variants with AC > 0 and keep string
  #monomorphic_variants <- !(current_variants$id %in% variants_with_ac_gt_threshold)
  n_ok <- sum(bool_ok)
  n_discarded <- sum(!bool_ok)
  n_total <- nrow(current_variants)
  #write(paste("Note:",args$phenotype,"had", n_discarded,"of",n_total,"variants with AC <", cutoff, "that were discarded."),stderr())
  write(paste("Note:",args$phenotype,"had", n_ok,"of",n_total," variants with AC >", cutoff),stderr())
  fwrite(data.table(x=variants_ok), args$outfile, row.names = FALSE, col.names = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--current_marker_file", default=NULL, required = TRUE, help = "list of currently used variants seperated by comma")
parser$add_argument("--allele_count_file", default=NULL, required = TRUE, help = "file to allele count by phenotypes")
parser$add_argument("--allele_count_threshold", default=0, required = FALSE, help = "Allele count threshold, greater than or equal '>='")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









