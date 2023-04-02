#!/usr/bin/env Rscript
devtools::load_all("utils/modules/R/gwastools")

library(argparse)

main <- function(args){
   phenotype <- args$phenotype
   use_bonf_corrected <- !(args$include_nominal_significant)
   all_phenotypes <- get_phenos_tested()
   stopifnot(phenotype %in% all_phenotypes)
   prs_phenotype <- get_phenos_prs(use_bonf_corrected)
   out <- ifelse(phenotype %in% prs_phenotype, 1 , 0)
   write(out, stdout())
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, help = "?")
parser$add_argument("--include_nominal_significant", default=FALSE, action="store_true", help = "?")
args <- parser$parse_args()

main(args)

