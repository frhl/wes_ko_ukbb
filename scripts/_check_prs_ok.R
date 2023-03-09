#!/usr/bin/env Rscript
devtools::load_all("utils/modules/R/gwastools")

library(argparse)

main <- function(args){
   phenotype <- args$phenotype
   all_phenotypes <- gwastools::get_phenos_tested()
   stopifnot(phenotype %in% all_phenotypes)
   prs_phenotype <- fread(gwastools::get_phenos_prs_path(), header=FALSE)$V1
   out <- ifelse(phenotype %in% prs_phenotype, 1 , 0)
   write(out, stdout())
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

