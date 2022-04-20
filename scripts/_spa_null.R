#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){
    
  print(args)
  stopifnot(file.exists(args$prsfile))
  stopifnot(file.exists(args$phenofile))
  stopifnot(args$covariates != "")
  stopifnot(args$phenotype != "")

  
  # load PRS and phenotypes
  prs <- fread(args$prsfile)
  colnames(prs)[1] <- 'eid'
  d <- fread(args$phenofile)
  covars <- unlist(strsplit(args$covariates, split = ","))
  
  # check that they are present in file
  stopifnot(all(covars %in% colnames(d)))
  stopifnot(args$phenotype %in% colnames(d))
  stopifnot("eid" %in% colnames(d))
  
  # combiune files and write
  d <- d[,c('eid',args$phenotype, covars), with = FALSE]
  mrg <- merge(d, prs, all.x = TRUE)
  fwrite(mrg, args$outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, help = "?")
parser$add_argument("--covariates", default=NULL, help = "?")
parser$add_argument("--prsfile", default=NULL, help = "?")
parser$add_argument("--phenofile", default=NULL, help = "?")
parser$add_argument("--outfile", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

