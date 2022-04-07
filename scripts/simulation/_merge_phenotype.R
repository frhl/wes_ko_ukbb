#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

  print(args)
  stopifnot(file.exists(args$input_path))
  stopifnot(file.exists(args$covars_keep))
  stopifnot(file.exists(args$real_phenotype_path))

  covars <- unlist(strsplit(args$covars_keep, split = ","))
  d1 <- fread(args$input_path)
  d2 <- fread(args$real_phenotype_path)
  d2 <- d2[,colnames(d2) %in% c("eid", covars)]
  d <- merge(d1, d2, by.x = 's', by.y = 'eid', all.x = TRUE)
  colnames(d)[colnames(d)=="s"] <- 'eid'
  fwrite(d, args$output_path, sep = args$delimiter)


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--real_phenotype_path", default=NULL, help = "?")
parser$add_argument("--covars_keep", default=NULL, help = "?")
parser$add_argument("--output_path", default=NULL, help = "?")
parser$add_argument("--delimiter", default="\t", help = "?")
args <- parser$parse_args()

main(args)

