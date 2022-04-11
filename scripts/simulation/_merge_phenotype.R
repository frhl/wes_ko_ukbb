#!/usr/bin/env Rscript

library(argparse)
library(stringr)
library(data.table)


main <- function(args){

  print(args)
  stopifnot(dir.exists(dirname(args$input_path)))
  stopifnot(file.exists(args$real_phenotype_path))

  # load phenotypes and covariates to keep
  covars <- unlist(strsplit(args$covars_keep, split = ","))
  d <- fread(args$real_phenotype_path)
  d <- d[,colnames(d) %in% c("eid", covars), with = FALSE] 
  stopifnot(ncol(d) > 0 )

  # load prefixes
  files <- list.files(dirname(args$input_path), pattern = basename(args$input_path), full.names = TRUE)
  files <- files[!grepl('_genes',files)]
  files <- files[!grepl('_phenos',files)]

  # merge prefixes
  lst <- lapply(files, function(f){
     id <- str_extract(basename(f),"[0-9]+.tsv.gz")
     id <- gsub('\\.tsv\\.gz','', id)
     d <- fread(f)
     colnames(d)[2:ncol(d)] <- paste0(colnames(d)[2:ncol(d)],'_', id)
     return(d)
  })

  # combine phentypes and each prefix file
  reduced <- Reduce(merge, lst)
  final <- merge(reduced, d, by.x = 's', by.y = 'eid', all.x = TRUE)

  # write out
  colnames(final)[colnames(final)=="s"] <- 'eid'
  fwrite(final, args$output_path, sep = args$delimiter)


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

