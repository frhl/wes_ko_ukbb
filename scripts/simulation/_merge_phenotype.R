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
  bname <- basename(args$input_path)
  dname <- dirname(args$input_path)
  files <- list.files(dname, pattern = bname, full.names = TRUE)
  if (length(files) == 0) stop(paste("no files with pattern",bname,"in directory",dname))
  files <- files[!grepl('_entries',files)]
  files <- files[!grepl('_phenos',files)]

  # merge prefixes
  lst <- lapply(files, function(f){
      id <- stringr::str_extract(basename(f),"[0-9]+_cols.tsv.gz")
      id <- gsub('_cols\\.tsv\\.gz','', id) 
      d <- fread(f)
      colnames(d)[2:ncol(d)] <- paste0(colnames(d)[2:ncol(d)],'_', id)
      return(d)
  })

  # combine phentypes and each prefix file
  reduced <- Reduce(merge, lst)
  final <- merge(reduced, d, by.x = 's', by.y = 'eid', all.x = TRUE)

  # add missing covariates
  final$age2 <- final$age^2
  final$sex_age <- final$age * final$sex

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

