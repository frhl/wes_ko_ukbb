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
  
  # calculate off-chrosome PRS contribution
  id_col <- colnames(prs)[1]
  chroms <- colnames(prs)[-1]
  stopifnot(length(chroms)==22)
  lst <- lapply(chroms, function(ch){
      dt_remove_chr <- prs[,-ch, with = FALSE]
      off_chrom_prs <- rowSums(dt_remove_chr[,-1])
      return(off_chrom_prs)
  })

  # create new data.frame
  loco_prs <- data.table(do.call(cbind, lst))
  colnames(loco_prs) <- paste0('loco_',chroms)
  loco_prs[[id_col]] <- prs[[id_col]]

  # Load phenotypes
  d <- fread(args$phenofile)
  covars <- unlist(strsplit(args$covariates, split = ","))
  
  print("covars:")
  print(covars) 
  # check that they are present in file
  stopifnot(all(covars %in% colnames(d)))
  stopifnot(args$phenotype %in% colnames(d))
  stopifnot("eid" %in% colnames(d))
  
  # combine files and write
  d <- d[,c('eid',args$phenotype, covars), with = FALSE]
  mrg <- merge(d, loco_prs, all.x = TRUE)
  
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

