#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# for performing residualizing
devtools::load_all('utils/modules/R/gwastools')

main <- function(args){
 
    print(args)  
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_path)))
    stopifnot(!is.null(args$covariates))
 
    # append new covariates
    dt <- fread(args$input_path)
    dt$age2 <- dt$age ^ 2
    dt$age3 <- dt$age ^ 3
    dt$sex_age <- dt$age * ifelse(dt$sex == 1, 1, -1)

    # only keep samples in which we have no missing covariates 
    covars <- unlist(strsplit(args$covariates, split = ","))
    missing <- rowSums(do.call(cbind, lapply(covars, function(col) is.na(dt[[col]]) ))) > 0
    n_missing <- sum(missing)
    if ((n_missing > 0) & (args$row_na_action == "remove")) {
        write(paste("Note: dropping", n_missing, "sample(s) with missing covariates."),stderr())
        dt <- dt[!missing, ]
    }

    # transform phenotypes (either RINT or INT)
    if (!is.null(args$transform_method) & !is.null(args$transform)){
       
        # only transform selected phenotypes 
        phenotypes <- unlist(strsplit(args$transform, split = ","))
        phenotypes <- phenotypes[phenotypes %in% colnames(dt)]
        stopifnot(length(phenotypes) > 0)

        # do transformation and aggregate
        f <- ifelse(args$transform_method == "int", gwastools::get_int, gwastools::get_rint)
        transformed_phenos <- lapply(phenotypes, function(ph){ return(f(dt[[ph]])) })
        transformed_phenos <- do.call(cbind, transformed_phenos)
        colnames(transformed_phenos) <- paste0(phenotypes, "_",args$transform_method)
        dt <- cbind(dt, transformed_phenos)
    }

    # Finally write file
    write(paste0("Done! Writing to",args$out_path), stdout())
    fwrite(dt, args$out_path, sep = "\t")
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_path", default=NULL, help = "Where should the file be written")
parser$add_argument("--input_path", default=NULL, help = "Curated phenotypes input path")
parser$add_argument("--transform_method", default=NULL, help = "Transformation of residuals, either INT or RINT")
parser$add_argument("--transform", default=NULL, help = "What phenotypes should be transformed?")
parser$add_argument("--row_na_action", default="keep", help = "set to 'remove' to delete rows with missing covariates")
parser$add_argument("--covariates", default=NULL, help = "comma seperated string of covariates")
args <- parser$parse_args()

main(args)

