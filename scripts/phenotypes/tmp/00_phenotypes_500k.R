#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# for performing residualizing
devtools::load_all('utils/modules/R/gwastools')

main <- function(args){
 
    print(args)  
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_path)))

    dt <- fread(args$input_path)
    rows_start <- nrow(dt)

    # remove invalid rows
    valid_sex <- as.logical(df$sex %in% c(1, 0)) 
    valid_age <- as.logical(!is.na(df$age))
    keep <- valid_sex & valid_age
    dt <- dt[keep, ]
    rows_final <- rows_start - nrow(dt)
    print(paste("Removed", rows_final, "invalid rows."))

    # build covariates
    dt$age2 <- dt$age ^ 2
    dt$age3 <- dt$age ^ 3
    dt$sex_age <- dt$age * ifelse(dt$sex == 1, 1, -1)
    
    # residualize cts traits
    if (!is.null(args$covariates)){
        # extract raw phenotypes
        skip_cols <- c("sequencing.batch", "age3", "WHR_adj_BMI_M", "WHR_adj_BMI_F")
        covars <- unlist(strsplit(args$covariate, split = ","))
        bool_keep <- as.logical(sapply(dt, function(x) is.numeric(x) & length(unique(x) > 100)) & 
            !(colnames(dt) %in% unique(c("eid", covars, skip_cols))) &
            !(grepl("residual",colnames(dt))) &
            !(grepl("PC", colnames(dt)))) 
        phenos <- colnames(dt)
        phenos_keep <- phenos[bool_keep]
        print(phenos_keep)
        # get rows with N/A
        na_covar_rows <- rowSums(is.na(dt[,covars, with = FALSE])) > 0
        total_rows <- nrow(dt)
        lst <- lapply(phenos_keep, function(ph){
           na_pheno_rows <- is.na(dt[[ph]])
           na_rows <- sum(na_covar_rows | na_pheno_rows)
           if (na_rows < total_rows) {
              residualize_trait_by_sex(ph, dt, "sex", covars)
           }
        })
        # append to data
        names(lst) <- paste0(phenos_keep, "_residual_sex_strat")
        stacked <- data.table(do.call(cbind, lst))
        dt <- cbind(dt, stacked)
    }

    # perform transformation
    if (!is.null(args$transform_method)){
        grep_phenos <- colnames(dt)
        grep_phenos <- grep_phenos[grepl('residual', grep_phenos)]
        f <- ifelse(args$transform_method == "int", gwastools::get_int, gwastools::get_rint)
        transformed_phenos <- lapply(grep_phenos, function(ph){ return(f(dt[[ph]])) })
        transformed_phenos <- do.call(cbind, transformed_phenos)
        colnames(transformed_phenos) <- paste0(grep_phenos, "_",args$transform_method)
        dt <- cbind(dt, transformed_phenos)
        write("Normalized residaulized phenotypes!", stdout())
    }


    write(paste0("Done! Writing to",args$out_path), stdout())
    fwrite(dt, args$out_path, sep = "\t")
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_path", default=NULL, help = "?")
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--transform_method", default=NULL, help = "?")
parser$add_argument("--covariates", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

