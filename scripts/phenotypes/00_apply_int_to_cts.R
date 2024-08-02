#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# for performing residualizing
devtools::load_all('utils/modules/R/gwastools')

main <- function(args){
 
    print(args)  
    stopifnot(file.exists(args$in_file))
    stopifnot(dir.exists(dirname(args$out_file)))
  
    dt <- fread(args$in_file)

    # only transform selected phenotypes 
    phenotypes <- unlist(strsplit(args$phenotypes_list, split = ","))
    phenotypes <- phenotypes[phenotypes %in% colnames(dt)]
    stopifnot(length(phenotypes) > 0)

    # Transformation function selection
    f <- ifelse(args$transform_method == "int", gwastools:::get_int, gwastools:::get_rint)

    # In-place transformation
    for (ph in phenotypes) {
        dt[, (ph) := lapply(.SD, f), .SDcols = ph]
    }

    # Finally write file
    write(paste0("Done! Writing to ",args$out_file), stdout())
    fwrite(dt, args$out_file, sep = "\t", quote = FALSE)
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_file", default=NULL, help = "Where should the file be written")
parser$add_argument("--in_file", default=NULL, help = "Curated phenotypes input path")
parser$add_argument("--transform_method", default="int", help = "Transformation of residuals, either 'int' or 'rint'. SAIGE uses 'int'.")
parser$add_argument("--phenotypes_list", default=NULL, help = "list of phenotypes to apply transformation to")
args <- parser$parse_args()

main(args)

