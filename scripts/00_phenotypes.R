#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){
 
    print(args)  
    stopifnot(file.exists(args$input_path))
    stopifnot(dir.exists(dirname(args$out_path)))

    dt <- fread(args$input_path)
    dt$age2 <- dt$age ^ 2
    dt$age3 <- dt$age ^ 3
    dt$sex_age <- dt$age * ifelse(dt$sex == 1, 1, -1)
    fwrite(dt, args$out_path, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_path", default=NULL, help = "?")
parser$add_argument("--input_path", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

