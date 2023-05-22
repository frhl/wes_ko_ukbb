#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")
source("scripts/post_hoc/utils.R")

library(argparse)
library(data.table)

main <- function(args){
    
    input_path <- args$input_path
    outfile <- args$outfile
    seed <- as.numeric(args$seed)
   
    d <- fread(input_path)

    # get current phenos to shuffle
    phenos_to_shuffle <- get_phenos_tested()
 
    set.seed(seed)
    for (pheno in phenos_to_shuffle){
        print(pheno)
        print(head(which(d[[pheno]])))
        shuffled_phenos <- sample(d[[pheno]])
        d[[pheno]] <- shuffled_phenos
        print(head(which(d[[pheno]])))
    }

    write(paste("writing to", outfile), stdout())   
    fwrite(d, outfile, quote = FALSE, sep = "\t") 
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--outfile", default=NULL, help = "")
parser$add_argument("--seed", default=NULL, help = "")
args <- parser$parse_args()

main(args)

