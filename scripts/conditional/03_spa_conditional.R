# module purge
# conda activate rpy

# setup paths and libs
library(argparse)
library(data.table)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

main <- function(args){

    d <- fread(args$input_path)
    cond <- "p.value_cond" %in% colnames(d)    
    if (!cond){
        d <- d[order(d$p.value),]
        

    } else {



    }

}



# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "what input should be read?")
parser$add_argument("--out_prefix", default=NULL, help = "Where should the results be written?")
args <- parser$parse_args()



