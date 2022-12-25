
library(data.table)
library(argparse)
library(pscl)

main <- function(args){
   
    fwrite( coeffecients, outfile, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


