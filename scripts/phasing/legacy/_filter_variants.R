library(argparse)
library(stringr)

main <- function(args){

  in_file1 <- args$in_file1
  in_file2 <- args$in_file2
  V1 <- fread(in_file1, header = FALSE)$V1
  V2 <- fread(in_file2, header = FALSE)$V1
  


     


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_file1", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--in_file2", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--extract", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
parser$add_argument("--out_file", default=NULL, required = TRUE, help = "Path to .bed file used for calculating LD panel")
args <- parser$parse_args()

main(args)


