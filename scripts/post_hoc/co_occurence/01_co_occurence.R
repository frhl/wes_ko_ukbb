
devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)


main <- function(args){

    # get things to keep
    chrom = args$chrom
    annotation = args$annotation
    out_prefix = args$out_prefix

    # read knockouts
    d <- read_ukb_wes_kos(annotation = annotation, chromosomes = chrom)
    types_to_keep <- c("Compound heterozygote (cis)", "Compound heterozygote", "Homozygote")
    d$is_cis <- d$knockout %in% "Compound heterozygote (cis)"
    d$is_chet <- d$knockout %in% "Compound heterozygote"
    d$is_hom <- d$knockout %in% "Homozygote"
    d <- d[(d$is_cis) | (d$is_chet) | (d$is_hom), ]
    
    # get co-occurence in almost linear time
    lib <- gwastools::cooccur_pack_lib(d)
    out <- gwastools::cooccur_unpack_lib(lib)
   
    outfile <- paste0(out_prefix,".txt.gz")
    fwrite(out, outfile, sep = "\t")
 
}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


