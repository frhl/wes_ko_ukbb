
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)
library(stringr)

main <- function(args){

   chrom = args$chrom
    annotation = args$annotation
    out_prefix = args$out_prefix
    path_phenotypes = args$path_phenotypes
    phenotype = args$phenotype

    # load phenotype data
    df_phenotype <- fread(path_phenotypes)
    stopifnot(phenotype %in% colnames(df_phenotype))

    # get cases and controls
    eid_cases <- df_phenotype$eid[df_phenotype[[phenotype]]]
    eid_controls <-  df_phenotype$eid[!df_phenotype[[phenotype]]]

    # read knockouts
    d <- read_ukb_wes_kos(annotation = annotation, chromosomes = chrom)
    types_to_keep <- c("Compound heterozygote (cis)", "Compound heterozygote", "Homozygote")
    d <- d[d$knockout %in% types_to_keep, ]
    d$is_cis <- d$knockout %in% "Compound heterozygote (cis)"
    d$is_chet <- d$knockout %in% "Compound heterozygote"
    d$is_hom <- d$knockout %in% "Homozygote"
    genes <- unique(d$gene_id)
    stopifnot(length(genes) > 0) 


    # iterate over phenotypes
    df_phenotype <- fread(path_phenotypes)
    
   #out <- out[out$N > 0,]
   #out <- out[!duplicated(out), ]
   #outfile <- paste0(out_prefix,".txt.gz")
   #write(paste("writing", outfile), stderr()) 
   #fwrite(out, outfile, sep = "\t")
     

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


