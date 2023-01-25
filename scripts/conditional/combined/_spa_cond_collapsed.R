#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

    phenotype <- args$phenotype
    min_mac <- as.numeric(args$min_mac)
    gene <- args$gene_id
    outfile <- args$outfile
    d_gene <- fread(args$path_markers)

    # subset to pseudo marker
    d_gene <- d_gene[d_gene$ensembl_gene_id %in% gene,]
    stopifnot(nrow(d_gene) > 0)
    markers <- d_gene$psuedo_marker
    fwrite(data.table(x=markers), outfile, row.names=FALSE, col.names=FALSE)


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--gene_id", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--path_markers", default=NULL, required = FALSE, help = "file containing the markers in the chromosome for the current VCF")
parser$add_argument("--min_mac", default=1, required = FALSE, help = "Allele count threshold, greater than or equal '>='")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









