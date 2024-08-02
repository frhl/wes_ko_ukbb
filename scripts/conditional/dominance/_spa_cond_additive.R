#!/usr/bin/env Rscript

devtools::load_all("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/gwastools")
library(argparse)
library(data.table)

main <- function(args){

    phenotype <- args$phenotype
    gene <- args$gene_id
    outfile <- args$outfile
    d_gene <- fread(args$path_markers)
    stopifnot(ncol(d_gene)>1)
    stopifnot(nrow(d_gene)>0)

    # gene needs to be the correct encoding
    print(head(d_gene))
    print(gene)


    # subset to pseudo marker
    d_gene <- d_gene[grepl(d_gene$ensembl_gene_id, pattern=tolower(gene)),]
    stopifnot(nrow(d_gene) > 0)
    markers <- d_gene$pseudo_marker
    fwrite(data.table(x=markers), outfile, row.names=FALSE, col.names=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--gene_id", default=NULL, required = TRUE, help = "phenotype, single string.")
parser$add_argument("--path_markers", default=NULL, required = FALSE, help = "file containing the markers in the chromosome for the current VCF")
parser$add_argument("--outfile", default=NULL, required = TRUE, help = "where should the subsetted markeres be written")
args <- parser$parse_args()

main(args)









