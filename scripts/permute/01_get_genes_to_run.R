#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R") # to get genes with cis hets
library(argparse)
library(data.table)

bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_hgnc <- bridge$hgnc_symbol
names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id


main <- function(args){

    kos <- read_ukb_wes_kos(annotation = "pLoF_damaging_missense", chrom = 1:22)
    kos$is_cis <- kos$knockout == "Compound heterozygote (cis)"
    kos$is_chet <- kos$knockout == "Compound heterozygote"
    genes_to_run <- Reduce(merge, 
       list(
           aggregate(is_cis ~ gene_id, data = kos, FUN=sum), 
           aggregate(is_chet ~ gene_id, data = kos, FUN=sum)
           )
       )
    colnames(genes_to_run)[1] <- "MarkerID"

    # subset by chets/cis counts
    min_cis <- as.numeric(args$min_cis)
    min_chet <- as.numeric(args$min_chet) 
    markers_keep <- genes_to_run$MarkerID[(genes_to_run$is_cis >= min_cis) & (genes_to_run$is_chet >= min_chet)]
    genes_to_run <- genes_to_run[genes_to_run$MarkerID %in% markers_keep, ]
    genes_to_run$CHR <- paste0("chr",ensembl_to_contig[genes_to_run$MarkerID])
    genes_to_run$hgnc_symbol <- ensembl_to_hgnc[genes_to_run$MarkerID]

    print(head(genes_to_run))
    # order as previous scripts
    my_cols <- c("MarkerID", "CHR", "hgnc_symbol","is_chet", "is_cis")
    genes_to_run <- genes_to_run[,colnames(genes_to_run) %in% my_cols] 
    genes_to_run <- genes_to_run[,c(1,4,5,2,3)]

    # write the file
    out_prefix_genes <- paste0(args$out_prefix, ".tsv.gz")
    fwrite(genes_to_run, out_prefix_genes, sep = '\t')
    write(paste("Note: wrote", out_prefix_genes),stderr())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--min_chet", default=NULL, required = TRUE, help = "")
parser$add_argument("--min_cis", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

