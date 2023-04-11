#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R") # to get genes with cis hets
library(argparse)
library(data.table)

main <- function(args){

    d <- fread(args$path_sig_hits)
    d$phenotype <- gsub("chr[0-9]+\\_","",d$phenotype)
    #d <- d[d$p.value < (0.05 / (311 * 1143)),]
    #d <- d[d$N_ko >= 10, ] # this is allele count and not KO (>4 KOs)
    #d <- d[d$N_ko_case >= 0]
    d$ensembl_gene_id <- d$MarkerID

    # get knockouts
    pLoF_damaging_missense <- read_ukb_wes_kos("pLoF_damaging_missense")
    cols_to_keep <- c("gene_id","s","knockout","annotation", "pKO", "chromosome", "transcript_id")
    dt <- pLoF_damaging_missense
    dt <- dt[!(dt$knockout %in% "Heterozygote"), ]
    dt <- dt[,colnames(dt) %in% cols_to_keep, with = FALSE]

    # subset to genes that are signifcant in primary analysis
    dt <- dt[dt$gene_id %in% unique(d$ensembl_gene_id),]

    # get is_chet/is_cis
    dt$is_chet <- dt$knockout %in% "Compound heterozygote"
    dt$is_cis <- dt$knockout %in% "Compound heterozygote (cis)"
    dt$is_hom <- dt$knockout %in% "Homozygote"
    dt$is_ko <- dt$is_hom | dt$is_chet
    dt$knockout <- NULL

    # cheeky rename
    colnames(dt)[colnames(dt) == "gene_id"] <- "ensembl_gene_id"
    colnames(dt)[colnames(dt) == "transcript_id"] <- "ensembl_transcript_id"

    # aggregate all by chet and cis
    aggr_cis <- aggregate(is_cis ~ ensembl_gene_id + ensembl_transcript_id, data = dt, FUN=sum)
    aggr_chet <- aggregate(is_chet ~ ensembl_gene_id + ensembl_transcript_id, data = dt, FUN=sum)
    aggr <- merge(aggr_cis, aggr_chet, by = c("ensembl_gene_id", "ensembl_transcript_id"), all = TRUE)
    aggr[is.na(aggr)] <- 0
    
    # subset to at least x/y cis chets
    aggr <- aggr[aggr$is_chet >= as.numeric(args$min_chet),]
    aggr <- aggr[aggr$is_cis >= as.numeric(args$min_cis),]
    mrg <- merge(aggr, d, by = "ensembl_gene_id")
    mrg <- mrg[,c("phenotype","ensembl_gene_id", "hgnc_symbol", "CHR", "is_cis", "is_chet")]

    # write
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(mrg, outfile, sep = "\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--path_sig_hits", default=NULL, required = TRUE, help = "")
parser$add_argument("--min_chet", default=NULL, required = TRUE, help = "")
parser$add_argument("--min_chet_case", default=0, required = FALSE, help = "")
parser$add_argument("--min_cis", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

