#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
source("scripts/post_hoc/utils.R") # to get genes with cis hets
library(argparse)
library(data.table)

main <- function(args){

    d <- fread(args$phenotypes_by_count)
    g <- fread(args$genes_by_count)
    p <- fread(args$phenotypes)
    stopifnot(nrow(d) > 0)
    stopifnot(nrow(g) > 0)
    stopifnot(nrow(p) > 0)

    # load knockouts
    kos <- read_ukb_wes_kos(annotation = "pLoF_damaging_missense", chrom = 1:22)
    kos$is_cis <- kos$knockout == "Compound heterozygote (cis)"
    kos$is_chet <- kos$knockout == "Compound heterozygote"
    kos$is_hom <- kos$knockout == "Homozygote"

    combos <- do.call(rbind, lapply(1:nrow(d), function(idx){
      
        # go from 
        phenotype <- d$phenotype[idx]
        gene_id_at_idx <- d$ensembl_gene_id[idx]
        hgnc_symbol <- d$hgnc_symbol[idx]
        eid_cases <- p$eid[p[[phenotype]]]

        # subset by gene
        kos_subset <- kos
        kos_subset <- kos_subset[kos_subset$gene_id == gene_id_at_idx,]
        kos_subset <- kos_subset[kos_subset$is_cis | kos_subset$is_chet | kos_subset$is_hom,]
        stopifnot(nrow(kos_subset)>0)

        # get case/control status
        d1 <- data.table(table(kos_subset[kos_subset$s %in% eid_cases,]$knockout))
        d1$cases <- TRUE
        d2 <- data.table(table(kos_subset[!kos_subset$s %in% eid_cases,]$knockout))
        d2$cases <- FALSE
        out <- rbind(d1, d2)
        out$hgnc_symbol <- hgnc_symbol
        out$gene <- gene_id_at_idx
        out$trait <- phenotype
        return(out)
        
    }))

    # we want at least two chets before testing a phenotype
    min_chet_cases <- as.numeric(args$min_chet_cases)
    combos <- combos[combos$V1 == "Compound heterozygote" & combos$cases == TRUE & combos$N >= min_chet_cases,]
    combos <- combos[,c("trait","gene","hgnc_symbol","N")]

    # write the file
    out_prefix <- paste0(args$out_prefix, ".tsv.gz")
    fwrite(combos, out_prefix, sep = '\t')
    write(paste("Note: wrote", out_prefix),stderr())

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotypes_by_count", default=NULL, required = TRUE, help = "")
parser$add_argument("--genes_by_count", default=NULL, required = TRUE, help = "")
parser$add_argument("--min_chet_cases", default=0, required = TRUE, help = "")
parser$add_argument("--phenotypes", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

