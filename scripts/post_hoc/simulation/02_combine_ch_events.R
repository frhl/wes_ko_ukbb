#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")
source("scripts/post_hoc/utils.R")
library(argparse)
library(data.table)

main <- function(args){

    

    in_dir <- args$in_dir
    in_regex <- args$in_regex

    files <- list.files(in_dir, pattern=in_regex, full.names=TRUE)
    combined <- do.call(rbind, lapply(files, fread))

    # aggregate homozygotes and compound hets
    aggr_homs <- setDT(aggregate(homs~gene_id, FUN=mean, data=combined))
    aggr_chets <- setDT(aggregate(chets~gene_id, FUN=mean, data=combined))
    setkeyv(aggr_homs, "gene_id")
    setkeyv(aggr_chets, "gene_id")

    # combine and summarize
    sim_counts <- merge(aggr_homs, aggr_chets)
    sim_counts$homs_per_chet <- sim_counts$homs / sim_counts$chets
    sim_counts <- sim_counts[(sim_counts$homs>0) | (sim_counts$chets>0),]
    sim_counts <- sim_counts[order(sim_counts$homs_per_chet),]
    setkeyv(sim_counts, "gene_id")

    # combine with observed counts
    mrg <- merge(sim_counts, obs_counts)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(mrg, outfile, sep="\n", quote=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "?")
parser$add_argument("--in_regex", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)



