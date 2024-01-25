#!/usr/bin/env Rscript

devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")
source("scripts/post_hoc/utils.R")
library(argparse)
library(data.table)
library(ggplot2)
library(ggpubr)

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

    ## get counts for actual observed CHs and HOMs
    d_plof <- read_ukb_wes_kos("pLoF")

    # count Homs
    hom_counts <- data.table(table(d_plof$gene_id[d_plof$knockout %in% "Homozygote"]))
    colnames(hom_counts) <- c("gene_id", "obs_homs")
    setkeyv(hom_counts, "gene_id")

    # count CHs
    ch_counts <- data.table(table(d_plof$gene_id[d_plof$knockout %in% "Compound heterozygote"]))
    colnames(ch_counts) <- c("gene_id", "obs_ch")
    setkeyv(ch_counts, "gene_id")

    # merge CH and homs
    obs_counts <- merge(hom_counts, ch_counts, all=TRUE)
    obs_counts[is.na(obs_counts)] <- 0
    obs_counts <- obs_counts[rev(order(obs_counts$obs_homs, obs_counts$obs_ch)),]
    obs_counts$obs_homs_per_ch <- obs_counts$obs_homs / obs_counts$obs_ch 

    # combine with observed counts
    mrg <- merge(sim_counts, obs_counts)
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(mrg, outfile, sep="\t", quote=FALSE)

    # plot the ratios
    p1 <- ggplot(mrg, aes(x=homs_per_chet, y=obs_homs_per_ch)) +
        geom_point() +
        geom_abline() +
        stat_cor(method="pearson") + 
        ggtitle("Ratio of homozygotes to CHs") +
        ylab("observed homozygotes per CHs") +
        xlab("simulated homozygotes per CHs") +
        theme_bw()

    mrg <- merge(sim_counts, obs_counts)
    p2 <- ggplot(mrg, aes(x=homs, y=obs_homs)) +
        geom_point() +
        geom_abline() +
        stat_cor(method="pearson", digits=4) +
        ggtitle("Homs Only") +
        ylab("n observed homozygotes") +
        xlab("n simulated homozygotes") +
        theme_bw()

    mrg <- merge(sim_counts, obs_counts)
    p3 <- ggplot(mrg, aes(x=chets, y=obs_ch)) +
        geom_point() +
        geom_abline() +
        stat_cor(method="pearson", digits=4) +
        ggtitle("CHs only") +
        ylab("n observed CHs") +
        xlab("n simulated CHs") +
        theme_bw()

    # plot the results
    ggsave(paste0(args$out_prefix,".ratio.pdf"), p1, width=5, height=5)
    ggsave(paste0(args$out_prefix,".homs.pdf"), p2, width=5, height=5)
    ggsave(paste0(args$out_prefix,".chs.pdf"), p3, width=5, height=5)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_dir", default=NULL, help = "?")
parser$add_argument("--in_regex", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)



