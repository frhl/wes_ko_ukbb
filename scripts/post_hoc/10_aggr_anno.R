
library(data.table)
library(argparse)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")

main <- function(args){

    # retrive them
    pLoF <- read_ukb_wes_kos("pLoF")
    damaging_missense <- read_ukb_wes_kos("damaging_missense")
    pLoF_damaging_missense <- read_ukb_wes_kos("pLoF_damaging_missense")
    other_missense <- read_ukb_wes_kos("other_missense")
    synonymous <- read_ukb_wes_kos("synonymous")

    # annotate
    pLoF$annotation <- "pLoF"
    damaging_missense$annotation <- "damaging_missense"
    pLoF_damaging_missense$annotation <- "pLoF_damaging_missense"
    other_missense$annotation <- "other_missense"
    synonymous$annotation <- "synonymous"

    # keep only relevant columns
    cols_to_keep <- c("gene_id","s","knockout","annotation", "pKO", "chromosome", "transcript_id")
    dt <- setDT(rbind(pLoF, damaging_missense, pLoF_damaging_missense))
    #dt <- dt[!(dt$knockout %in% "Heterozygote"), ]
    dt <- dt[,colnames(dt) %in% cols_to_keep, with = FALSE]

    dt$is_chet <- dt$knockout %in% "Compound heterozygote"
    dt$is_cis <- dt$knockout %in% "Compound heterozygote (cis)"
    dt$is_hom <- dt$knockout %in% "Homozygote"
    dt$is_het <- dt$knockout %in% "Heterozygote"
    dt$is_ko <- dt$is_hom | dt$is_chet
    dt$knockout <- NULL

    counts <- data.table(table(dt$is_ko, dt$gene_id))
    counts <- counts[counts$V1 == TRUE,]
    counts <- counts[rev(order(counts$N))]
    counts$excluded <- counts$N > 20000
    gene_id_exclude <- counts$V2[counts$N > 20000]
    dt <- dt[!(dt$gene_id %in% gene_id_exclude),]
    outfile <- paste0(args$out_prefix, ".excl.primary.count.txt.gz")
    fwrite(counts, outfile, sep = "\t")

    # get other_missense and synonymous
    dt2 <- setDT(rbind(other_missense, synonymous))
    #dt2 <- dt2[!(dt2$knockout %in% "Heterozygote"), ]
    dt2 <- dt2[,colnames(dt2) %in% cols_to_keep, with = FALSE]

    dt2$is_chet <- dt2$knockout %in% "Compound heterozygote"
    dt2$is_cis <- dt2$knockout %in% "Compound heterozygote (cis)"
    dt2$is_hom <- dt2$knockout %in% "Homozygote"
    dt2$is_het <- dt2$knockout %in% "Heterozygote"
    dt2$is_ko <- dt2$is_hom | dt2$is_chet
    dt2$knockout <- NULL

    counts2 <- data.table(table(dt2$is_ko, dt2$gene_id))
    counts2 <- counts2[counts2$V1 == TRUE,]
    counts2 <- counts2[rev(order(counts2$N))]
    counts2$excluded <- counts2$N > 100000
    gene_id_exclude2 <- counts2$V2[counts2$N > 100000]
    #dt2 <- dt2[!(dt2$gene_id %in% gene_id_exclude2),]
    #outfile <- paste0(args$out_prefix, ".excl.other.count.txt.gz")
    #fwrite(counts2, outfile, sep = "\t")

    # write file containing all annotations 
    dt <- rbind(dt, dt2)
    colnames(dt)[colnames(dt) == "gene_id"] <- "ensembl_gene_id"
    colnames(dt)[colnames(dt) == "transcript_id"] <- "ensembl_transcript_id"
    #outfile <- paste0(args$out_prefix, ".txt.gz")
    #fwrite(dt, outfile, sep = "\t")

    # aggregate counts
    cols_aggr <- c("ensembl_gene_id", "ensembl_transcript_id", "annotation")
    aggr_cis <- aggregate(is_cis ~ ensembl_gene_id + ensembl_transcript_id + annotation, data = dt, FUN=sum)
    aggr_chet <- aggregate(is_chet ~ ensembl_gene_id + ensembl_transcript_id + annotation, data = dt, FUN=sum)
    aggr_hom <- aggregate(is_hom ~ ensembl_gene_id + ensembl_transcript_id + annotation, data = dt, FUN=sum)
    aggr_het <- aggregate(is_het ~ ensembl_gene_id + ensembl_transcript_id + annotation, data = dt, FUN=sum)
    aggr_mrg <- merge(aggr_chet, aggr_hom, by = cols_aggr, all = TRUE)
    aggr_mrg <- merge(aggr_mrg, aggr_cis, by = cols_aggr, all = TRUE)
    aggr_mrg <- merge(aggr_mrg, aggr_het, by = cols_aggr, all = TRUE)

    outfile <- paste0(args$out_prefix, ".counts.txt.gz")
    fwrite(aggr_mrg, outfile, sep = "\t")



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


