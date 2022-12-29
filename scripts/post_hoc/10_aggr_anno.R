
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
    dt <- dt[!(dt$knockout %in% "Heterozygote"), ]
    dt <- dt[,colnames(dt) %in% cols_to_keep, with = FALSE]

    dt$is_chet <- dt$knockout %in% "Compound heterozygote"
    dt$is_cis <- dt$knockout %in% "Compound heterozygote (cis)"
    dt$is_hom <- dt$knockout %in% "Homozygote"
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
    dt2 <- dt2[!(dt2$knockout %in% "Heterozygote"), ]
    dt2 <- dt2[,colnames(dt2) %in% cols_to_keep, with = FALSE]

    dt2$is_chet <- dt2$knockout %in% "Compound heterozygote"
    dt2$is_cis <- dt2$knockout %in% "Compound heterozygote (cis)"
    dt2$is_hom <- dt2$knockout %in% "Homozygote"
    dt2$is_ko <- dt2$is_hom | dt2$is_chet
    dt2$knockout <- NULL

    counts2 <- data.table(table(dt2$is_ko, dt2$gene_id))
    counts2 <- counts2[counts2$V1 == TRUE,]
    counts2 <- counts2[rev(order(counts2$N))]
    counts2$excluded <- counts2$N > 50000
    gene_id_exclude2 <- counts2$V2[counts2$N > 50000]
    dt2 <- dt2[!(dt2$gene_id %in% gene_id_exclude2),]
    outfile <- paste0(args$out_prefix, ".excl.other.count.txt.gz")
    fwrite(counts2, outfile, sep = "\t")

    # write file containing all annotations 
    dt <- rbind(dt, dt2)
    print(table(dt$annotation))
    outfile <- paste0(args$out_prefix, ".nohets.txt.gz")
    fwrite(dt, outfile, sep = "\t")

    # aggregate all by chet and cis
    full <- dt[dt$is_chet | dt$is_cis, ]
    aggr_ko <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=sum)
    aggr_full <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=length)
    aggr <- merge(aggr_ko, aggr_full, by = c("gene_id","transcript_id", "annotation"))
    colnames(aggr) <- c("ensembl_gene_id", "ensembl_transcript_id", "annotation", "n","total")
    aggr$subset <- "chet+cis"
    aggr_chet_cis <- aggr

    # aggregate all by hom and cis
    full <- dt[dt$is_hom | dt$is_cis, ]
    aggr_ko <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=sum)
    aggr_full <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=length)
    aggr <- merge(aggr_ko, aggr_full, by = c("gene_id","transcript_id", "annotation"))
    colnames(aggr) <- c("ensembl_gene_id", "ensembl_transcript_id", "annotation", "n","total")
    aggr$subset <- "hom+cis"
    aggr_hom_cis <- aggr

    # aggregate all by knockout and cis
    full <- dt[dt$is_ko | dt$is_cis, ]
    aggr_ko <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=sum)
    aggr_full <- aggregate(pKO ~ gene_id + transcript_id + annotation, data = full, FUN=length)
    aggr <- merge(aggr_ko, aggr_full, by = c("gene_id","transcript_id", "annotation"))
    colnames(aggr) <- c("ensembl_gene_id", "ensembl_transcript_id", "annotation", "n","total")
    aggr$subset <- "ko+cis"
    aggr_ko_cis <- aggr

    # get transcript lengths and mapping
    transcript <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/221216_enstid_ensgid_lengths.txt.gz")
    transcript <- transcript[,c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "length")]
    transcript$norm_length <- (transcript$length-mean(transcript$length))/sd(transcript$length)
    transcript$norm_length2 <- transcript$norm_length ^ 2

    # merge with transcript
    combined <- rbind(aggr_chet_cis, aggr_hom_cis, aggr_ko_cis)
    out <- merge(combined, transcript, all.x = TRUE)
    stopifnot(nrow(combined) == nrow(out))

    # write file containing all annotations 
    outfile <- paste0(args$out_prefix, ".counts.txt.gz")
    print(table(out$annotation))
    fwrite(out, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


