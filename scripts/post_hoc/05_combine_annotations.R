
library(data.table)
library(argparse)
source("scripts/post_hoc/utils.R")

main <- function(args){
   
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
    synonymous$annotation <- "other_missense"
    
    # keep only relevant columns
    cols_to_keep <- c("gene_id","s","knockout","annotation", "pKO", "chromosome", "transcript_id")
    dt <- setDT(rbind(pLoF, damaging_missense, pLoF_damaging_missense))
    dt <- dt[!(dt$knockout %in% "Heterozygote"), ]
    dt <- dt[,colnames(dt) %in% cols_to_keep, with = FALSE]

    # format
    dt$is_chet <- dt$knockout %in% "Compound heterozygote"
    dt$is_cis <- dt$knockout %in% "Compound heterozygote (cis)"
    dt$is_hom <- dt$knockout %in% "Homozygote"
    dt$is_het <- dt$knockout %in% "Heterozygote"
    dt$is_ko <- dt$is_hom | dt$is_chet
    dt$knockout <- NULL

    # no exclusion
    outfile <- paste0(args$out_prefix,"all.txt.gz")
    write(paste("writing", outfile), stderr())
    fwrite(dt_out, outfile)

    # exclude hets
    annotations_minimal <- c("pLoF", "pLoF_damaging_missense", "damaging_missense")
    dt <- dt[!dt$is_het, ]
    dt <- dt[dt$annotation %in% annotations_minimal, ]

    # exclude common knockoiuts
    counts <- data.table(table(dt$is_ko, dt$gene_id))
    counts <- counts[counts$V1 == TRUE,]
    counts <- counts[rev(order(counts$N))]
    gene_id_exclude <- counts$V2[counts$N > 10000]
    dt <- dt[!(dt$gene_id %in% gene_id_exclude),]
    message <- paste("excluded", gene_id_exclude)
    print(message)

    # writ efile
    outfile <- paste0(args$out_prefix,"nohets.nocommonkos.txt.gz")
    write(paste("writing", outfile), stderr())
    fwrite(dt_out, outfile)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


