# module purge
# conda activate rpy
# Rscript 04_knockouts.R --in_dir derived/knockouts/all/211013_ptv --out_prefix derived/summary

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
#parser <- ArgumentParser()
#parser$add_argument("--in_dir", default=NULL, help = "What directory should files be loaded from?")
#args <- parser$parse_args()

# read in .tsv.bgz files
zcat <- function(path,...){
    cmd = paste("zcat", path)
    return(fread(cmd = cmd,...))
}

# how many of the variants are in gnomAD and GB
summarize_external_db <- function(d, as_pct = TRUE){
    
    stopifnot('inGnomAD' %in% colnames(d))
    stopifnot('gb_AC' %in% colnames(d))
    stopifnot('locus' %in% colnames(d))
    chrom <- unique(unlist(strsplit(d$locus, split = ':'))[1])
    total <- nrow(d)
    in_gnomad <- sum(d$inGnomAD)
    in_gb <- sum(!is.na(d$gb_AC))
    in_both <- sum(d$inGnomAD & !is.na(d$gb_AC))
    if (as_pct){
        in_gnomad <- in_gnomad / total
        in_gb <- in_gb / total
        in_both <- in_both / total
    }
    return(data.frame(chrom, in_gnomad, in_gb, in_both, total))
}


# read in external files
files_external <- sort(list.files('data/variants/',pattern = 'ukb_wes_200k_external_qc_chr[0-9]+.tsv.bgz', full.names = TRUE))
data_external <- do.call(rbind, lapply(files_external, function(file){
   d <- zcat(file)
   d <- summarize_external_db(d)
   return(d)
}))
print(data_external)
write.table(data_external, 'derived/external_qc_comparison.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# Sample-level stats
files_sample <- list.files('data/qc/',pattern = 'ukb_wes_200k_chr[0-9]+_variants_summary_phased.tsv.bgz', full.names = TRUE)
d_sample <- do.call(rbind,lapply(files_sample, function(file) zcat(file)))
d_sample_aggr <- aggregate(.~s, d_sample, sum)
write.table(d_sample_aggr, 'derived/tables/ukb_phased_sample_level_urv_count.tsv', quote = FALSE, row.names = FALSE) 

# Gene-level stats
files_genes <- list.files('data/qc/',pattern = 'ukb_wes_200k_chr[0-9]+_urv_by_genes_phased.tsv.bgz', full.names = TRUE)
d_genes <- do.call(rbind,lapply(files_genes, function(file) zcat(file)))
d_genes$s <- NULL
d_genes_aggr <- aggregate(.~gene_id, d_genes, sum)
write.table(d_genes_aggr, 'derived/tables/ukb_phased_gene_level_urv_count.tsv', quote = FALSE, row.names = FALSE) 





