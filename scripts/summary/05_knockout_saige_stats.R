# module purge
# conda activate rpy

# setup paths and libs
library(argparse)
library(devtools)
library(data.table)
library(ggplot2)
library(dplyr)


devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')
setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# input
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "out prefix for files")
args <- parser$parse_args()

# static params
knockouts <- list.files('derived/knockouts/211111/', full.names = TRUE)
pheno_binary <- unlist(strsplit(readLines('data/phenotypes/UKBB_WES200k_binary_phenotypes_header.txt'), split = '\t'))
pheno_cts <- unlist(strsplit(readLines('data/phenotypes/UKBB_WES200k_cts_phenotypes_header.txt'), split = '\t'))
saige_binary <- list.files('data/saige/output/combined/binary/step2/211111/', full.names = TRUE)

# thresholds
mutations <- c('ptv','ptv_damaging_missense','synonymous')
mafs <- c('00_01','01_50','00_50')


# Combine phenotypes/mutations chromosome wise and write
for (phenotype in unique(pheno_binary)){
    print(phenotype)
    outlist <- list()
    for (maf in mafs){
        print(maf)
        for (mutation in mutations){
            id <- paste0(phenotype,maf,mutation)
            dts <- load_saige_bundle(saige_binary, maf, mutation, phenotype)
            kos <- load_knockout_bundle(knockouts, maf, mutation)
            mrg <- merge(dts, kos, by = 'gene_id', all.x = TRUE)
            if (nrow(kos) > 0 & nrow(dts) > 0){
		mrg$phenotype <- phenotype
		mrg$maf <- maf 
		mrg$mutation <- mutation
		outlist[[id]] <- mrg
	    } else {
		print(paste0('skipping',phenotype))
        }
    }   
    combined <- as.data.table(do.call(rbind, outlist))
    outfile = paste0(args$out_prefix,'_', phenotype,'.txt.gz')
    print(outfile)
    fwrite(combined, outfile, sep = '\t', quote = FALSE, row.names = FALSE)
}







