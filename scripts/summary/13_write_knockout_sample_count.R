# module purge
# conda activate rpy
# Rscript 13_write_knockout_sample_count.R --out_prefix derived/tables/knockouts/ukb_wes200k

# setup paths and libs
library(argparse)
library(data.table)
library(ggplot2)

devtools::load_all('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/utils/modules/R/ukbtools')

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, help = "Where should the result be written to?")
args <- parser$parse_args()


files <- list.files('derived/knockouts/211111/', full.names = TRUE)
stopifnot(dir.exists(dirname(args$out_prefix)))
mutations <- c('ptv','ptv_damaging_missense','synonymous')
mafs <- c('00_01','01_50','00_50')

# Count number of knockouts stratified by the three Maf thresholds
for (mutation in mutations){

	dt00_01 <- load_knockout_bundle(files, '00_01',mutation)
	colnames(dt00_01)[2:4] <- paste0('00_01.',colnames(dt00_01)[2:4])
	dt01_50 <- load_knockout_bundle(files, '01_50',mutation)
	colnames(dt01_50)[2:4] <- paste0('01_50.',colnames(dt01_50)[2:4])
	dt00_50 <- load_knockout_bundle(files, '00_50',mutation)
	colnames(dt00_50)[2:4] <- paste0('00_50.',colnames(dt00_50)[2:4])
	mrg <- merge(merge(dt00_01, dt01_50), dt00_50)
	outfile <- paste0(args$out_prefix, '_gene_knockout_by_maf_count_',mutation,'.tsv.gz')
	print(paste('writing',outfile))
	fwrite(mrg, outfile, quote = FALSE, sep = '\t', row.names = FALSE)
}

# write out alllele specific mutations and count them by number of individuals carrying them
for (mutation in mutations){
	for (maf in mafs){
		dt <- get_knockout_sample_counts(files, maf, mutation)
		outfile <- paste0(args$out_prefix, '_knockout_alleles_in_maf', maf, '_',mutation,'.tsv.gz')
		print(paste('writing',outfile))
		fwrite(dt, outfile, quote = FALSE, sep = '\t', row.names = FALSE)
	}
}


