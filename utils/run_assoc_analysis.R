#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

## load libs
library(data.table)
devtools::load_all('/well/lindgren/flassen/code/our/')

## load data
pheno_df = fread(args[1]) #fread(IN_PHENO)
ko_df = fread(args[2], sep = ',') #fread(IN_KO) # REMEBER TO SUBSET BY TRANS HERE!!
goi = args[3]
in_sge_task_id = as.numeric(args[4])
out_dir = args[5]
in_covars = args[6]
in_analysis = args[7]
in_encoding = args[8]

# setup params
analysis = in_analysis
encoding = in_encoding

# check files
stopifnot(dir.exists(out_dir))

# check input
if (ncol(ko_df) > 6) stop('too many columns were inputted!')
if (ncol(ko_df) < 6) stop('too few columns were inputted!')

# extract data
colnames(ko_df) <- c('ID','gene','knockout','strand','chrom', 'genotype')

# combine phenotype data and knockout data
merged = ukbb_format_assoc(pheno_df, ko_df, goi, analysis = analysis, encode_knockout = encoding)

# extract data
#colnames(ko_df) <- c('ID','gene','knockout','strand','chrom', 'genotype')
#merged = ukbb_format_assoc(pheno_df, ko_df, goi, analysis = analysis)

# since we are merging knockouted individuals with non-knockouted indivuals,
# and avoid iterating over other individuals, assume that these are non-knockouts
merged$data[[goi]][is.na(merged$data[[goi]])] <- 0

# since we are using tereas phenotypes, which are slightly different QCed. merge
# with out QCWB file to ensure that the right individuals are included.
#merged$data <- merge(merged$data$ID, df_qcwb, by = 'ID')
merged$data <- merged$data[merged$data$white.british == 1,]

# setup trait
trait = merged$traits[as.numeric(in_sge_task_id)]
stopifnot(!is.na(trait))

# fit model
#covars <- c('age','sex', 'genotyping.array', 'ukbb.centre', paste0('PC',1:4))
covars <- unlist(strsplit(readLines(in_covars), split = ' '))
fit <- ukbb_run_assoc(merged$data, trait, goi, covars, model = in_analysis)

# write out
outfile_full = file.path(out_dir, paste0('full.',goi,'_',trait,'.',analysis,'.assoc'))
outfile_minimal = file.path(out_dir, paste0('minimal.',goi,'_',trait,'.',analysis,'.assoc'))
write.table(fit$full, file = outfile_full, quote = F)
write.table(fit$minimal, file = outfile_minimal, quote = F)



