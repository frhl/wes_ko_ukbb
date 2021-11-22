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
parser$add_argument("--out_prefix", default='derived/tables/worst_csq_for_variant_canonical/summary', help = "Dest dir")
args <- parser$parse_args()

# Get count of worst consequence for each variant
files <- list.files('derived/tables/worst_csq_for_variant_canonical/', pattern = 'ukb_wes_200k_unphased_chr[0-9]+_category\\.tsv', full.names = TRUE)
files_singletons <- list.files('derived/tables/worst_csq_for_variant_canonical/', pattern = 'ukb_wes_200k_unphased_chr[0-9]+_category_singletons\\.tsv', full.names = TRUE)
print(files)

# Aggregate each category across chromosomes
d1 <- aggregate(n ~ consequence_category, data = do.call(rbind,lapply(files, fread)), FUN = sum)
d2 <- aggregate(n ~ consequence_category, data = do.call(rbind,lapply(files_singletons, fread)), FUN = sum)

# combine the two datasets
d <- cbind(d1, d2[,2])
colnames(d) <- c('category','n','singletons')
d$total <- d$n + d$singletons

# make it look nice
d$singletons_pct <- round((d$singletons / d$total)*100, 2)

d_out <- data.frame(
    category = d$category,
    total = formatC(d$total, format = "f", big.mark = ",", drop0trailing = TRUE),
    singletons = paste0(d$singletons_pct, '%')
)

# Write table
rownames(d_out) <- NULL
d_out <- d_out[rev(order(d$total)),]
print(d_out)

outfile = paste0(args$out_prefix, '_.tsv')
fwrite(d_out, outfile, quote = FALSE, row.names = FALSE, sep = '\t')

## repeat but for each variant category ##
files <- list.files('derived/tables/worst_csq_for_variant_canonical/', pattern = 'ukb_wes_200k_unphased_chr[0-9]+.tsv', full.names = TRUE)
files_singletons <- list.files('derived/tables/worst_csq_for_variant_canonical/', pattern = 'ukb_wes_200k_unphased_chr[0-9]+_singletons.tsv', full.names = TRUE)

# aggregate by consequence
d1 <- aggregate(n ~ most_severe_consequence, data = do.call(rbind,lapply(files, fread)), FUN = sum)
d2 <- aggregate(n ~ most_severe_consequence, data = do.call(rbind,lapply(files_singletons, fread)), FUN = sum)

# combine the two datasets
d <- cbind(d1, d2[,2])
colnames(d) <- c('category','n','singletons')
d$total <- d$n + d$singletons
d$singletons_pct <- round((d$singletons / d$total)*100, 2)
d_out <- data.frame(
    category = d$category,
    total = formatC(d$total, format = "f", big.mark = ",", drop0trailing = TRUE),
    singletons = paste0(d$singletons_pct, '%')
)
rownames(d_out) <- NULL
d_out <- d_out[rev(order(d$total)),]
print(d_out)

outfile = paste0(args$out_prefix, '_vep_category.tsv')
fwrite(d_out, outfile, quote = FALSE, row.names = FALSE, sep = '\t')






