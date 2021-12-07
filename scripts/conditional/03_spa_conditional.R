# module purge
# conda activate rpy

# setup paths and libs
library(argparse)
library(data.table)

setwd('/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/')

# add arguments
parser <- ArgumentParser()
parser$add_argument("--prefix", default=NULL, help = "Where should the results be written?")
args <- parser$parse_args()

# subet to chromosomal files
prefix = args$prefix
files = list.files(dirname(prefix), pattern = basename(prefix), full.names = TRUE)
files_chr = files[grepl('chr', files)]
print(files_chr)

# combine into single table
d <- do.call(rbind, lapply(files_chr, fread))
print(prefix)
print(head(d))
print(colnames(d))

# if cond_p col does not exist, add empty columns
conds = c('Tstat_cond','p.value_cond','varT_cond','BETA_cond','SE_cond')
bool = conds %in% colnames(d)
print(bool)

if (!all(conds %in% colnames(d))){
    ndf = data.table(Tstat_cond = NA, p.value_cond = NA, varT_cond = NA, BETA_cond = NA, SE_cond = NA)
    d <- cbind(d, ndf)
}

print(colnames(d))
bool = conds %in% colnames(d)
print(bool)
# re-organize columns in table (relevant for interations>1)
d <- cbind(
    d[,-conds, with = FALSE],
    d[,conds, with = FALSE]
)

# Select P-value
conditioning = all(is.na(d$p.value.cond)) & all(is.na(d$p.value))
write(paste0('Conditioning markers included: ', conditioning),stdout())
d$Tstat_new = unlist(ifelse(conditioning, list(d$Tstat_cond), list(d$Tstat)))
d$p.value_new = unlist(ifelse(conditioning, list(d$p.value_cond), list(d$p.value)))

# save combined table
d <- d[order(d$p.value_new)]
outfile = paste0(args$prefix, '.txt')
fwrite(d, outfile, quote = FALSE, row.names = FALSE, sep = ' ')


