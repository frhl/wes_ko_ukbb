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

# if cond_p col does not exist, add empty columns
conds = c('Tstat_cond','p.value_cond','varT_cond','BETA_cond','SE_cond')

if (!all(conds %in% colnames(d))){
    ndf = data.table(Tstat_cond = NA, p.value_cond = NA, varT_cond = NA, BETA_cond = NA, SE_cond = NA)
    d <- cbind(d, ndf)
    print("Adding columns")
}

# re-organize columns in table (relevant for interations>1)
d <- cbind(
    d[,-conds, with = FALSE],
    d[,conds, with = FALSE]
)

#print(head(d))

# Check rows
nrows = nrow(d)
na_values = sum(is.na(d$p.value_cond) & is.na(d$p.value))
if (na_values > 0){
    print(paste0(sum(is.na(d$p.value_cond)),"/", nrows," conditional P-values are NAs."))
    print(paste0(sum(is.na(d$p.value)),"/", nrows, " regular P-values are NAs."))
    stop('NAs found among P-values!')
}

# Select P-value
#print(head(d$p.value_cond))
#print(sum(!is.na(d$p.value_cond)))
conditioning = as.logical(sum(!is.na(d$p.value_cond)) > 0)
write(paste0('Conditioning markers included: ', conditioning),stdout())
d$Tstat_new = unlist(ifelse(conditioning, list(d$Tstat_cond), list(d$Tstat)))
d$p.value_new = unlist(ifelse(conditioning, list(d$p.value_cond), list(d$p.value)))

# save combined table
d <- d[order(d$p.value_new)]
print(head(d))
print(paste('cols:',ncol(d)))
outfile = paste0(args$prefix, '.txt')
fwrite(d, outfile, quote = FALSE, row.names = FALSE, sep = ' ', na = NA)


