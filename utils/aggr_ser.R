# script for aggregating switch errors by site


library(data.table)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--switch-error-file", default=NULL, help = "Input file generated with bcftools +trio-switch-errors")
parser$add_argument("--variant-file", default=NULL, help = "File containing the variants in the original file used for SER")
parser$add_argument("--outfile", default=NULL, help = "Destioning to output file")
args <- parser$parse_args()

file_switch_error <- args$switch_error_file
file_variants <- args$variant_file
outfile <- args$outfile

stopifnot(!is.null(file_switch_error))
stopifnot(!is.null(file_variants))
stopifnot(!is.null(outfile))

# deal with switch errors and aggregate
cmd <- paste("grep -v '^#'", file_switch_error)
dt <- fread(cmd = cmd)
dt <- dt[order(dt$Pos, dt$TRIO),]

colnames(dt) <- c('trio_id','index', 'switch', 'switches')
dt$index <- dt$index + 1 # go from 0 to 1-index

# deal with variants and combine with swtich erros
vcf <- fread(file_variants)
if (ncol(vcf) == 10) colnames(vcf) <- c('ID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'AF', 'AC', 'AN', 'HWE')
if (ncol(vcf) == 8) colnames(vcf) <- c('ID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'AF', 'HWE')
if (ncol(vcf) == 7) colnames(vcf) <- c('ID', 'CHR', 'POS', 'REF', 'ALT', 'MAF', 'AF')
if (ncol(vcf) == 6) colnames(vcf) <- c('ID', 'CHR', 'POS', 'REF', 'ALT', 'MAF')
ndt <- cbind(vcf[dt$index,], dt)

# write out
stopifnot(nrow(ndt) == nrow(dt))
fwrite(ndt, outfile, sep = '\t')

# 



