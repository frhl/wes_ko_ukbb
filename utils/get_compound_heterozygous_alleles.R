#!/usr/bin/env Rscript

# an r-script that calls a tested package that allows us to find compound heterozygoys variation. 

args = commandArgs(trailingOnly=TRUE)

# check input
if (length(args)<3) {
  stop("At least three argument must be supplied (input/input/out files)", call.=FALSE)
} 
print(args)

# load package
library(devtools)
library(data.table)
devtools::load_all('/well/lindgren/flassen/code/our/')

# setup arguments
IN_FILE1 = args[1] # VCF, phased data
IN_FILE2 = args[2] # VEP file with the variants
OUT_FILE = args[3] # where should the resulting file be written

# read data
vcf = fread(IN_FILE1, sep = '\t', skip = "#", header = T) # note, header must be true!
ref = fread(IN_FILE2, sep = ' ', header = F)

# run the script
result = get_compound_hetz_ko(vcf = vcf, vcf.col.chr = 1, vcf.col.rsid = 3,
                                 ref = ref, ref.col.chr = 1, ref.col.rsid = 2, ref.col.gene = 3, ref.col.info = 4,
                                 verbose = T)

# write out
fwrite(result, file = OUT_FILE)
write("Writing compound heterozygous alleles/samples..",stdout())


