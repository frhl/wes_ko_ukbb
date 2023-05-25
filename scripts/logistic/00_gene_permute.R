#!/usr/bin/env Rscript

library(argparse)
library(data.table)
# for read_ukb_wes_kos(annotation = "pLoF_damaging_missense", chrom = 1)
source("scripts/post_hoc/utils.R")

# read VCF file
fread_vcf <- function(path){
    stopifnot(file.exists(path))
    command <- paste('zcat ', path, '| grep -v "##" ')
    vcf <- fread(cmd=command, sep = '\t')
    cols <- colnames(vcf)
    final_id_col <- suppressWarnings(max(which(is.na(as.numeric(colnames(vcf))))))
    metadata <- vcf[,1:final_id_col]
    genotypes <- vcf[,(final_id_col+1):ncol(vcf)]
    return(list(metadata=metadata, genotypes=genotypes))
}

# standardize genotypes
genotypes_to_integer <- function(genotypes){
    genotypes[genotypes=="."] <- NA
    genotypes[genotypes=="0.00000"] <- 0
    genotypes[genotypes=="1.00000"] <- 1
    genotypes[genotypes=="2.00000"] <- 2
    genotypes <- sapply(genotypes, as.numeric) # matrix
    return(genotypes)
}

# make header of VCF file
make_vcf_dosage_header <- function(chrom){
    vcf_format <- '##fileformat=VCFv4.2'
    vcf_entry <-  '##FORMAT=<ID=DS,Number=1,Type=Float,Description="">'
    vcf_filter <- '##FILTER=<ID=PASS,Description="All filters passed">"'
    vcf_i1 <- '##INFO=<ID=AC,Number=A,Type=Integer,Description="Knockout count multiplied by two">'
    vcf_i2 <- '##INFO=<ID=HASH,Number=A,Type=String,Description="Hash function applied to dosages">'
    vcf_contig <- paste0('##contig=<ID=',chrom,',length=81195210>')
    vcf_out <- paste(vcf_format, vcf_entry, vcf_filter, vcf_i1, vcf_i2, vcf_contig, sep = '\n')
    return(vcf_out)
}

main <- function(args){

    autosomes <- paste0("chr",1:22)
    stopifnot(file.exists(args$input_path))
    stopifnot(args$chrom %in% autosomes)

    # seed for reproducibility
    seed <- as.numeric(args$seed)
    set.seed(seed)

    # read file and prep shuffle
    d <- fread_vcf(args$input_path)
    d$genotypes <- genotypes_to_integer(d$genotypes)
    n_snps <- length(d$metadata$POS)

    # shuffle each snp seperately
    shuffled <- list()
    for (i in 1:n_snps) {
        shuffled_snv <- sample(d$genotypes[i,])
        shuffled[[i]] <- data.table(shuffled_snv)
    }

    # combien and transpose to get correct order
    shuffled_gts <- data.table(t(do.call(cbind, shuffled)))
    colnames(shuffled_gts) <- colnames(d$genotypes)
    final <- cbind(d$metadata, shuffled_gts)
    
    # (1) write header of VCF
    vcf_out = make_vcf_dosage_header(args$chrom)
    outfile = paste0(args$out_prefix, ".vcf")
    writeLines(text = vcf_out, outfile)

    # (2) append with permuted data
    fwrite(final, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, help = "chromosome")
parser$add_argument("--input_path", default=NULL, help = "path a file containing sample ordering to be used")
parser$add_argument("--seed", default=NULL, help = "seed for randomizer")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

