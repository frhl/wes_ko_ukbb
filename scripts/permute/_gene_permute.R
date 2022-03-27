#!/usr/bin/env Rscript

library(argparse)
library(data.table)

# simple method to shuffle knockouts
shuffle_knockouts <- function(d){
    d$KO <- rbinom(n=nrow(d), size=1, prob = d$pTKO)
    d$pKO <- ifelse((d$KO == 1), 1,
             ifelse((d$phased_het == 1 & d$unphased_het > 0), 1 - 1*(1/2)^d$unphased_het,
             ifelse((d$phased_het == 0 & d$unphased_het > 1), 1 - 2*(1/2)^d$unphased_het, 0)))
    return(d$pKO)
}

main <- function(args){
    
    print(args)
    stopifnot(file.exists(args$input_path))
    stopifnot(!is.na(as.numeric(args$permutations)))

    # seed for reproducibility
    seed <- as.numeric(args$seed)
    set.seed(seed)
    
    # replicate knockout 
    n <- as.numeric(args$permutations)
    d <- fread(args$input_path)
    reps <- replicate(n, shuffle_knockouts(d))
    rownames(reps) <- d$s
    reps <- data.table(t(reps))
    
    # knockout ceiling
    #ko_count <- ceiling(apply(reps, 1, sum))
    
    # synthethic row
    row <- data.table(
      "#CHROM" = args$chrom, 
      POS = 1:n,
      ID = args$vcf_id,
      REF = '0',
      ALT = '1',
      INFO = '.',
      FORMAT = 'DS'
    )

    # combine synthethic rows with knockout matrix
    M <- cbind(row, reps)

    # create synthethic VCF output
    vcf_format <- '##fileformat=VCFv4.2'
    vcf_entry <-  '##FORMAT=<ID=DS,Number=1,Type=Float,Description="">'
    vcf_filter <- '##FILTER=<ID=PASS,Description="All filters passed">"'
    vcf_contig <- paste0('##contig=<ID=',args$chr,',length=81195210>')
    vcf_out <- paste(vcf_format, vcf_entry, vcf_filter, vcf_contig, sep = '\n')

    # write header of vcf
    outfile = paste0(args$out_prefix, ".vcf")
    writeLines(text = vcf_out, outfile)
    fwrite(M, outfile, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--chrom", default=NULL, help = "chromosome")
parser$add_argument("--input_path", default=NULL, help = "path to the input")
parser$add_argument("--permutations", default=NULL, help = "number of times the gene should be permuted")
parser$add_argument("--seed", default=NULL, help = "seed for randomizer")
parser$add_argument("--vcf_id", default="GENE", help = "Substitute for rsid")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

