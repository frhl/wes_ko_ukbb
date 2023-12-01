#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

    # get samples
    samples_qced <- fread(args$sample_file_qc, header=FALSE)$V1
    samples_vcf <- fread(args$sample_file_vcf, header=FALSE)$V1
    pedigree <- fread(args$pedigree_file)
    n_samples <- as.numeric(args$n_samples)
    set.seed(as.numeric(args$seed))

    # get samples vcf and pedigree file
    stopifnot(n_samples>1)
    stopifnot(length(samples_vcf)>0)
    stopifnot(length(samples_qced)>0)
    stopifnot(nrow(pedigree)>0)

    # ensure full trios are there
    pedigree <- pedigree[(pedigree$IID != 0)]
    pedigree <- pedigree[(pedigree$MAT != 0)]
    pedigree <- pedigree[(pedigree$PAT != 0)]
    if (nrow(pedigree)==0) stop("Pedigree has no samples left after subsetting to full trios")

    # subset pedigree to samples in VCF
    pedigree <- pedigree[(pedigree$IID %in% samples_vcf)]
    pedigree <- pedigree[(pedigree$PAT %in% samples_vcf)]
    pedigree <- pedigree[(pedigree$MAT %in% samples_vcf)]

    if (nrow(pedigree)==0) stop("Pedigree has no samples left after subsetting to VCF")

    # subset to qced samples in VCF. Note, that some (all) of the parents are not in teh QCed file.
    # we want to keep them nonetheless.
    trio_samples_in_vcf <- sort(unique(c(pedigree$IID, pedigree$PAT, pedigree$MAT)))
    final_samples_to_keep <- unique(c(samples_qced, trio_samples_in_vcf))
    samples_vcf <- samples_vcf[samples_vcf %in% final_samples_to_keep]
    stopifnot(length(samples_vcf)>0)

    # get samples not in trios and shuffle
    samples_not_in_trio <- samples_vcf[!samples_vcf %in% trio_samples_in_vcf]

    # randomly sample those that are
    if (n_samples < (length(trio_samples_in_vcf))) stop(paste("Can not include all trios at current 'n_samples'! Expecting >=", length(trio_samples_in_vcf)))
    sample_size <- n_samples - length(trio_samples_in_vcf)
    write(paste("Sampling size to use:",sample_size), stderr())
    new_samples <- sample(samples_not_in_trio, size=sample_size, replace=FALSE)

    # now append trios back in
    final_samples <- sort(c(new_samples, trio_samples_in_vcf))

    # write files
    outfile <- paste0(args$out_prefix, ".samples")
    fwrite(data.table(s=final_samples), outfile, sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

}

parser <- ArgumentParser()
parser$add_argument("--sample_file_qc", default=NULL, required = TRUE, help = "Input_file")
parser$add_argument("--sample_file_vcf", default=NULL, required = TRUE, help = "Input_file")
parser$add_argument("--pedigree_file", default=NULL, required = TRUE, help = "sample_file")
parser$add_argument("--seed", default=NULL, required = TRUE, help = "sample_file")
parser$add_argument("--n_samples", default=NULL, required = TRUE, help = "sample_file")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()


main(args)

