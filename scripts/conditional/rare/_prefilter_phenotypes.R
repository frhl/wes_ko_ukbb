#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(digest)

main <- function(args){

    stopifnot(file.exists(args$in_vcf))
    stopifnot(file.exists(args$pheno_file))
    stopifnot(file.exists(args$covariates))
    stopifnot(dir.exists(dirname(args$out_prefix)))

    # read in VCF
    cmd <- paste0("zcat ", args$in_vcf, " | grep -v '##' | head -n 500" )
    d <- fread(cmd = cmd)

    # get columns with genotypes / metaid
    genotype_cols <- suppressWarnings(!is.na(as.numeric(colnames(d))))
    id_cols <- suppressWarnings(is.na(as.numeric(colnames(d))))

    # Need to ensure that G is all numerics
    G <- d[,genotype_cols, with = FALSE]
    G[, names(G) := lapply(.SD, as.numeric)]

    # read in phenotypes
    phenotype <- args$phenotype
    pheno_df <- fread(args$pheno_file)
    stopifnot(phenotype %in% colnames(pheno_df))

    # are there samples with missing covariates?
    col_cov <- unlist(strsplit(readLines(args$covariates), split = ","))
    lst <- lapply(col_cov, function(col){row_ok <- is.na(pheno_df[[col]])})
    missing_cov <- rowSums(do.call(cbind, lst)) > 0
    pheno_df <- pheno_df[!missing_cov, ]
    msg <- paste("Note: Removed", sum(missing_cov),"samples with missing covariates.")
    write(msg, stdout())

    # get defined phenotypes
    defined_phenos <- !is.na(pheno_df[[phenotype]])
    eid_with_defined_phenos <- pheno_df$eid[defined_phenos]
    stopifnot(length(eid_with_defined_phenos) > 0)

    # get subset of G which contains defiend phenotypes
    G_with_defined_phenos <- colnames(G) %in% eid_with_defined_phenos
    G_subset <- G[,G_with_defined_phenos, with = FALSE]

    # Get allele count and hash for genotypes
    G_subset_AC <- rowSums(G_subset, na.rm = TRUE)
    G_subset_hash <- hash_genotypes(G_subset, "xxhash32")

    id <- d[,id_cols, with = FALSE]
    out <- data.table(
        chr = id$`#CHROM`, 
        pos = id$POS, 
        id = id$ID, 
        ref = id$REF,
        alt = id$ALT
    )

    # generate allele count outfile
    outfile_ac <- paste0(args$out_prefix,"_AC.txt.gz")
    out[[phenotype]] <- G_subset_AC
    fwrite(out, outfile_ac, sep = '\t')

    # generate allele count outfile
    outfile_hash <- paste0(args$out_prefix,"_hash.txt.gz")
    out[[phenotype]] <- G_subset_hash
    fwrite(out, outfile_hash, sep = '\t')

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--pheno_file", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--phenotype", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--covariates", default=NULL, required = TRUE, help = "Path to vcf")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Out prefix path")
args <- parser$parse_args()

main(args)









