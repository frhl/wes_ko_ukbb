
library(data.table)
library(argparse)

main <- function(args){

    ldsc <- fread(args$summary_ldsc)
    bin <- fread(args$summary_prs_bin)
    icd <- fread(args$icd)
   
    # convert to long 
    ldsc_h2 <- ldsc[ldsc$coef == "h2",]
    ldsc_intercept <- ldsc[ldsc$coef == "intercept",]
    ldsc_h2$coef <- NULL
    ldsc_intercept$coef <- NULL
    ldsc_intercept$n <- NULL
    ldsc_intercept$n_eff <- NULL
    ldsc_intercept$n_snps <- NULL

    # change columns names
    cols_to_change <- c("estimate", "std_error", "zstat", "pvalue")
    colnames(ldsc_intercept)[colnames(ldsc_intercept) %in% cols_to_change] <- paste0("intercept_",cols_to_change)
    colnames(ldsc_h2)[colnames(ldsc_h2) %in% cols_to_change] <- paste0("h2_",cols_to_change)

    # combine h2 and intercept
    ldsc_combined <- merge(ldsc_intercept, ldsc_h2, by = "phenotype")

    # combine with AUC
    ldsc_final <- merge(ldsc_combined, bin, all.x = TRUE)
    colnames(ldsc_final)[colnames(ldsc_final) == "phenotype"] <- "unix_code"

    # combine with ICD
    ldsc_auc_icd <- merge(icd[,c('phenotype','unix_code')], ldsc_final, by = "unix_code")
    stopifnot(!any(duplicated(ldsc_auc_icd$unix_code)))

    # oder by heritabiltiy p-value
    ldsc_auc_icd <- ldsc_auc_icd[order(ldsc_auc_icd$h2_pvalue),]

    # write file containing all annotations 
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(ldsc_auc_icd, outfile, sep = "\t", na = "NA")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--summary_ldsc", default=NULL, required = TRUE, help = "--")
parser$add_argument("--summary_prs_bin", default=NULL, required = TRUE, help = "--")
parser$add_argument("--icd_path", default=NULL, required = TRUE, help = "--")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


