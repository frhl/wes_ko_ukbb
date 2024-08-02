#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(argparse)

main <- function(args) {

    vep_path <- args$vep_path
    AC_path <- args$AC_path
    out <- args$out

    dt_brava_annot <- fread(vep_path)
    stopifnot(nrow(dt_brava_annot)>0)

    # re-name columns for later matching
    colnames(dt_brava_annot)[colnames(dt_brava_annot) == "varid"] <- "SNP_ID"
    colnames(dt_brava_annot)[colnames(dt_brava_annot) == "brava_csqs"] <- "annotation"
    colnames(dt_brava_annot)[colnames(dt_brava_annot) == "csqs"] <- "CSQ"
    setkeyv(dt_brava_annot, c("SNP_ID"))

    dt_AC <- fread(AC_path)
    stopifnot(nrow(dt_AC)>0)

    dt_AC[, SNP_ID := SNP]
    setkey(dt_AC, "SNP_ID")

    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, check := `C(HOM A1)` + `C(HOM A2)` + `C(HET)` + `C(MISSING)`]
    dt_AC <- dt_AC[MAC > 0]
    n_samples <- dt_AC$check[1]
    # dt_AC[, MAF := MAC/(2*n_samples)]
    dt_AC[, MAF := MAC/(AC_A1 + AC_A2)]
    dt_AC[, MAF_bin := cut(
        MAF,
        breaks=c(0, 0.001, 0.01, 0.05, 0.5, 1),
        labels=c("<0.1%", "0.1-1%", "1-5%", ">5%", ">50%"),
        include.lowest=TRUE)
    ]
    dt_AC[, MAC_bin := cut(
        MAC,
        breaks=c(0, 1, 5, 10, 100, 1000, 10000, Inf),
        labels=c("singletons", "(1,5]", "(5,10]", "(10,100]", "(100,1,000]",
            "(1,000,10,000]", ">10,000"))]
    cols <- c("SNP_ID","MAC", "MAF_bin", "MAC_bin")

    # drop old freqx columns
    cols_keep <- c("SNP_ID", "AC_A1", "AC_A2", "MAC", "check", "MAF", "MAF_bin","MAC_bin")
    dt_AC <- dt_AC[,..cols_keep]

    # merge with VEP
    dt_variant_gene_AC <- merge(dt_AC, dt_brava_annot)

    write(paste("writing to", out), stdout())
    fwrite(dt_variant_gene_AC, out, quote=FALSE, sep="\t")

}
# Add arguments
parser <- ArgumentParser()
parser$add_argument("--AC_path", default=NULL, required=TRUE,
    help="Path to allele count information, output from plink --freqx")
parser$add_argument("--vep_path", default=NULL, required=TRUE,
    help=paste0("Path to the processed VEP file"))
parser$add_argument("--out", default=NULL, required=TRUE,
    help="Output filepath")
args <- parser$parse_args()

main(args)

