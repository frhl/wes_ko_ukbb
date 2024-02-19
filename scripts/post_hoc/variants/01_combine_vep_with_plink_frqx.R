#!/usr/bin/env Rscript

library(argparse)
library(ggsci)
library(ggplot2)
library(reshape2)
library(data.table)


read_frqx_path <- function(frqx_path){
    dt_AC <- fread(frqx_path)
    stopifnot(nrow(dt_AC)>0)
    
    # get allele frequencies
    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, AN := AC_A1 + AC_A2]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, MAF := MAC / AN]
    dt_AC$POS <- as.numeric(sub(".*:(\\d+):.*", "\\1", dt_AC$SNP))
    
    # need to recreate SNP ID as some might have I/D in name.
    dt_AC <- dt_AC[!duplicated(dt_AC),]
    cols_to_keep <- c("CHR", "POS", "A1", "A2", "MAC", "MAF", "AN")
    dt_AC <- dt_AC[,..cols_to_keep]
    setkeyv(dt_AC, cols_to_keep[1:4])
    return(dt_AC)
}



main <- function(args){

    AC_path_before_pp_filter <- args$AC_path_before_pp_filter
    AC_path_after_pp_filter <- args$AC_path_after_pp_filter
    vep_path <- args$vep_path
    out_prefix <- args$out_prefix

    AUTOSOMES <- 1:22
    lst_after <- lapply(AUTOSOMES, function(chr){

        AC_path_before_pp_filter <- gsub("CHR", chr, AC_path_before_pp_filter)
        AC_path_after_pp_filter <- gsub("CHR", chr, AC_path_after_pp_filter)
        vep_path <- gsub("CHR", chr, vep_path)

        stopifnot(file.exists(AC_path_before_pp_filter))
        stopifnot(file.exists(AC_path_after_pp_filter))
        stopifnot(file.exists(vep_path))

        # rename to brava convention
        dt_brava_annot <- fread(vep_path)
        colnames(dt_brava_annot)[colnames(dt_brava_annot) == "varid"] <- "SNP_ID"
        colnames(dt_brava_annot)[colnames(dt_brava_annot) == "brava_csqs"] <- "annotation"
        colnames(dt_brava_annot)[colnames(dt_brava_annot) == "csqs"] <- "CSQ"

        # set keys
        setkeyv(dt_brava_annot, c("SNP_ID"))
        dt_brava_annot[, SNP_ID := gsub("chr", "", SNP_ID)]

        # read AC before and after filtering by PP
        dt_AC_before_pp_filter <- read_frqx_path(AC_path_before_pp_filter)
        dt_AC_after_pp_filter <- read_frqx_path(AC_path_after_pp_filter)
        dt_AC <- merge(dt_AC_before_pp_filter, dt_AC_after_pp_filter, suffixes = c(".before_pp", ".after_pp"), all=TRUE)
        dt_AC$MAC.after_pp[is.na(dt_AC$MAC.after_pp)] <- 0

        # assign SNP ID
        dt_AC$SNP_ID <- NA

        # need for chr:pos:a2:a1 and chr:pos:a1:a2
        SNP_ID_A2_A1 <- paste0(dt_AC$CHR, ":", dt_AC$POS, ":", dt_AC$A2, ":", dt_AC$A1)
        SNP_ID_A1_A2 <- paste0(dt_AC$CHR, ":", dt_AC$POS, ":", dt_AC$A1, ":", dt_AC$A2)
        A2_A1_in_brava_annot <- SNP_ID_A2_A1 %in% dt_brava_annot$SNP_ID
        A1_A2_in_brava_annot <- SNP_ID_A1_A2 %in% dt_brava_annot$SNP_ID

        # we set SNP ID to chr:pos:a2:a1, unless that does not match 
        # brava annotation in which case we flip it
        dt_AC$SNP_ID <- SNP_ID_A2_A1
        dt_AC$SNP_ID[!(A2_A1_in_brava_annot) & (A1_A2_in_brava_annot)] <- SNP_ID_A1_A2[!(A2_A1_in_brava_annot) & (A1_A2_in_brava_annot)]
        setkey(dt_AC, "SNP_ID")

        # now we can merge annotations
        mrg <- merge(dt_AC, dt_brava_annot, all.x=TRUE)
        return(mrg)

    })

    d <- rbindlist(lst_after)
    d$locus <- NULL
    d$alleles <- NULL
    outfile <- paste0(out_prefix, ".txt.gz")
    write(paste("writing to", outfile), stdout())
    fwrite(d, outfile, sep="\t", quote=FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--AC_path_before_pp_filter", default=NULL, help = "?")
parser$add_argument("--AC_path_after_pp_filter", default=NULL, help = "?")
parser$add_argument("--vep_path", default=NULL, help = "?")
parser$add_argument("--out_prefix", default=NULL, help = "?")
args <- parser$parse_args()

main(args)


