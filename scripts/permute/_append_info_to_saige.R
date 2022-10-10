#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){
    
    stopifnot(file.exists(args$input_vcf))
    stopifnot(file.exists(args$input_saige))

    incmd <- paste("zcat ", args$input_vcf, " | grep -v '##'")
    d <- fread(cmd=incmd)
    # subset to columns with meta information
    d <- d[,is.na(suppressWarnings(as.numeric(colnames(d)))), with = FALSE]
    # make artificial merge header
    d$mrg_id <- paste0(d$POS,"_",d$ID,"_",d$REF,"_",d$ALT)
    
    # setup info
    INFO <- gsub("(AC\\=)|(HASH\\=)","",d$INFO)
    d_info <- data.table(do.call(rbind, strsplit(INFO, split = ';')))
    colnames(d_info) <- c("AC_VCF","HASH_VCF")
    d_info$mrg_id <- d$mrg_id

    # read saige
    saige <- fread(args$input_saige)
    saige$mrg_id <- paste0(saige$POS,"_",saige$MarkerID,"_",saige$Allele1,"_",saige$Allele2)

    #saige
    mrg <- merge(saige, d_info, all.x = TRUE, by = "mrg_id")
    mrg$mrg_id <- NULL
    
    # write to output
    fwrite(mrg, args$output, quote = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_saige", default=NULL, help = "path to the input")
parser$add_argument("--input_vcf", default=NULL, help = "path to the input")
parser$add_argument("--output", default=1, help = "path to the input")
args <- parser$parse_args()

main(args)

