#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

    eid_tte_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                                 sep = "\t", header = T, stringsAsFactors = F)
    colnames(eid_tte_matrix) <- gsub("^X", "", colnames(eid_tte_matrix))

    # Disease dictionary
    dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                             sep = "\t", header = T, stringsAsFactors = F, 
                             quote = "", comment.char = "$")
    DISEASES <- dictionary$phenotype[match(colnames(eid_tte_matrix)[c(-1:-3)],
                                           dictionary$unique_code)]
    colnames(eid_tte_matrix)[c(-1:-3)] <- DISEASES
    d <- eid_tte_matrix

    d <- cbind(eid_tte_matrix$eid, eid_tte_matrix[,colnames(eid_tte_matrix) %in% DISEASES])
    colnames(d)[1] <- "eid"

    # rename samvida's phenotypes
    encoding <- fread("data/phenotypes/phenotype_icd_chapter.txt")
    stopifnot(any(colnames(d) %in% encoding$phenotype))

    # match with spiro phenotyoes
    thematch <- encoding$unix_code[match(colnames(d), encoding$phenotype)]
    not_na <- which(!is.na(thematch))
    colnames(d)[not_na] <- thematch[not_na]

    # clean up
    d <- d[,!grepl(pattern="NA",colnames(d))]
    d <- d[,-2]
    outfile <- "data/phenotypes/jan23_tte_spiro_phenotypes.txt.gz"
    fwrite(d, outfile, sep = '\t')

    # setup boolean matrix (i.e. discard tte)
    d_bool <- data.table(!is.na(d[,-1]))
    d_bool$eid <- d$eid
    d_bool <- d_bool[,c("eid",colnames(d_bool)[-ncol(d_bool)]), with =FALSE]
    outfile <- "data/phenotypes/jan23_tte_spiro_phenotypes_logical.txt.gz"
    fwrite(d_bool, outfile, sep = '\t')


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--tbd", default=NULL, help = "")
args <- parser$parse_args()

main(args)

