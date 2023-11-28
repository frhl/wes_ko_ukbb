#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

    # setup dictionary
    dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary_with_unix_codes.txt",
                         sep = "\t", header = T, stringsAsFactors = F, 
                         quote = "", comment.char = "$")

    # In this matrix, columns 6:349 are the Spiros and Duncan disease codes
    # Column names are UNIX names from the dictionary
    tte_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/all_phenos_time_to_event.txt",
                                     sep = "\t", header = T, stringsAsFactors = F,
                                     quote = "")

    # In this matrix, columns 2:308 are the Spiros disease codes
    # Columns need to be renamed to UNIX names from the dictionary
    #pheno_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_phenotype_matrix.txt", 
    #                               sep = "\t", quote = "",
    #                               header = T, stringsAsFactors = F)
    #rownames(pheno_matrix) <- pheno_matrix$eid
    #pheno_matrix <- pheno_matrix[, -1]
    #colnames(pheno_matrix) <- dictionary$unix_code

    # subset tte_matrix
    samples_to_keep <- fread(args$samples, header = FALSE)$V1
    tte_matrix <- tte_matrix[tte_matrix$eid %in% samples_to_keep,]
    
    # rename columns based on current definitions
    header <- fread(args$path_header, header = FALSE)$V1

    # remap to universal spiro naming in the project
    spiro_to <- header[grepl(pattern='spiro', header)]
    spiro_from <- gsub("spiro_","",spiro_to)
    spiro_map <- spiro_to
    names(spiro_map) <- spiro_from

    # remap to unversal palmer naming in the project
    palmer_to <- header[!grepl(pattern='spiro', header)]
    palmer_to <- palmer_to[grepl("combined", palmer_to)]
    palmer_to <- palmer_to[!grepl("primary_care", palmer_to)]
    palmer_from <- gsub("_combined(_primary_care)?","",palmer_to)
    palmer_map <- palmer_to
    names(palmer_map) <- palmer_from 

    # change spiro columns
    cnames_to_change_spiro <- colnames(tte_matrix)[colnames(tte_matrix) %in% spiro_from]
    colnames(tte_matrix)[match(cnames_to_change_spiro, colnames(tte_matrix))] <- spiro_map[cnames_to_change_spiro]

    # change palmer columns
    cnames_to_change_palmer <- colnames(tte_matrix)[colnames(tte_matrix) %in% palmer_from]
    colnames(tte_matrix)[match(cnames_to_change_palmer, colnames(tte_matrix))] <- palmer_map[cnames_to_change_palmer]

    # force all to numeric (need to be a data.frame!)
    tte_matrix <- data.frame(tte_matrix)
    cols <- which(colnames(tte_matrix) %in% header)
    tte_matrix[cols] <- lapply(tte_matrix[cols], function(x) suppressWarnings(as.numeric(x)))

    write(paste("writing to", args$out_path), stdout())
    fwrite(tte_matrix, args$out_path)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_path", default=NULL, help = "")
parser$add_argument("--path_header", default=NULL, help = "")
parser$add_argument("--samples", default=NULL, help = "")
args <- parser$parse_args()

main(args)

