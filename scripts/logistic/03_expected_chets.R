#!/usr/bin/env Rscript

library(argparse)
library(data.table)
# for read_ukb_wes_kos(annotation = "pLoF_damaging_missense", chrom = 1)
source("scripts/post_hoc/utils.R")
devtools::load_all("utils/modules/R/gwastools")
devtools::load_all("utils/modules/R/prstools")

flip_main_diag <- function(mat) {
      # Transpose the matrix and reverse the columns
      flipped_mat <- t(mat)[, ncol(mat):1]
  return(flipped_mat)
}

get_trait <- function(trait, pheno_df){
    
    stopifnot(trait %in% colnames(pheno_df))
    defined <- !is.na(pheno_df[[trait]])
    subset_df <- pheno_df[defined, ]
    bool_case <- subset_df[[trait]]
    eid_case <- subset_df$eid[bool_case]
    eid_ctrl <- subset_df$eid[!bool_case]
    return(list(cases=eid_case, controls=eid_ctrl))
    
}


main <- function(args){

    path_phenotypes <- args$path_phenotypes
    out_prefix <- args$out_prefix
    ko_definition <- args$ko_definition
    variant_annotation <- args$variant_annotation

    # phenotypes we are using
    phenotypes <- gwastools::get_phenos_tested()
    
    # allow running chets+homs seperately or together
    ko_definition <- unlist(strsplit(ko_definition, split=","))
    ko_definition <- gsub("_", " ", ko_definition)
    stopifnot(all(ko_definition %in% c("Homozygote", "Compound heterozygote")))
    
    # cutoffs
    kos_cutoffs <- c(0, 5, 10, 15, 20, 25, 30, 50)

    # get phenotypes and variants
    pLoF_damaging_missense <- read_ukb_wes_kos(variant_annotation)
    pheno_df <- fread(path_phenotypes)

    final <- list()
    
    for (phenotype in phenotypes){
    
        # get case / controls
        trait <- get_trait(phenotype, pheno_df)
        final[[phenotype]] <- list()
        
        for (kos_cutoff in kos_cutoffs) {
            
            # mutations 
            mut <- mut[mut$s %in% c(trait$cases, trait$controls)]
            mut$ko <- mut$knockout %in% c("Compound heterozygote")

            # what genes can be tested (for phenotypes that are defined)
            genes_to_test <- data.table(table(mut$ko, mut$gene_id))
            genes_to_test <- genes_to_test[genes_to_test$V1 == TRUE, c(2,3)]
            colnames(genes_to_test) <- c('gene_id','kos')
            genes_to_test <- genes_to_test[genes_to_test$kos >= kos_cutoff, ]

            # subset to genes that are testable and stratify by cases and controls
            mut <- mut[mut$gene_id %in% genes_to_test$gene_id,]
            mut_cases <- mut[mut$s %in% trait$cases,]
            mut_controls <- mut[mut$s %in% trait$controls]

            # get cases that are carriers and non-carriers
            tbl_case <- data.table(table(mut_cases$ko, mut_cases$gene_id))
            tbl_case <- dcast(V2~V1, data=tbl_case, value.var="N")
            colnames(tbl_case) <- c("gene_id", "case_not_ko", "case_and_ko")

            # get controls that are carriers and non-carriers
            tbl_control <- data.table(table(mut_controls$ko, mut_controls$gene_id))
            tbl_control <- dcast(V2~V1, data=tbl_control, value.var="N")
            colnames(tbl_control) <- c("gene_id", "control_not_ko", "control_and_ko")
            fisher <- merge(tbl_case, tbl_control)
            fisher$case_or_control_and_ko <- fisher$case_and_ko + fisher$control_and_ko

            # get expected number
            fisher$n_cases <- length(trait$cases)
            fisher$n_controls <- length(trait$controls)
            fisher$p_cases <- (fisher$n_cases/(fisher$n_cases+fisher$n_controls)) 
            fisher$expected <- fisher$case_or_control_and_ko * fisher$p_cases
            fisher$kos_cutoff <- kos_cutoff
            
            # append some stats
            fisher <- merge(genes_to_test, fisher, by="gene_id")
            fisher <- cbind(phenotype, fisher)
            final[[phenotype]][[as.character(kos_cutoff)]] <- fisher
            
        }
        
    }

    # combine tables from list
    out <- rbindlist(lapply(final, rbindlist)) 
    outfile <- paste0(out_prefix, ".txt.gz")
    fwrite(out, outfile, sep="\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ko_definition", default=NULL, help = "what constitutes a 'knockout'")
parser$add_argument("--variant_annotation", default=NULL, help = "either 'pLoF' or 'pLoF_damaging_missense'")
parser$add_argument("--path_phenotypes", default=NULL, help = "path to a file containing binary phenotypes")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

