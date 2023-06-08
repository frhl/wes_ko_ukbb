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
    path_unrelated <- args$path_unrelated
    kos_cutoff <- as.numeric(args$kos_cutoff)
    out_prefix <- args$out_prefix
    ko_definition <- args$ko_definition
    variant_annotation <- args$variant_annotation

    # phenotypes we are using
    phenotypes <- gwastools::get_phenos_tested()
    
    # allow running chets+homs seperately or together
    ko_definition <- unlist(strsplit(ko_definition, split=","))
    ko_definition <- gsub("_", " ", ko_definition)
    stopifnot(all(ko_definition %in% c("Homozygote", "Compound heterozygote")))

    # get phenotypes and variants
    pLoF_damaging_missense <- read_ukb_wes_kos(variant_annotation)
    pheno_df <- fread(path_phenotypes)
    print(pLoF_damaging_missense)

    final <- list()

    for (phenotype in phenotypes){

        #phenotype <- "CC_combined"
        # get case / controls
        trait <- get_trait(phenotype, pheno_df)
        write(phenotype, stderr())

        # get knockouts for unrelated samples
        unrelated <- fread(path_unrelated)
        mut <- pLoF_damaging_missense[pLoF_damaging_missense$s %in% unrelated$s,]
        mut <- mut[mut$s %in% c(trait$cases, trait$controls)]
        mut$ko <- mut$knockout %in% ko_definition

        # what genes can be tested (for phenotypes that are defined)
        genes_to_test <- data.table(table(mut$ko, mut$gene_id))
        genes_to_test <- genes_to_test[genes_to_test$V1 == TRUE, c(2,3)]
        colnames(genes_to_test) <- c('gene_id','kos')
        genes_to_test <- genes_to_test[genes_to_test$kos >= kos_cutoff, ]
        # subset to genes that are testable and stratify by cases and controls
        mut <- mut[mut$gene_id %in% genes_to_test$gene_id,]
        mut_cases <- mut[mut$s %in% trait$cases,]
        mut_controls <- mut[mut$s %in% trait$controls]

        n_cases_ko <- sum(mut_cases$ko)
        n_controls_ko <- sum(mut_controls$ko)
        if (n_cases_ko > 0) {

            # get cases that are carriers and non-carriers
            tbl_case <- data.table(table(mut_cases$ko, mut_cases$gene_id))
            tbl_case <- dcast(V2~V1, data=tbl_case, value.var="N")
            colnames(tbl_case) <- c("gene_id", "case_not_ko", "case_and_ko")
            print(head(tbl_case))

            # get controls that are carriers and non-carriers
            tbl_control <- data.table(table(mut_controls$ko, mut_controls$gene_id))
            tbl_control <- dcast(V2~V1, data=tbl_control, value.var="N")
            colnames(tbl_control) <- c("gene_id", "control_not_ko", "control_and_ko")
            fisher <- merge(tbl_case, tbl_control)

            # run fisher's exact test
            exact_by_trait <- rbindlist(lapply(1:nrow(fisher), function(row_idx){
                row <- fisher[row_idx, ]
                # we need to flip the matrix to get the OR for increased risk
                fisher_matrix <- flip_main_diag(matrix(as.numeric(row[,-1]), nrow=2, ncol=2))
                res <- fisher.test(fisher_matrix)
                row$p.value <- res$p.value
                row$conf.int.lower <- res$conf.int[1]
                row$conf.int.upper <- res$conf.int[2]
                row$estimate <- res$estimate
                row$ko_cutoff <- kos_cutoff
                row$null.value <- res$null.value
                row$alternative <- res$alternative
                row$pop.cases <- length(trait$cases)
                row$pop.controls <- length(trait$controls)
                return(row)
            }))

            # order by P-value
            exact_by_trait <- merge(genes_to_test, exact_by_trait, by ="gene_id")
            exact_by_trait <- exact_by_trait[order(exact_by_trait$p.value)]
            exact_by_trait <- cbind(phenotype, exact_by_trait)
            final[[phenotype]] <- setDT(exact_by_trait)
        }
    }

    out <- rbindlist(final)
    outfile <- paste0(out_prefix, ".txt.gz")
    fwrite(out, outfile, sep="\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--ko_definition", default=NULL, help = "what constitutes a 'knockout'")
parser$add_argument("--kos_cutoff", default=NULL, help = "How many occurences in the population before we test?")
parser$add_argument("--variant_annotation", default=NULL, help = "either 'pLoF' or 'pLoF_damaging_missense'")
parser$add_argument("--path_unrelated", default=NULL, help = "path to a file containing unrelated samples")
parser$add_argument("--path_phenotypes", default=NULL, help = "path to a file containing binary phenotypes")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

