#!/usr/bin/env Rscript

library(argparse)
library(data.table)

bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_hgnc <- bridge$hgnc_symbol
names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id

# calcualte empircal P-value
get_emp_p <- function(d){
    true_p <- d$p[d$is_permuted == 0]
    pvalue <- d$p[(d$is_permuted == 1) & (grepl("ENSG", d$out_marker))]
    #emp_p <- sum(true_p > pvalue)/length(pvalue)
    emp_p <- sum(true_p >= pvalue)/length(pvalue)
    return(emp_p)
}

main <- function(args){
    permute_dir <- args$permute_dir
    path_trait_genes <- args$path_trait_genes
    stopifnot(dir.exists(permute_dir))
    stopifnot(file.exists(path_trait_genes))
    trait_gene_combos <- fread(path_trait_genes)
    phenotypes <- unique(trait_gene_combos$trait)
    stopifnot(length(phenotypes)>1)
    for (trait in phenotypes){
        pattern <- paste0(trait,".pvalues")
        files <- list.files(permute_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
        # aggregate permuted and observed pvalues into a single file
        d <- rbindlist(lapply(files, function(f){
            d <- fread(f)
            gene <- stringr::str_extract(basename(f), pattern = "ENSG[0-9]+")
            d$hgnc_symbol <- ensembl_to_hgnc[gene]
            d$ensembl_gene_id <- d$marker_out
            d$marker_out <- NULL
            return(d)
        }))
        out_prefix_pheno <- paste0(args$out_prefix,".permuted.p.",trait,".txt.gz")
        fwrite(d, out_prefix_pheno, sep = '\t')
        
        # calculate empirical P-value 
        d <- rbindlist(lapply(files, function(f){
            d <- fread(f)
            gene <- stringr::str_extract(basename(f), pattern = "ENSG[0-9]+")
            hgnc <- ensembl_to_hgnc[gene]
            permutations <- sum(d$is_permuted)
            return(data.table(gene = gene, hgnc_symbol=hgnc, p.obs=get_emp_p(d), permutations=permutations))
        }))
        # write the resulting file
        out_prefix_pheno <- paste0(args$out_prefix,".empirical.p.",trait,".txt.gz")
        fwrite(d, out_prefix_pheno, sep = '\t')
    }    
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--permute_dir", default=NULL, required = TRUE, help = "")
parser$add_argument("--path_trait_genes", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)

