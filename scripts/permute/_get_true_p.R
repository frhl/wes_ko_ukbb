#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){
 
    # input
    true_p_path <- args$true_p_path
    gene <- args$gene
    phenotype <- args$phenotype
    annotation <- args$annotation
    use_prs <- args$use_prs
    target <- args$target

    stopifnot(file.exists(true_p_path))
    d <- fread(true_p_path)
 
    # subset by phenotype and annotation
    bool_phenotype <- d$phenotype %in% phenotype
    bool_gene <- d$MarkerID %in% gene
    bool_annotation <- d$annotation %in% annotation

    # first check
    if (!any(bool_phenotype)) write(paste0("No P-values for phenotype: ", phenotype), stderr())
    if (!any(bool_phenotype)) write(paste0("No P-values for gene: ", gene), stderr())
    if (!any(bool_phenotype)) write(paste0("No P-values for Annotation: ", annotation), stderr())
    
    # subset by PRS
    d1 <- d[(bool_phenotype & bool_annotation & bool_gene),]
    bool_prs <- (!is.na(d1$prs))
    if ((use_prs == "1") & sum(bool_prs) > 0){
        d1 <- d1[bool_prs,]
        write(paste0("PRS found for ", phenotype), stderr())
    } 

    # check that we have at least one row.
    if (nrow(d1) > 1) write("Ambigious input specifications! Too many markers found.", stderr())
    if (nrow(d1) == 0) write(paste("No markers found for combination:", phenotype, gene, annotation), stderr())

    # onlt return non-NA value if fields are formatted correctly
    if (nrow(d1) == 1) {
      if (target == "p"){
          write(d1$pvalue, stdout())
      } else if (target == "t"){
          write(d1$tstat, stdout())
      } else {
          stop(paste(target, "is not a valid parameter to return!"))
      } 
    } else {
        write("NA", stdout())
    }

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--true_p_path", default=NULL, help = "A value")
parser$add_argument("--gene", default=NULL, help = "A value")
parser$add_argument("--phenotype", default=NULL, help = "A value")
parser$add_argument("--annotation", default=NULL, help = "B value")
parser$add_argument("--target", default=NULL, help = "B value")
parser$add_argument("--use_prs", default="0",  help = "operator")
args <- parser$parse_args()

main(args)

