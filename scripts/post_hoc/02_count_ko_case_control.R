#!/usr/bin/env Rscript

library(argparse)
library(data.table)


main <- function(args){

  stopifnot(file.exists(args$phenotypes))
 
  # read in phenotypes
  phenotypes <- fread(args$phenotypes)
  binary_cols <- unlist(sapply(phenotypes, class)) == 'logical'
  binary_names <- names(binary_cols)[binary_cols]
  keep_cols <- c("eid","sex")
  phenotypes <- phenotypes[,binary_cols | colnames(phenotypes) %in%  keep_cols, with = FALSE] 
 
  # read on all knockouts
  ko <- do.call(rbind, lapply(1:22, function(chr){fread(gsub("CHR",chr,args$ko_file))}))
  ko <- ko[ko$pKO >= 0.5]
  ko <- ko[,c('gene_id',"s","knockout")]
  ko <- merge(ko, phenotypes, by.x = "s", by.y = "eid", all.x = TRUE)

  by_samples <- do.call(rbind, lapply(binary_names, function(p){
    
    # count number of samples that are knockouts
    ko_tmp <- ko[,c('s','knockout',p), with = FALSE]
    f <- as.formula(paste0("s~knockout"))
    ko_tmp <- dcast(f, data = ko_tmp, fun.aggregate = sum, value.var = p)
    ko_tmp[is.na(ko_tmp)] <- 0
    names <- colnames(ko_tmp)[-1]
    mat <- data.frame(t(matrix(colSums(ko_tmp[,-1]))))
    colnames(mat) <- names

    # annotate with all
    all_phenotypes <- phenotypes[[p]]
    mat$phenos_na <- sum(is.na(all_phenotypes))
    mat$cases_all <- sum(all_phenotypes == 1, na.rm = TRUE)
    mat$controls_all <- sum(all_phenotypes == 0, na.rm = TRUE)
    
    # annotate with non finnish europeans
    bool_nfe <- phenotypes$genetic.eur.no.fin.oct2021 == 1
    nfe_phenotypes <- phenotypes[[p]][bool_nfe]
    mat$cases_nfe <- sum(nfe_phenotypes == 1, na.rm = TRUE)
    mat$controls_nfe <- sum(nfe_phenotypes == 0, na.rm = TRUE)
    
    # clean up and move phenotype to first column
    mat$phenotype <- p
    mat <- mat[,c(ncol(mat), 1:(ncol(mat)-1))]
    return(mat)
  }))

  # write samples file
  by_samples$infile <- args$ko_file
  outfile_samples <- paste0(out_prefix, "_by_samples.txt.gz")
  fwrite(by_samples, outfile_samples, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
 
  by_genes <- do.call(rbind, lapply(binary_names, function(p){
    ko_tmp <- ko[,c('gene_id','knockout',p), with = FALSE]
    f <- as.formula(paste0("gene_id~knockout"))
    ko_tmp <- dcast(f, data = ko_tmp, fun.aggregate = sum, value.var = p)
    ko_tmp[is.na(ko_tmp)] <- 0
    names <- colnames(ko_tmp)[-1]
    mat <- data.frame(t(matrix(colSums(ko_tmp[,-1]))))
    colnames(mat) <- names

    # annotate with all
    all_phenotypes <- phenotypes[[p]]
    mat$phenos_na <- sum(is.na(all_phenotypes))
    mat$cases_all <- sum(all_phenotypes == 1, na.rm = TRUE)
    mat$controls_all <- sum(all_phenotypes == 0, na.rm = TRUE)
    
    # annotate with non finnish europeans
    bool_nfe <- phenotypes$genetic.eur.no.fin.oct2021 == 1
    nfe_phenotypes <- phenotypes[[p]][bool_nfe]
    mat$cases_nfe <- sum(nfe_phenotypes == 1, na.rm = TRUE)
    mat$controls_nfe <- sum(nfe_phenotypes == 0, na.rm = TRUE)
    
    # clean up and move phenotype to first column
    mat$phenotype <- p
    mat <- mat[,c(ncol(mat), 1:(ncol(mat)-1))]
    return(mat)
  }))
  
  # write samples file
  by_genes$infile <- args$ko_file
  outfile_genes <- paste0(out_prefix, "_by_genes.txt.gz")
  fwrite(by_genes, outfile_genes, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
 
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotypes", default=NULL, help = "chromosome")
parser$add_argument("--ko_file", default=NULL, help = "chromosome")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
args <- parser$parse_args()

main(args)

