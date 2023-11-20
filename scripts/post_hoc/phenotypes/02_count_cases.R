#!/usr/bin/env Rscript

source("scripts/post_hoc/utils.R")
library(argparse)
library(data.table)

labels <- c(
    "Heterozygote" = "HET",
    "Homozygote" = "HOM",
    "Compound heterozygote" = "CHET-TRANS",
    "Compound heterozygote (cis)" = "CHET-CIS",
    "Possible Compound heterozygote" = "CHET-UNKNOWN"
)

# bridge for mapping esngid -> hgnc
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensg_to_hgnc <- as.vector(bridge$hgnc_symbol)
names(ensg_to_hgnc) <- bridge$ensembl_gene_id

main <- function(args){

  stopifnot(file.exists(args$phenotypes))
  stopifnot(file.exists(args$phenos_path))
  stopifnot(file.exists(args$samples))
 
  # read in phenotypes
  phenotypes <- fread(args$phenotypes)
  
  # subset columns
  keep_cols <- c("eid","sex")
  pheno_list <- fread(args$phenos_path, header=FALSE)$V1
  phenotypes <- phenotypes[,colnames(phenotypes) %in% c(keep_cols, pheno_list),with=FALSE]
  pheno_list <- pheno_list[pheno_list %in% colnames(phenotypes)]

  # read samples to keep
  samples_to_keep <- fread(args$samples, header = FALSE)$V1
  phenotypes <- phenotypes[phenotypes$eid %in% samples_to_keep,]

  # read all knockouts
  ko <- read_ukb_wes_kos(args$annotation, chromosomes = 1:22)
  ko <- ko[,c('gene_id',"s","knockout")]
  ko$knockout <- labels[ko$knockout]
  ko <- merge(ko, phenotypes, by.x = "s", by.y = "eid", all.x = TRUE)

  by_samples <- do.call(rbind, lapply(pheno_list, function(p){
    
    # count number of samples that are knockouts
    ko_tmp <- ko[,c('s','knockout',p), with = FALSE]
    ko_tmp$p <- as.logical(ko_tmp[[p]])
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
  #by_samples$infile <- args$ko_file
  outfile_samples <- paste0(args$out_prefix, "_by_phenotypes.txt.gz")
  fwrite(by_samples, outfile_samples, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)

  # knockouts by sample
  ko_tmp <- ko[,c('s','knockout'), with = FALSE]
  ko_samples <- dcast(s ~ knockout, data = ko_tmp, fun.aggregate = length)
  
  #ko_samples$infile <- args$ko_file
  outfile_ko_samples <- paste0(args$out_prefix, "_by_samples.txt.gz")
  fwrite(ko_samples, outfile_ko_samples, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
 
  # count number of genes that are knockouts
  ko_tmp <- ko[,c('gene_id','knockout'), with = FALSE]
  ko_genes <- dcast(gene_id ~ knockout, data = ko_tmp, fun.aggregate = length)
  colnames(ko_genes)[1] <- "ensembl_gene_id"
  
  # write ko_genes
  ko_genes$hgnc_symbol <- ensg_to_hgnc[ko_genes$ensembl_gene_id]
  ko_genes$hgnc_symbol[is.na(ko_genes$hgnc_symbol) | ko_genes$hgnc_symbol == ""] <- NA
  last_col <- ncol(ko_genes)
  ko_genes <- ko_genes[,c(last_col,1:(last_col-1)), with = FALSE]

  #ko_genes$infile <- args$ko_file
  outfile_ko_genes <- paste0(args$out_prefix, "_by_genes.txt.gz")
  fwrite(ko_genes, outfile_ko_genes, append = TRUE, quote = FALSE, sep = '\t', col.names = TRUE)
 
}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotypes", default=NULL, help = "path to data.frame for phenotype status")
parser$add_argument("--phenos_path", default=NULL, help = "(path to) phenotypes to actually run")
parser$add_argument("--samples", default=NULL, help = "file (without header) containg samples to subset by")
parser$add_argument("--out_prefix", default=NULL, help = "prefix for out file")
parser$add_argument("--annotation", default="pLoF_damaging_missense", help = "what annotation to be extracted")
args <- parser$parse_args()

main(args)

