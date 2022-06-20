#!/usr/bin/env Rscript

library(argparse)
library(data.table)

null_omit <- function(lst) {
    lst[!vapply(lst, is.null, logical(1))]
}

main <- function(args){

    lines <- lapply(1:22, function(chr) readLines(paste0('data/mt/vep/ukb_eur_wes_200k_csqs_chr',chr,'.tsv.gz')))
    d <- as.data.table(do.call(rbind, lines))
    #d <- fread(args$input_path)
    M <- d[,c("csqs.gene_id", "varid", "consequence_category")]
    colnames(M) <- c('gene','variant','anno')

    genes <- unique(M$gene)
    out <- lapply(genes, function(g){
      variants <- M$variant[M$gene %in% g]
      annotations <- M$anno[M$gene %in% g]
      nas <- (is.na(variants) | is.na(annotations))
      variants <- variants[!nas]
      annotations <- annotations[!nas]
      accepted <- annotations %in% c('pLoF','damaging_missense', "synonymous")
      if (sum(accepted) > 0){
          variants <- variants[accepted]
          annotations <- annotations[accepted]
          row1 <- paste(c(g,'var',variants), collapse = " ")
          row2 <- paste(c(g, 'anno', annotations), collapse = " ")
          return(paste0(c(row1, row2), collapse = '\n'))
      }
    }) 

    out <- null_omit(out)
    nout <- list()
    for (g1 in names(out)){
        for (g2 in names(out)){
            if (g1 != g2){
                gene1 <- unlist(strsplit(out[[g1]], split = "\n"))
                gene1_row1 <- unlist(strsplit(gene1[1], split = " "))
                gene1_row2 <- unlist(strsplit(gene1[2],  split = " "))
                gene1_variants <- gene1_row1[3:length(gene1_row1)]
                gene1_annotations <- gene1_row2[3:length(gene1_row2)]
                gene2 <- unlist(strsplit(out[[g2]], split = "\n"))
                gene2_row1 <- unlist(strsplit(gene2[1], split = " "))
                gene2_row2 <- unlist(strsplit(gene2[2],  split = " "))
                gene2_variants <- gene2_row1[3:length(gene2_row1)]
                gene2_annotations <- gene2_row2[3:length(gene2_row2)]
                ng <- paste0(g1,"_",g2)
                row1 <- paste(c(ng,'var', gene1_variants, gene2_variants), collapse = " ")
                row2 <- paste(c(ng, 'anno', gene1_annotations, gene2_annotations), collapse = " ")
                nout[[ng]] <- paste0(c(row1, row2), collapse = '\n')
            }
        }
    }
    writeLines(paste(nout, collapse = '\n'), args$output_path)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--output_path", default=NULL, help = "?")
parser$add_argument("--delimiter", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

