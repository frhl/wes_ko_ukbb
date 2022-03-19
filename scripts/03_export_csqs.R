#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    d <- fread(args$input_path)
    M <- d[,c("csqs.gene_id", "varid", "consequence_category")]
    colnames(M) <- c('gene','variant','anno')

    genes <- unique(M$gene)

    out <- lapply(genes, function(g){
      variants <- M$variant[M$gene %in% g]
      annotations <- M$anno[M$gene %in% g]
      nas <- (is.na(variants) | is.na(annotations))
      variants <- variants[!nas]
      annotations <- annotations[!nas]
      row1 <- paste(c(g,'var',variants), collapse = args$delimiter)
      row2 <- paste(c(g, 'anno', annotations), collapse = args$delimiter)
      return(paste0(c(row1, row2), collapse = '\n'))
    })

    writeLines(paste(out, collapse = '\n'), args$output_path)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--output_path", default=NULL, help = "?")
parser$add_argument("--delimiter", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

