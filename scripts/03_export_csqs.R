#!/usr/bin/env Rscript

library(argparse)
library(data.table)

main <- function(args){

    d <- fread(args$input_path)
    M <- d[,c("csqs.gene_id", "varid", "consequence_category")]
    colnames(M) <- c('gene','variant','anno')

    genes <- unique(M$gene)

    out <- do.call(rbind, lapply(genes, function(g){
        variants <- M$variant[M$gene %in% g]
        annotations <- M$anno[M$gene %in% g]
        variants[is.na(variants)] <- "N/A"
        variants[is.na(annotations)] <- "N/A"
        row1 <- c(g,'var',variants)
        row2 <- c(g, 'anno', annotations)
        stopifnot(length(row1) == length(row2))
        as.data.frame(t(data.frame(
            r1 = paste0(row1, collapse = '\t'), 
            r2 = paste0(row2, collapse = '\t'))))
    }))

    fwrite(out, args$output_path, sep = '\t', col.names = FALSE, quote = FALSE)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", default=NULL, help = "?")
parser$add_argument("--output_path", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

