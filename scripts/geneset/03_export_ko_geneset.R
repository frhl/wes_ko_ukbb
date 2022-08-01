#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(genoppi) # for msigdb geneset

null_omit <- function(lst) {
    lst[!vapply(lst, is.null, logical(1))]
}

main <- function(args){

    print(args)
    stopifnot(file.exists(args$in_vcf))
    stopifnot(file.exists(args$bridge))
    
    command <- paste0('zcat ', args$in_vcf, '| grep -v "#" | cut -f1-3')
    d <- fread(cmd=command)
    n <- nrow(d)
    colnames(d) <- c("chrom","id","ensembl_gene_id")

    # load hgnc to ensemble mapping
    bridge <- fread(args$bridge)
    bridge <- bridge[,c(1,2,5,6,7)]
    bridge <- bridge[!duplicated(bridge)]
    dt <- merge(d, bridge, all.x = TRUE)

    # setup geneset here!
    geneset <- msigdb_h_table
    colnames(geneset) <- c("hgnc_symbol","geneset")
    mrg <- merge(dt, geneset, by = "hgnc_symbol", all.x = TRUE)
    mrg <- mrg[!is.na(mrg$geneset),]
    n_mrg <- nrow(mrg)

    # combine geneset
    genesets <- unique(mrg$geneset)
    out <- lapply(genesets, function(g){
      genes <- mrg$ensembl_gene_id[mrg$geneset %in% g]
      nas <- is.na(genes)
      genes <- genes[!nas]
      n <- length(genes)
      if (n > 0){
        row1 <- paste(c(g,'var', genes), collapse = " ")
        row2 <- paste(c(g, 'anno', rep(g, n)), collapse = " ")
        return(paste0(c(row1, row2), collapse = '\n'))
      }
    })

    # write to out file
    out <- null_omit(out)
    writeLines(paste(out, collapse = '\n'), args$output_path)

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_vcf", default=NULL, help = "?")
parser$add_argument("--bridge", default=NULL, help = "?")
parser$add_argument("--output_path", default=NULL, help = "?")
args <- parser$parse_args()

main(args)

