
source("scripts/post_hoc/utils.R")
library(data.table)
library(argparse)
library(stringr)

main <- function(args){

    # load gtex and get top 10% specific
    d <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv")
    d <- cbind(gene_id = d$ENSGID, data.table(apply(d[,-1], 2, function(x) x > quantile(x, probs = 0.90))))
    colnames(d) <- gsub("\\_$","",gsub("\\_+","\\_",gsub("(\\()|(\\)|(\\-))","_",colnames(d))))                               

    # map from ensembl to hgnc
    bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
    ensembl_to_hgnc <- bridge$hgnc_symbol
    names(ensembl_to_hgnc) <- bridge$ensembl_gene_id

    # get columns
    cols <- colnames(d)
    tissue_specific_dt <- do.call(rbind, lapply(1:nrow(d), function(idx){
        gene <- d$gene_id[idx]
        row <- unlist(d[idx, 2:ncol(d)])
        tissue_specific <- names(row[row])
        out <- data.table(gene_id=gene, tissues=paste0(tissue_specific, collapse = ";"))
        return(out)
    }))

    # write to file
    tissue_specific_dt <- tissue_specific_dt[tissue_specific_dt$tissue != "",]
    tissue_specific_dt$hgnc_symbol <- ensembl_to_hgnc[tissue_specific_dt$gene_id]
    tissue_specific_dt <- tissue_specific_dt[,c("gene_id","hgnc_symbol","tissues")]
    outfile <- paste0(args$out_prefix, ".txt.gz")
    fwrite(tissue_specific_dt, outfile, sep = "\t")

}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Where should the results be written?")
args <- parser$parse_args()

main(args)


