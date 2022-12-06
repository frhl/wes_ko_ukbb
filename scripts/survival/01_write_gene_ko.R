library(argparse)
library(data.table)

format_from_hail <- function(x) gsub('(\\[)|(\\])|(")|(\\})|(\\{)','',x)

# load auxillary files
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensg_to_hgnc <- as.vector(bridge$hgnc_symbol)
names(ensg_to_hgnc) <- bridge$ensembl_gene_id
ensg_to_chrom <- as.vector(bridge$chromosome_name)
names(ensg_to_chrom) <- bridge$ensembl_gene_id

# new labels
labels <- c(
    "Heterozygote" = "het",
    "Homozygote" = "hom",
    "Compound heterozygote" = "chet_trans",
    "Compound heterozygote (cis)" = "chet_cis",
    "Possible Compound heterozygote" = "chet_unknown"
)



main <- function(args){

    stopifnot(file.exists(args$in_file))
    d <- fread(args$in_file)

    # annotate from aux files
    d$hgnc_symbol <- ensg_to_hgnc[d$gene_id]
    d$chromosome <- ensg_to_chrom[d$gene_id]

    # subset and clean
    cols_to_keep <- c("s","chromosome","gene_id","hgnc_symbol","knockout","varid")
    d <- d[,cols_to_keep, with = FALSE]
    colnames(d)[colnames(d) == "s"] <- "eid"
    colnames(d)[colnames(d) == "gene_id"] <- "ensembl_gene_id"
    colnames(d)[colnames(d) == "knockout"] <- "annotation"
    colnames(d)[colnames(d) == "varid"] <- "variants"

    # re-annotate for nicer output
    d$annotation <- labels[d$annotation]
    d$variants <- format_from_hail(d$variants)    

    # get genes
    unique_genes <- unique(d$ensembl_gene_id)
    for (gene in unique_genes){
        out_gene_prefix <- paste0(args$out_prefix, "_", gene, ".txt.gz")
        rows_to_keep <- d$ensembl_gene_id %in% gene
        d_out <- d[rows_to_keep, ]
        if (nrow(d_out) > 0){
            fwrite(d_out, out_gene_prefix, sep = "\t", quote = FALSE)
            #break
        }
    }


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_file", default=NULL, required = TRUE, help = "Directory in which to search for knockouts")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "Directory for which to output results (will create a file for each gene)")
args <- parser$parse_args()

main(args)











