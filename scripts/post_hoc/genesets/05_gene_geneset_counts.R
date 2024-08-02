
library(ggplot2)
library(data.table)
library(argparse)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")
hgnc_to_ensembl <- get_mapping_hgnc_to_ensembl()

# new nice labels and order for levels
mapping <- c(
    "Heterozygote" = "Heterozygote",
    'Homozygote' = 'Homozygote',
    "Compound heterozygote" = "Compound heterozygote",
    "Compound heterozygote (cis)" = "Two-hit (Cis)",
    "Knockout" = "Knockout"
)

# go over a geneset and count
check_geneset <- function(dt, geneset, to_real_table = TRUE){
    df <- dt
    df$geneset <- df$gene_id %in% geneset
    gene_df <- aggregate(gene_id ~ geneset + knockout, data = df, FUN = length)
    gene_df <- gene_df[gene_df$geneset == TRUE,]
    colnames(gene_df)[3] <- "n"
    gene_df$geneset_size <- length(unique(geneset))
    gene_df$pct_of_geneset <- round((gene_df$n / gene_df$geneset_size) * 100,1)
    gene_df$geneset <- NULL
    rownames(gene_df) <- NULL
    if (to_real_table){
        gene_df$geneset_size <- NULL
        gene_df$n <- paste0(gene_df$n, ' (', gene_df$pct_of_geneset, '%)')
        gene_df$pct_of_geneset <- NULL
    }
    return(gene_df)
}


main <- function(args){

    # get knockouts 
    annotation <- args$annotation
    dt <- read_ukb_wes_kos(annotation)
    dt <- dt[,c("gene_id", "knockout",'pKO')]
    dt <- dt[!duplicated(dt),]
    ndt <- dt[(dt$pKO == 1)]
    ndt$knockout <- "Knockout"
    ndt <- ndt[!duplicated(ndt$gene_id),]
    dt <- rbind(dt, ndt)

    # protein coding
    coding <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart//protein_coding_genes.tsv") 

    # omim genes
    omim <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/omim//211103_morbidmap_by_gene.txt")
    omim_hgnc <- unique(omim$gene)
    omim_ensgid <- unique(omim$ensgid)

    # gtex brain
    gtex <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gtex/GTEx.tstat.tsv")
    gtex_brain <- gtex$ENSGID[gtex$Brain_Cortex >= quantile(gtex$Brain_Cortex, probs = 0.9)]

    # gnomad
    gnomad_genes <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/gnomad/karczewski2020/supplementary_dataset_13_gene_lists.tsv")
    gnomad_genes$gene_id <- hgnc_to_ensembl[gnomad_genes$gene]
    genes_auto_dom <- gnomad_genes$gene_id[gnomad_genes$gene_list %in% "Autosomal Dominant"]
    genes_auto_rec <- gnomad_genes$gene_id[gnomad_genes$gene_list %in% "Autosomal Recessive"]
    genes_olf <- gnomad_genes$gene_id[gnomad_genes$gene_list %in% "Olfactory Genes"]
    genes_hap <- gnomad_genes$gene_id[gnomad_genes$gene_list %in% "Haploinsufficient"]

     
    genes <- check_geneset(dt, unique(dt$gene_id))
    colnames(genes)[2] <- "Genes n (% of geneset)"
    protein_coding <- check_geneset(dt, unique(coding$ensembl_gene_id))
    colnames(protein_coding)[2] <- "Protein Coding n (% of geneset)"
    omim <- check_geneset(dt, omim_ensgid)
    colnames(omim)[2] <- "OMIM n (% of geneset)"
    dom <- check_geneset(dt, genes_auto_dom)
    colnames(dom)[2] <- "Autosomal Dominant n (% of geneset)"
    rec <- check_geneset(dt, genes_auto_rec)
    colnames(rec)[2] <- "Autosomal Recessive n (% of geneset)"
    olf <- check_geneset(dt, genes_olf)
    colnames(olf)[2] <- "Olfactory Genes n (% of geneset)"
    hap <- check_geneset(dt, genes_hap)
    colnames(hap)[2] <- "Haploinsufficient n (% of geneset)"

    final <- Reduce(merge, list(genes, protein_coding, omim, dom, rec, olf))
    final$knockout <- mapping[final$knockout]
    final <- final[c(3,2,1,4,5),]
    rownames(final) <- NULL

    # get gtex categories
    outfile <- paste0(args$out_prefix,".", annotation, ".txt")
    fwrite(final, outfile, sep = "\t")

}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
parser$add_argument("--annotation", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


