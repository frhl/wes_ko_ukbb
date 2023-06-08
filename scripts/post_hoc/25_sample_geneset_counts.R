
library(ggplot2)
library(data.table)
library(argparse)

# holds all the functions to load knockouts
source("scripts/post_hoc/utils.R")

# new nice labels and order for levels
mapping <- c(
    "Heterozygote" = "Heterozygote",
    'Homozygote' = 'Homozygote',
    "Compound heterozygote" = "Compound heterozygote",
    "Compound heterozygote (cis)" = "Two-hit (Cis)",
    "Knockout" = "Knockout"
)


# a function to find the samples with at least one knockout in gene and fraction of total samples
get_sample_ko_table <- function(dt, categories, geneset, total_samples, knockout_annotation){

    stopifnot(!is.null(geneset))
    stopifnot(total_samples > 100)
    stopifnot(length(knockout_annotation) > 1)
    
    # get individual categories
    cat_lst <- lapply(categories, function(category){
        carry_ge_one_ko <- dt[dt$knockout == category,]
        carry_ge_one_geneset_ko <- dt[dt$knockout == category & dt$gene_id %in% geneset, ]
        n_all <- length(unique(carry_ge_one_ko$s))
        n_geneset <- length(unique(carry_ge_one_geneset_ko$s))
        res <- data.frame(
            category = category,
            n_all = n_all,
            n_pct_all = round( (n_all / total_samples) * 100, 2),
            n_geneset = n_geneset,
            n_pct_geneset = round( (n_geneset / total_samples) * 100, 2)
        ) 
        return(res)
    })

    cat_dt <- do.call(rbind, cat_lst)

    # go over combined category and get count
    carry_ge_one_ko <- dt[dt$knockout %in% knockout_annotation]
    carry_ge_one_geneset_ko <- dt[dt$knockout %in% knockout_annotation & dt$gene_id %in% geneset ]
    n_all <- length(unique(carry_ge_one_ko$s))
    n_geneset <- length(unique(carry_ge_one_geneset_ko$s))
    res <- data.frame(
            category = "Knockout",
            n_all = n_all,
            n_pct_all = round( (n_all / total_samples) * 100, 2),
            n_geneset = n_geneset,
            n_pct_geneset = round( (n_geneset / total_samples) * 100, 2)
    ) 

    final <- rbind(cat_dt, res)
    return(final)
}

# evaluate number of samples with at least one gene affected by a category
eval_wes_ko_geneset <- function(geneset, total_samples, exclude_ensgid = NULL){

    knockout_annotation <- c("Compound heterozygote","Homozygote")
    # read combined pLoF & damaging missense
    d <- read_ukb_wes_kos("pLoF_damaging_missense")
    if (!is.null(exclude_ensgid)) d <- d[!(d$gene_id %in% exclude_ensgid),]
    categories <- unique(d$knockout)
    d_plof_damaging_missense <- get_sample_ko_table(d, categories, geneset, total_samples, knockout_annotation)
    # read pLoF
    d <- read_ukb_wes_kos("pLoF")
    if (!is.null(exclude_ensgid)) d <- d[!(d$gene_id %in% exclude_ensgid),]
    categories <- unique(d$knockout)
    d_plof <- get_sample_ko_table(d, categories, geneset, total_samples, knockout_annotation)
    
    # read damaging missense
    d <- read_ukb_wes_kos("damaging_missense")
    if (!is.null(exclude_ensgid)) d <- d[!(d$gene_id %in% exclude_ensgid),]
    categories <- unique(d$knockout)
    d_damaging_missense <- get_sample_ko_table(d, categories, geneset, total_samples, knockout_annotation)
    
    # synonymous
    d <- read_ukb_wes_kos("synonymous")
    if (!is.null(exclude_ensgid)) d <- d[!(d$gene_id %in% exclude_ensgid),]
    categories <- unique(d$knockout)
    d_synonymous <- get_sample_ko_table(d, categories, geneset, total_samples, knockout_annotation)
    
    # other_missense
    d <- read_ukb_wes_kos("other_missense")
    if (!is.null(exclude_ensgid)) d <- d[!(d$gene_id %in% exclude_ensgid),]
    categories <- unique(d$knockout)
    d_other_missense <- get_sample_ko_table(d, categories, geneset, total_samples, knockout_annotation)
    
    # clean up names
    d_plof$category <- mapping[d_plof$category]
    d_plof_damaging_missense$category <- mapping[d_plof_damaging_missense$category]
    d_damaging_missense$category <- mapping[d_damaging_missense$category]
    d_synonymous$category <- mapping[d_synonymous$category]
    d_other_missense$category <- mapping[d_other_missense$category]
    
    # match the columns and check match
    d_plof_damaging_missense <- d_plof_damaging_missense[match(d_plof$category, d_plof_damaging_missense$category),]
    d_damaging_missense <- d_damaging_missense[match(d_plof$category, d_damaging_missense$category),]
    d_synonymous <- d_synonymous[match(d_plof$category, d_synonymous$category),]
    d_other_missense <- d_other_missense[match(d_plof$category, d_other_missense$category),]

    # check that all categories are matching    
    stopifnot(d_plof_damaging_missense$category == d_plof$category)
    stopifnot(d_plof_damaging_missense$category == d_damaging_missense$category)
    stopifnot(d_plof_damaging_missense$category == d_synonymous$category)
    stopifnot(d_plof_damaging_missense$category == d_other_missense$category)
    
    # collapse categories
    plof_damaging_missense <- collapse_categories(d_plof_damaging_missense, "pLoF_damaging_missense")
    damaging_missense <- plof <- collapse_categories(d_damaging_missense, "damaging_missense")
    plof <- collapse_categories(d_plof, "pLoF")
    synonymous <- collapse_categories(d_synonymous, "synonymous")
    other_missense <- collapse_categories(d_synonymous, "other_missense")
    
    # combine files
    final <- Reduce(merge, list(plof, damaging_missense, plof_damaging_missense, synonymous, other_missense))
    final <- final[match(mapping, final$category),]
    rownames(final) <- NULL
    return(final)
    
}



main <- function(args){

    # omim genes
    omim <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/omim//211103_morbidmap_by_gene.txt")
    omim_hgnc <- unique(omim$gene)
    omim_ensgid <- unique(omim$ensgid)    

    # protein coding genes in ensembl
    coding <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart//protein_coding_genes.tsv") 
    coding_ensgid <- unique(coding$ensembl_gene_id)

    # total samples
    total_samples = 176587

    # almost all samples have these "knockouts". We filter them out:
    common_kos_to_exclude <- c('ENSG00000094963','ENSG00000178917','ENSG00000188163')

    # evaluate omim
    dt_omim <- eval_wes_ko_geneset(omim_ensgid,  total_samples, common_kos_to_exclude)   
    outfile <- paste0(args$out_prefix,".omim.txt")
    fwrite(dt_omim, outfile, sep = "\t")

    # evalauate ensembl protein coding
    dt_omim <- eval_wes_ko_geneset(coding_ensgid,  total_samples, common_kos_to_exclude)   
    outfile <- paste0(args$out_prefix,".ensembl.protein_coding.txt")
    fwrite(dt_omim, outfile, sep = "\t")



}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


