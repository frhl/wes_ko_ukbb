
# 
devtools::load_all("utils/modules/R/gwastools")
library(data.table)
library(argparse)

# map from ENSEMBL to HGNC
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id

sig_gene_trait <- function(x, sig_T){
    x <- x[x$p.value <= sig_T,]
    x <- x[x$N_ko_case >= 2]
    x <- x[x$N_ko >= 5]
    x <- x$gene_trait
    return(x)
}

main <- function(args){

    n_tested_genes <- 958
    n_tested_phenos <- 311
    significance_T <- 0.05 / (n_tested_genes * n_tested_phenos)

    f0 <- args$excl_prs
    f1 <- args$prefer_prs
    f2 <- args$cond_full
    f3 <- args$chet_only
    f4 <- args$hom_only

    stopifnot(file.exists(f0))
    stopifnot(file.exists(f1))
    stopifnot(file.exists(f2))

    d0 <- fread(f0)
    d1 <- fread(f1)
    d2 <- fread(f2)
    d2$phenotype <- gsub("chr[0-9]+\\_","",d2$phenotype)

    stats <- d0[,c("phenotype","MarkerID","hgnc_symbol","CHR",'N_case', 'N_ctrl', 'N_ko_case', 'N_ko_ctrl', 'N_ko')]

    d0$gene_trait <- paste0(d0$phenotype,":",d0$MarkerID)
    d1$gene_trait <- paste0(d1$phenotype,":",d1$MarkerID)
    d2$gene_trait <- gsub("chr[0-9]+_","",paste0(d2$phenotype,":",d2$MarkerID))

    sig_gene_traits <- unique(c(
        sig_gene_trait(d0, significance_T), 
        sig_gene_trait(d1, significance_T), 
        sig_gene_trait(d2, significance_T)
     ))

    keys <- c("phenotype","MarkerID","hgnc_symbol","CHR")
    setkeyv(d0, keys)
    setkeyv(d1, keys)
    setkeyv(d2, keys)
    setkeyv(stats, keys) 
        
    cols_to_keep <- c("phenotype","MarkerID","hgnc_symbol","CHR","BETA","SE","p.value")
    # combine prs and no prs analysis
    mrg <- merge(d0[,..cols_to_keep], d1[,..cols_to_keep], all=TRUE, suffixes = c(".nocond",".prs"))
    # combine with cond analysis
    mrg <- merge(mrg, d2[,..cols_to_keep], all = TRUE)
    colnames(mrg)[colnames(mrg) == "p.value"] <- "p.value.fullcond"
    colnames(mrg)[colnames(mrg) == "BETA"] <- "BETA.fullcond"
    colnames(mrg)[colnames(mrg) == "SE"] <- "SE.fullcond"
    # subset to hits that are significant once across all analyses
    mrg$gene_trait <-  paste0(mrg$phenotype,":",d0$MarkerID)
    mrg <- mrg[mrg$gene_trait %in% sig_gene_traits,]
    mrg <- merge(stats, mrg)
    mrg <- mrg[order(mrg$p.value.nocond),]

    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(mrg, outfile, sep = "\t")

    # append hom / chet only information
    stopifnot(file.exists(f3))
    stopifnot(file.exists(f4))
 
    # get hom only analysis
    d3 <- fread(f3)
    d4 <- fread(f4) 
    
    # get counts for hom chets
    stat_chet <- d3[,c("phenotype","MarkerID","hgnc_symbol", "CHR", 'N_ko_case', 'N_ko')]
    stat_hom <- d4[,c("phenotype","MarkerID","hgnc_symbol", "CHR", 'N_ko_case', 'N_ko')]

    setkeyv(stat_hom, keys)
    setkeyv(stat_chet, keys)

    mrg_chet_hom <- merge(d3[,..cols_to_keep], d4[,..cols_to_keep], by = keys, suffixes=c(".chetonly",".homonly"))
    mrg_chet_hom_stat <- merge(stat_chet, stat_hom, by = keys, suffixes=c(".chetonly",".homonly"))

    mrg <- merge(mrg_chet_hom_stat, mrg, by=keys, all.y=TRUE)
    mrg <- merge(mrg, mrg_chet_hom, all.x=TRUE)

    outfile <- paste0(args$out_prefix, ".extended.txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(mrg, outfile, sep = "\t")


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--excl_prs", default=NULL, required = TRUE, help = "")
parser$add_argument("--prefer_prs", default=NULL, required = TRUE, help = "")
parser$add_argument("--cond_full", default=NULL, required = TRUE, help = "")
parser$add_argument("--chet_only", default=NULL, required = TRUE, help = "")
parser$add_argument("--hom_only", default=NULL, required = TRUE, help = "")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


