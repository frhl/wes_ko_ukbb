
# 
devtools::load_all("utils/modules/R/gwastools")
library(data.table)
library(argparse)


phenos_with_prs <- get_phenos_prs()

# map from ENSEMBL to HGNC
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id

sig_gene_trait <- function(x, sig_T){
    x <- x[x$p.value <= sig_T,]
    x <- x[x$N_ko >= 5]
    x <- x$gene_trait
    return(x)
}

# icd codes for phenotypes
icd <- fread("data/phenotypes/phenotype_icd_chapter.txt")
icd <- icd[order(icd$ICD_chapter),]
na_chapter <- "None of the above"
icd[is.na(icd$ICD_chapter_desc_short)]$ICD_chapter_desc_short <- na_chapter 
chapters <- unique(icd$ICD_chapter_desc_short)


main <- function(args){

    print(args)
    n_tested_genes <- 952
    n_tested_phenos <- 311
    #significance_T <- 0.05 / (n_tested_genes * n_tested_phenos)
    significance_T <- 0.05 / (n_tested_genes) # nominal P

    f0 <- args$excl_prs
    f1 <- args$prefer_prs
    f2 <- args$cond_full
    f3 <- args$additive
    f4 <- args$cond_additive
    f5 <- args$common
    f6 <- args$chet_only
    f7 <- args$hom_only
    f8 <- args$only_singletons
    f9 <- args$exclude_singletons
    
    stopifnot(file.exists(f0))
    stopifnot(file.exists(f1))
    stopifnot(file.exists(f2))
    stopifnot(file.exists(f3))
    stopifnot(file.exists(f4))
    stopifnot(file.exists(f5))

    d0 <- fread(f0)
    d1 <- fread(f1)
    d2 <- fread(f2)
    d3 <- fread(f3)
    d4 <- fread(f4)
    d5 <- fread(f5)

    stopifnot(nrow(d0) > 0)
    stopifnot(nrow(d1) > 0)
    stopifnot(nrow(d2) > 0)
    stopifnot(nrow(d3) > 0)
    stopifnot(nrow(d4) > 0)
    stopifnot(nrow(d5) > 0)

    # cond hits need phenotype column to be cleaned
    d2$phenotype <- gsub("chr[0-9]+\\_","",d2$phenotype)
    d4$phenotype <- gsub("chr[0-9]+\\_","",d4$phenotype)
    d5$phenotype <- gsub("chr[0-9]+\\_","",d5$phenotype)

    stats <- d0[,c("phenotype","MarkerID","hgnc_symbol","CHR",'N_case', 'N_ctrl', 'N_ko_case', 'N_ko_ctrl', 'N_ko')]

    d0$gene_trait <- paste0(d0$phenotype,":",d0$MarkerID)
    d1$gene_trait <- paste0(d1$phenotype,":",d1$MarkerID)
    d2$gene_trait <- paste0(d2$phenotype,":",d2$MarkerID)
    d3$gene_trait <- paste0(d3$phenotype,":",d3$MarkerID)
    d4$gene_trait <- paste0(d4$phenotype,":",d4$MarkerID)
    d5$gene_trait <- paste0(d5$phenotype,":",d5$MarkerID)

    sig_gene_traits <- unique(c(
        sig_gene_trait(d0, significance_T), 
        sig_gene_trait(d1, significance_T), 
        sig_gene_trait(d2, significance_T)
     ))

    keys <- c("phenotype","MarkerID","hgnc_symbol","CHR")
    setkeyv(d0, keys)
    setkeyv(d1, keys)
    setkeyv(d2, keys)
    setkeyv(d3, keys)
    setkeyv(d4, keys)
    setkeyv(d5, keys)
    setkeyv(stats, keys) 
        
    cols_to_keep <- c("phenotype","MarkerID","hgnc_symbol","CHR","BETA","SE","p.value")
    # combine prs and no prs analysis
    d0 <- d0[,..cols_to_keep]
    d1 <- d1[,..cols_to_keep]
    d2 <- d2[,..cols_to_keep]

    # combine standard and prs analysis
    mrg <- merge(d0, d1, all=TRUE, suffixes = c(".nocond",".prs"))
    
    # combine with cond analysis
    mrg <- merge(mrg, d2, all = TRUE)
    colnames(mrg)[colnames(mrg) == "p.value"] <- "p.value.fullcond"
    colnames(mrg)[colnames(mrg) == "BETA"] <- "BETA.fullcond"
    colnames(mrg)[colnames(mrg) == "SE"] <- "SE.fullcond"
    
    # subset to hits that are significant once across all analyses
    mrg$gene_trait <-  paste0(mrg$phenotype,":",d0$MarkerID)
    mrg <- mrg[mrg$gene_trait %in% sig_gene_traits,]
    mrg <- merge(stats, mrg)
    
    # append additive analysis hits
    d3 <- d3[d3$gene_trait %in% sig_gene_traits,] 
    d3 <- d3[,..cols_to_keep]
    colnames(d3)[colnames(d3) == "p.value"] <- "p.value.additive"
    colnames(d3)[colnames(d3) == "BETA"] <- "BETA.additive"
    colnames(d3)[colnames(d3) == "SE"] <- "SE.additive"

    # aappend cond additive analysis hits
    d4 <- d4[d4$gene_trait %in% sig_gene_traits,] 
    d4 <- d4[,..cols_to_keep]
    colnames(d4)[colnames(d4) == "p.value"] <- "p.value.additivecond"
    colnames(d4)[colnames(d4) == "BETA"] <- "BETA.additivecond"
    colnames(d4)[colnames(d4) == "SE"] <- "SE.additiveccond"
   
    # aappend cond additive analysis hits
    d5 <- d5[d5$gene_trait %in% sig_gene_traits,] 
    d5 <- d5[,..cols_to_keep]
    colnames(d5)[colnames(d5) == "p.value"] <- "p.value.commoncond"
    colnames(d5)[colnames(d5) == "BETA"] <- "BETA.commoncond"
    colnames(d5)[colnames(d5) == "SE"] <- "SE.commoncond"
    
    # merge additive hits onto main table
    mrg <- merge(mrg, d5, by = keys, all.x=TRUE)
    mrg <- merge(mrg, d3, by = keys, all.x=TRUE)
    mrg <- merge(mrg, d4, by = keys, all.x=TRUE)
    
    # order by full cond
    #mrg <- mrg[order(mrg$p.value.nocond),]
    mrg <- mrg[order(mrg$p.value.fullcond),]

    # write main table
    outfile <- paste0(args$out_prefix, ".txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(mrg, outfile, sep = "\t", na="NA")

    # append hom / chet only information
    stopifnot(file.exists(f6))
    stopifnot(file.exists(f7))
 
    # get hom only analysis
    d6 <- fread(f6)
    d7 <- fread(f7) 
    
    stopifnot(nrow(d6) > 0)
    stopifnot(nrow(d7) > 0)
    
    # get counts for hom chets
    stat_chet <- d6[,c("phenotype","MarkerID","hgnc_symbol", "CHR", 'N_ko_case', 'N_ko')]
    stat_hom <- d7[,c("phenotype","MarkerID","hgnc_symbol", "CHR", 'N_ko_case', 'N_ko')]

    setkeyv(stat_hom, keys)
    setkeyv(stat_chet, keys)

    mrg_chet_hom <- merge(d6[,..cols_to_keep], d7[,..cols_to_keep], by = keys, suffixes=c(".chetonly",".homonly"))
    mrg_chet_hom_stat <- merge(stat_chet, stat_hom, by = keys, suffixes=c(".chetonly",".homonly"))

    mrg <- merge(mrg_chet_hom_stat, mrg, by=keys, all.y=TRUE)
    mrg <- merge(mrg, mrg_chet_hom, all.x=TRUE)

    # combine data with ICD mapping
    mrg <- merge(icd, mrg, by.x = "unix_code", by.y = "phenotype", all.y=TRUE)
    mrg$ICD_chapter_desc[is.na(mrg$ICD_chapter_desc)] <- na_chapter
    mrg$ICD_chapter_desc_short[is.na(mrg$ICD_chapter_desc_short)] <- na_chapter

    mrg <- mrg[order(mrg$p.value.fullcond),]
    outfile <- paste0(args$out_prefix, ".extended.txt.gz")
    write(paste0("writing ", outfile), stdout())
    fwrite(mrg, outfile, sep = "\t", na="NA")

    # extend with co_table
    a <- mrg
    d <- fread(args$co_table)
    colnames(d)[1] <- "MarkerID"
    d$gene_trait <- paste0(d$phenotype,":",d$MarkerID)
    d <- d[d$gene_trait %in% a$gene_trait,]

    # Manually populate missing values for saige traits that could not be run
    # beacuse of min_mac thresholding (for example chet_only, hom_only)
    for (i in 1:nrow(a)){
        na_row <- (is.na(a$N_ko_case.chetonly[i]) & is.na(a$N_ko_case.homonly[i]))
        if (na_row){
            pheno <- a$phenotype[i]
            gene <- a$MarkerID[i]
            row_gene_trait <- a$gene_trait[i]
            ndf <- d[d$gene_trait %in% row_gene_trait,]
            # replace values
            a$N_ko_case.chetonly[i] <- ndf$N_ko_case.chetonly
            a$N_ko_case.homonly[i] <- ndf$N_ko_case.homonly
            a$N_ko.chetonly[i] <- ndf$N_ko.chetonly
            a$N_ko.homonly[i] <- ndf$N_ko.homonly
        }
    }

    # overwrite extended file
    a <- a[order(a$p.value.fullcond),]
    outfile <- paste0(args$out_prefix, ".extended.txt.gz")
    write(paste0("Rewriting ", outfile, " with NAs replaced."), stdout())
    fwrite(a, outfile, sep = "\t", na="NA")

    #


}

# add arguments
parser <- ArgumentParser()
parser$add_argument("--excl_prs", default=NULL, required = TRUE, help = "Path to file without PRS conditioning")
parser$add_argument("--prefer_prs", default=NULL, required = TRUE, help = "Path to file with PRS conditioning when available")
parser$add_argument("--cond_full", default=NULL, required = TRUE, help = "Path to fully conditioned hits (prs, common, rare)")
parser$add_argument("--additive", default=NULL, required = FALSE, help = "Path to file with additive encoding)")
parser$add_argument("--common", default=NULL, required = FALSE, help = "Path to file with common encoding)")
parser$add_argument("--cond_additive", default=NULL, required = FALSE, help = "Path to file to recessive analysis conditioned on additive encoding")
parser$add_argument("--chet_only", default=NULL, required = FALSE, help = "Path to chet-only analysis")
parser$add_argument("--hom_only", default=NULL, required = FALSE, help = "Path to hom-only analysis")
parser$add_argument("--only_singletons", default=NULL, required = FALSE, help = "Path to hom-only analysis")
parser$add_argument("--exclude_singletons", default=NULL, required = FALSE, help = "Path to hom-only analysis")
parser$add_argument("--co_table", default=NULL, required = FALSE, help = "Path to collapsed co-occurence table to populate NAs from chet/hom only")
parser$add_argument("--out_prefix", default=NULL, required = TRUE, help = "")
args <- parser$parse_args()

main(args)


