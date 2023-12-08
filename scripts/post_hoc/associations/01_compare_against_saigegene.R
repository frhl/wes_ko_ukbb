library(data.table)
library(argparse)
library(ggplot2)
library(MASS)


# map from ENSEMBL to HGNC
bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
ensembl_to_hgnc <- bridge$hgnc_symbol
names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
hgnc_to_ensembl <- bridge$ensembl_gene_id
names(hgnc_to_ensembl) <- bridge$hgnc_symbol
ensembl_to_pos <- (bridge$start_position + bridge$end_position)/2
names(ensembl_to_pos) <- bridge$ensembl_gene_id
ensembl_to_contig <- bridge$chromosome_name
names(ensembl_to_contig) <- bridge$ensembl_gene_id

group_dir_no_pp <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_saige_group_no_pp/min_mac4/"
group_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_saige_group/min_mac4"
ko_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/data/saige/output/binary/step2_encoding_012/min_mac4"

gene_traits <- fread("/well/lindgren/flassen/downloads/cell_gen_S15.csv")

gene_traits$ID <- paste0(gene_traits$`Phenotype code`,'_',gene_traits$`ensembl gene id`)
nrow(gene_traits)


mrgs <- list()
mrgs_table <- list()
traits <- readLines("data/phenotypes/dec22_phenotypes_binary_200k_header_additive_to_try.txt")
for (trait in traits){

    # assume PRS
    group_file_no_pp <- file.path(group_dir_no_pp, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense_locoprs.vep95.txt.gz"))
    group_file <- file.path(group_dir, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense_locoprs.vep95.txt.gz"))
    ko_file <- file.path(ko_dir, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense_locoprs.txt.gz"))

    # check for non-prs
    missing_group_file_no_pp <- !file.exists(group_file_no_pp)
    missing_group_file <- !file.exists(group_file)
    missing_ko_file <- !file.exists(ko_file)

    if ((missing_group_file) | (missing_ko_file) | (missing_group_file_no_pp)){
        group_file_no_pp <- file.path(group_dir_no_pp, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense.vep95.txt.gz"))
        group_file <- file.path(group_dir, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense.vep95.txt.gz"))
        ko_file <- file.path(ko_dir, paste0("ukb_eur_wes_200k_",trait,"_pLoF_damaging_missense.txt.gz"))
    }

    # check for non-prs
    missing_group_file_no_pp <- !file.exists(group_file_no_pp)
    missing_group_file <- !file.exists(group_file)
    missing_ko_file <- !file.exists(ko_file)

    if ((!missing_group_file) & (!missing_ko_file) & (!missing_group_file_no_pp)){

        # read files
        group_no_pp <- fread(group_file_no_pp)
        group <- fread(group_file)
        ko <- fread(ko_file)

        # merge files
        group$MarkerID <- group$Region
        group_no_pp$MarkerID <- group_no_pp$Region

        # subset
        cols_to_keep <- c("MarkerID", "Pvalue_Burden", "Pvalue")
        group <- group[,..cols_to_keep]
        group_no_pp <- group_no_pp[,..cols_to_keep]
        colnames(group)[2:3] <- c("Burden", "SKATO")
        colnames(group_no_pp)[2:3] <- c("Burden", "SKATO")

        # original file
        ko <- ko[,c("MarkerID","p.value")] #,"BETA","SE", "N_case_het", "N_case_hom")]
        colnames(ko)[2] <- "Pvalue_KO"

        # combine
        mrg <- merge(group, group_no_pp, by="MarkerID", all=TRUE, suffixes = c(".PP_ge_90",".PP_ge_50"))
        mrg <- merge(mrg, ko, by="MarkerID", all=TRUE)
        mrg$hgnc_symbol <- ensembl_to_hgnc[mrg$MarkerID]

        # melta nd combine
        melted <- melt(mrg, id.vars = c("MarkerID", "hgnc_symbol", "Pvalue_KO"))
        melted <- melted[!is.na(melted$Pvalue_KO),]
        colnames(melted) <- c("gene_id", "hgnc_symbol", "pvalue_ko", "variable", "pvalue_saige_gene")
        
        # create label
        melted$label <- melted$hgnc_symbol
        melted$label[melted$p.value > 1e-5] <- NA
        melted$trait <- trait

        # keep id
        melted$id <- paste0(melted$trait, "_", melted$gene_id)
        melted$to_keep <- melted$id %in% gene_traits$ID
        mrgs[[trait]] <- melted

        # geneate merge table
        mrg$trait <- trait
        mrg$id <- paste0(mrg$trait, "_", mrg$MarkerID)
        mrg$to_keep <- mrg$id %in% gene_traits$ID
        mrgs_table[[trait]] <- mrg


    } else {
        print(paste("Skipping missing trait", trait))
        if (missing_group_file) print("- missing group file")
        if (missing_ko_file) print("- missing ko file")
        if (missing_group_file_no_pp) print(paste("- missing group_no_pp file", group_file_no_pp))

    }

}

# combine and subset to relevant gene-traits
final <- rbindlist(mrgs_table)
final <- final[final$to_keep,]
final$to_keep <- NULL
length(unique(final$id))
final$id <- NULL

# do various comparisons
cols_for_min_p <- c("Burden.PP_ge_90", "SKATO.PP_ge_90", "Burden.PP_ge_50", "SKATO.PP_ge_50")
cols_for_min_p_pp90 <- c("Burden.PP_ge_90", "SKATO.PP_ge_90")
cols_for_min_p_pp50 <- c("Burden.PP_ge_50", "SKATO.PP_ge_50")
final$saige_min_p <- apply(final[,..cols_for_min_p], 1, min, na.rm=TRUE)
final$saige_min_p_p90 <- apply(final[,..cols_for_min_p_pp90], 1, min, na.rm=TRUE)
final$saige_min_p_p50 <- apply(final[,..cols_for_min_p_pp50], 1, min, na.rm=TRUE)
final$p_ko_lt_p_saigegene <- final$Pvalue_KO < final$saige_min_p
final$p_ko_lt_p_saigegene_pp90 <- final$Pvalue_KO < final$saige_min_p_p90
final$p_ko_lt_p_saigegene_pp50 <- final$Pvalue_KO < final$saige_min_p_p50
cols_from <- c('trait', 'hgnc_symbol', 'MarkerID', 'p_ko_lt_p_saigegene','saige_min_p','Pvalue_KO','Burden.PP_ge_90','SKATO.PP_ge_90','saige_min_p_p90','p_ko_lt_p_saigegene_pp90','Burden.PP_ge_50','SKATO.PP_ge_50', 'saige_min_p_p50', 'p_ko_lt_p_saigegene_pp50')
final <- final[,..cols_from]
final$saige_min_p_p90[is.infinite(final$saige_min_p_p90)] <- NA

outfile="derived/tables/231208_compare_against_saigegene.txt"
write(paste("writing to", outfile), stdout())
fwrite(final, outfile, sep="\t", na="NA")



######### plotting #########

cols_to_keep <- c("trait","hgnc_symbol", "Pvalue_KO",
                  "Burden.PP_ge_50", "SKATO.PP_ge_50",
                  "Burden.PP_ge_90", "SKATO.PP_ge_90")

d <- final[,..cols_to_keep]
melted_d <- melt(d, id.vars = c("trait", "hgnc_symbol", "Pvalue_KO"))
melted_d$pp_set <- factor(ifelse(grepl(melted_d$variable, pattern="50"), "PP>=50", "PP>=0.90"))
melted_d$model <- factor(ifelse(grepl(melted_d$variable, pattern="Burden"), "Burden", "SKATO"))

options(repr.plot.width=12, repr.plot.height=10)
p <- ggplot(melted_d, aes(x=-log10(value),y=-log10(Pvalue_KO), label=hgnc_symbol)) +
    geom_point(size=2) +
    geom_text_repel(color="black", max.overlaps = 70, point.padding = 0.15, box.padding = 1) +
    geom_abline(linetype="dashed") +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=5)) +
    xlab(expression(paste(-log[10],'(Haplotype Encoding P-value)' ))) +
    ylab(expression(paste(-log[10],'(SAIGE-GENE+ P-value)' ))) +
    labs(color="(Knockout*2) count cutoff") +
    theme_bw() +
    theme(
        axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.title.x = element_text(margin=ggplot2::margin(t=16)),
        axis.title.y = element_text(margin=ggplot2::margin(r=16)),
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.text=element_text(size=15),
        legend.position="top"
    ) +
    facet_grid(model~pp_set)
p

outfile <- file.path("/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb/derived/plots/231208_saigegene_vs_additive.pdf")
ggsave(outfile, width=12, height=10)





