
library(data.table)
library(pscl)


# load all consequences
plof_csqs = c("transcript_ablation", "splice_acceptor_variant",
             "splice_donor_variant", "stop_gained", "frameshift_variant")

missense_csqs = c("stop_lost", "start_lost", "transcript_amplification",
                 "inframe_insertion", "inframe_deletion", "missense_variant")

synonymous_csqs = c("stop_retained_variant", "synonymous_variant")

other_csqs = c("mature_miRNA_variant", "5_prime_UTR_variant",
              "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
              "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant",
              "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
              "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
              "regulatory_region_variant", "feature_truncation", "intergenic_variant")



# get path to full file of ukb wes knockouts
ukb_wes_ko_path <- function(annotation = "pLoF_damaging_missense", chr = "21"){
    thedir <- paste0("data/knockouts/alt/pp90/recoded")
    thefile <- paste0("ukb_eur_wes_200k_chr",chr,".pp90.recoded.",annotation,".txt.gz")
    path <- file.path(thedir, thefile)
    if (!file.exists(path)) stop(paste(path, "does not exist!"))
    return(path)
}

# get path to full file of ukb wes knockouts
ukb_wes_syn_path <- function(annotation = "synonymous", chr = "21"){
    #thedir <- paste0("data/knockouts/alt/pp90/recoded_tmp")
    thedir <- paste0("data/knockouts/alt/pp90/recoded_other_by_merge")
    thefile <- paste0("ukb_eur_wes_200k_chr",chr,".pp90.recoded.",annotation,".txt.gz")
    path <- file.path(thedir, thefile)
    if (!file.exists(path)) stop(paste(path, "does not exist!"))
    return(path)
}

# read full paths
read_ukb_wes_kos <- function(annotation, chromosomes=1:22, allow_hets = TRUE){
    d <- do.call(rbind, lapply(chromosomes, function(chr){
        if (annotation %in% c("pLoF", "damaging_missense", "pLoF_damaging_missense")){
            path <- ukb_wes_ko_path(annotation, chr)
            d <- fread(path)
            
        } else if (annotation %in% c("synonymous", "other_missense")){
            path <- ukb_wes_syn_path(annotation, chr)
            d <- fread(path)
        }
    }))
    if (!allow_hets) d <- d[!d$knockout %in% "Heterozygote", ]
    return(d)
}


get_mapping_hgnc_to_ensembl <- function(){
    bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
    hgnc_to_ensembl <- bridge$ensembl_gene_id
    names(hgnc_to_ensembl) <- bridge$hgnc_symbol
    return(hgnc_to_ensembl)
}


get_mapping_ensembl_to_hgnc <- function(){
    bridge <- fread("/well/lindgren/flassen/ressources/genesets/genesets/data/biomart/220524_hgnc_ensg_enst_chr_pos.txt.gz")
    ensembl_to_hgnc <- bridge$hgnc_symbol
    names(ensembl_to_hgnc) <- bridge$ensembl_gene_id
    return(ensembl_to_hgnc)
}


# collapse categories into a single column
collapse_categories <- function(dt, annotation){
    dt_out <- data.frame(
        category = dt$category,
        n_all = paste0(dt$n_all, " (", dt$n_pct_all,"%)"),
        n_geneset = paste0(dt$n_geneset, " (", dt$n_pct_geneset,"%)")
    )
    colnames(dt_out)[2:3] <- paste0(annotation,"_",colnames(dt_out)[2:3])
    return(dt_out)
}



